#!/usr/bin/env Rscript

# This script will load netcdfs and create two different final datasets:
# Dataset 1 contains the highres aridity index of presupplied data (Prec/PotEvap), EVI and EVI-quality. All at 0.05 degree, which gives 22475 cells that are indexed.
# Dataset 2 contains the modeled variables, including an aridity variable (Prec - PotEvap) as each model estimates it's own PotEvap. Some extra corrections take place because the sign-conventions are not strictly followed on the dataportal. The averaged EVI is also present. All at 0.25 degree, which gives 899 cells that are indexed.
# Two other DEM datasets. One at the highest available resolution (to supply elevation values to the groundwater observations when they are lacking.)

# Notes:
# The workflow depends on a hierarchy of functions. One function opens netcdfs, reads timestamps and locations so it can be converted to a tidy format. Then it calls upon other functions for further operations on the values. In case of Dataset 2 the whole operation of variable reading and transformation needs to be done multiple times (to produce raw, normalized, and empirical cdf values) so there is a meta function above it all.
# Inherits the working directory from the bash script
# Whether the addition of cellnrs is consistent over the model data depends on the data storage in the netcdf file (inversion of latitudes, which is corrected in the download_e2o_script.sh)
# Normalisation happens per cell and per month, this means that there is a maximum of 14 values to compute mean and stdev with. Additionally the assumption of a normal distribution to derive probabilities does not seem appropriate for all variables. Then defenition of drought in terms of percentiles becomes problematic. For visual inspection it is still nice though.
# Therefore it is better not to assume a parametric shape and transform values to the associated probabilities 

require(raster)
require(tidyr)
require(tibble)
require(lubridate)
require(dplyr)
require(ggplot2)

# ===========================
# Part 1 function definitions
# ===========================

# High level function that loads a 3d (lon,lat,time) netcdf set, and extracts one variable + timestamps + (optionally) lat/lon/cellnrs. 
# It also has the option to apply an intermediate function to the values (such as normalisation or the conversion to Quality bitsequence)
# It calls upon the mygather function below for a different evaluation of the arguments.
extract_var <- function(filepath = "./stack_p05deg.nc", variable = "evi", timeandloc = TRUE, FUN=NULL, add.args=NULL) {
  
  # The stack() function from raster can correctly read the netcdfs because we ensured uniform time-axes and orientation. It recognizes the Fill values as stored in the variable attributes and places NA:
  capture <- as.data.frame(stack(x = filepath, varname = variable), xy=timeandloc)
  
  if (timeandloc) {
    capture$cellnr <- 1:nrow(capture) # Cellnumber is added for reference as a last entry (starts north-west/top-left, increases row wise)
    melt <- mygather(capture, "timestamp", variable, 3:(ncol(capture)-1)) # This gathering creates that the first ncell rows have the first timestamp.
  } else {
    melt <- mygather(capture, "timestamp", variable, 1:ncol(capture))
  }
  
  # To convert the timestamp string to something useful, and transfer the whole thing to tibble
  melt$time <- as_date(substring(melt$timestamp, first = 2, last = 11))
  result <- as_tibble(melt[,which(names(melt) != "timestamp")])
  
  # Call upon the intermediate function (which is supplied the whole tibble at once and should output a vector)
  if (is.null(FUN)) {
    return(result)
  } else {
    result[[eval(variable)]] <- do.call(what = FUN, args = c(list(result), add.args))
    return(result)
  }
}

# Adapted gather function to work with column names contained in variables. As the original gather function cannot be programmed with. Used inside extract_var() and create_set2()
mygather <- function(mydata, key.col, val.col, gather.cols) {
  new.data <- gather_(data = mydata,
                      key_col = key.col,
                      value_col = val.col,
                      gather_cols = colnames(mydata)[gather.cols])
  return(new.data)    
}

# Intermediate function for conversion of the MODIS quality variable: in the netCDF it is stored as an integer (see download_MODIS_script.sh) and to read it it must become a Quality bitsequence
# The meaning of the bits can be found in https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod13_user_guide.pdf table 5 and beyond
qualfun <- function(frame, name = "quality", noBits = 16) {
  # Function that converts a single integer into a character bitsequence:
  number2binary = function(number, noBits) {
    binary_vector = as.integer(intToBits(number))
    paste(binary_vector[1:noBits], collapse = "")
  }
  # apply it to the whole vector (vapply is fastest)
  return(vapply(X = frame[[name]], FUN = number2binary, FUN.VALUE = "character", USE.NAMES = FALSE , noBits = noBits))
}

# Intermediate function for evi scaling.
evifun <- function(frame, name = "evi", fctor = 0.0001) {
  return(frame[[name]]*fctor)
}

# Intermediate function for monthly normalisation of a timeseries, uses the 'time' column and can be provided the name of the variable to be normalized. Separates gridcells in the mean and sd computation. Has an option for sign reversion. Because there is a maximum of 14 cells, which is not much for the computation of the centralised moments, all should be present. So no na.rm = TRUE in the scaling.
normalisation <- function(frame, name = "PotEvap", reverse = FALSE) {
  if (reverse) {
    frame[[name]] <- -frame[[name]]
  }
  # Grouping by months and cellnr. Compute mean, stdev for each group and normalise (all in scale() function)
  frame <- frame %>% mutate(month = month(time)) # month() is the function from lubridate, not the name
  frame %>% group_by(month, cellnr) %>% mutate_at(funs(scale(.) %>% as.vector), .vars = c(name)) %>% pull(name)
}

# Intermediate function for transformation of values (per month and per cell) to an estimated empirical non-exceedence probability (e.g. when p[X <= x ] is very low this x value can be characterised as drought). Calls upon low-level function below. Has an option for sign reversion.
probabilisation <- function(frame, name = "PotEvap", reverse = FALSE) {
  if (reverse) {
    frame[[name]] <- -frame[[name]]
  }
  # Grouping by months and cellnr.
  frame <- frame %>% mutate(month = month(time)) # second month is the function from lubridate
  frame %>% group_by(month, cellnr) %>% mutate_at(funs(cdf_est(., a = 0.5, na.tolerate = FALSE)), .vars = c(name)) %>% pull(name) 
}

# Low level function that takes in a vector and estimates the cumulative probability for all values in this set (non-exceedence). It uses an estimator of choice instead of the native ecdf() of R, because the latter does not permit occurences outside the supplied range.
# Just like the normalisation, also here the number of observations can be low. So it has the option to tolerate the absence of values, or to return nothing for the whole vector.
cdf_est <- function(vals, a = 0.5, na.tolerate = FALSE) {
  if ((na.tolerate & any(is.na(vals))) | (!any(is.na(vals)))) {
    ranks <- rank(x = vals, na.last = "keep", ties.method = "random") # Random resolvement of ties is a choice, however not likely to occur with values from numerical simulation with a lot of decimals.
    prob <- (ranks - a) / (sum(!is.na(ranks)) + 1 - 2*a) # P = (r-a)/(N+1-2a), empirical estimation function. if a = 1 then it matches the ecdf() function of R.
    return(prob)
  } else {
    return(rep(as.numeric(NA), length.out = length(vals)))
  }
}

# Intermediate function that works allows for a similar datastream (tibble to numeric) as the other intermediate function, but its goal is to do nothing except for a possible reversion
notransformation <- function(frame, name = "PotEvap", reverse = FALSE) {
  if (reverse) {
    frame[[name]] <- -frame[[name]]
  }
  frame %>% pull(name)
}

# Intermediate function to compute the average value per cell over the whole timeseries.
tempavg <- function(frame, name = "PotEvap") {
  frame %>% group_by(cellnr) %>% mutate_at(funs(mean(.)), .vars = c(name)) %>% pull(name)
}

# ===========================
# Part 2 Build dataset 1
# ===========================

# Building of the first dataset. Combines evi, quality, and the highres aridity index at 0.05 degree resolution.
print("START build of dataset 1 at 0.05 degree")

# EVI is only scaled to its original range [-0.2, 1] and not normalised. For quality, the bit conversion qual_fun() is called.
evi_p05deg <- extract_var(filepath = "stack_p05deg.nc", variable = "evi", timeandloc = TRUE, FUN = evifun, add.args = list(name = "evi"))
qual_p05deg <- extract_var(filepath = "stack_p05deg.nc", variable = "quality", timeandloc = FALSE, FUN = qualfun, add.args = list(name = "quality", noBits=16))

# Construction of the Aridity Index from highres sets: Precip/PotEvap (because of sign conventions (towards surface) PotEvap has a negative value in the data). Precip has 596 Missing values. It is a climatological aridity, so the timeseries are averaged per location, by function tempavg()
avg_precip_p05deg <- extract_var(filepath = "Precip_2000-2014_p05deg.nc", variable = "Precip", timeandloc = TRUE, FUN = tempavg, add.args = list(name = "Precip"))
avg_potevap_p05deg <- extract_var(filepath = "PotEvap_2000-2014_p05deg.nc", variable = "PotEvap", timeandloc = TRUE, FUN = tempavg, add.args = list(name = "PotEvap"))
# actual computation in mutate()
arid_p05deg <- left_join(x = avg_potevap_p05deg, y=avg_precip_p05deg, by= c("cellnr", "time")) %>% mutate(aridityindex = Precip / (- PotEvap)) %>% dplyr::select(time, cellnr, aridityindex) # Regular select() might be masked by Raster.

# The aridity datasets are shorter (only until 2014-12-01), therefore a left_join to the evi/quality (till 2017-09-01) set is appropriate.
Set1_p05deg <- left_join(x = evi_p05deg %>% mutate(quality = qual_p05deg$quality), y=arid_p05deg)

# The groundwater data are not timeseries, while the time dimension is very present in the current Set1_p05deg. Therefore the gw data will be added only in the analysis where EVI is also converted to summary statistics. Otherwise for one cell the same gw entry would appear each month, which is inefficient storage.

# Save set and clean memory.
save(Set1_p05deg, file = "./Set1_p05deg.RData")
rm(list = ls(pattern = "*_p05deg"))
gc()

# ===========================
# Part 3 Build dataset 2
# ===========================

# This is wrapped in a function because we need to be able to write the normalized timeseries, but also do the same for a probabilisation of the timeseries (conversion of values to empirical non-exceedence probabilities) and also to output them as untransformed variables.
# The name_of_values options are: "norm", "prob", "untr"
create_set2 <- function(FUN=normalisation, name_of_values = "norm") {

  print(paste("START build of dataset 2 at 0.25 degree with", as.list(match.call())[[2]] ,"of the timeseries", sep = " "))
  
  # Create list of models and of variables. And apply the extraction over them to create the datasets.
  insts <- c("anu", "polytechfr", "ecmwf", "uu") # leads to models: W3RA, ORCHIDEE, H-TESSEL, PCR-GLOBWB (ceh, nerc)
  vars <- c("PotEvap", "Precip", "TVeg", "TotMoist")
  
  # Create a dummy tibble of the right length to put each extracted variables in. They can just be added and do not need a merge by unique location/time combinations as the data is stored in netcdfs of exactly equal dimensions.
  extracted <- extract_var(filepath = "e2o_anu_PotEvap_2000-2014_p25deg.nc", variable = vars[1], timeandloc = TRUE) %>% select(x,y,cellnr,time)
  
  # Collection of the variables, plus extra corrections where sign-conventions were not followed (goddamnit): anu-TVeg, polytechfr-TVeg and polytechfr-PotEvap
  for (inst in insts) {
    for (var in vars) {
      filepath <- paste0("e2o_", inst, "_", var, "_2000-2014_p25deg.nc")
      print(filepath)
      if (identical(character(0), list.files(pattern = filepath))) {
        print("combination does not exist")
      } else if (var == vars[4]) { # For TotMoist all institutions have properly stored the data. Only transformation might be needed
        extracted[[paste(inst, var, sep = ".")]] <- extract_var(filepath = filepath, variable = var, timeandloc = TRUE, FUN = FUN, add.args = list(name = var, reverse = FALSE))[[var]]
      } else if (var == vars[2]) { # For Precip all is properly stored (positive downward). No transformation in any case, because later on the aridity variable is constructed from: Precip - PotET
        extracted[[paste(inst, var, sep = ".")]] <- extract_var(filepath = filepath, variable = var, timeandloc = TRUE)[[var]]
      } else if (var == vars[1]) { # For PotEvap, only polytechfr has stored it unconventionally. The rest stored it as positive downward. Here these three are directed upward for the aridity variable: Precip - PotET. No transformation yet.
        if ( ! inst == insts[2]) {
          extracted[[paste(inst, var, sep = ".")]] <- - extract_var(filepath = filepath, variable = var, timeandloc = TRUE)[[var]]
        } else {
          extracted[[paste(inst, var, sep = ".")]] <- extract_var(filepath = filepath, variable = var, timeandloc = TRUE)[[var]]
        }
      } else { # For TVeg, both anu and polytechfr have stored it unconventionally. The rest stored it as positive downward. Here the last two are directed upward (with the reverse-option in the supplied FUN function)
        if ( ! inst %in% insts[1:2]) {
          extracted[[paste(inst, var, sep = ".")]] <- extract_var(filepath = filepath, variable = var, timeandloc = TRUE, FUN = FUN, add.args = list(name = var, reverse = TRUE))[[var]]
        } else {
          extracted[[paste(inst, var, sep = ".")]] <- extract_var(filepath = filepath, variable = var, timeandloc = TRUE, FUN = FUN, add.args = list(name = var, reverse = FALSE))[[var]]
        }
      }
    }
  }

  # Now the aridity variable Precip-PotEvap can be created. Remember that originally the precipitation datasets of polytechfr and anu have not been present on the portal. And that in the download_e2o_script.sh we have labeled an averaged highres MSWEP-precipitation as their input (as this also was the prescribed input for all these tier 2 model runs). The aridity variable is possibly transformed in the second line of this loop.
  for (inst in insts) {
    temp <- extracted[[paste(inst, "Precip", sep = ".")]] - extracted[[paste(inst, "PotEvap", sep = ".")]] # Computation of the variable Prec-PET 
    extracted[[paste(inst, "Precip-PotEvap", sep = ".")]] <- FUN(frame = extracted %>% mutate(temporary = temp), name = "temporary")
  }
  
  # Now we want to further tidy the data and select variables for the final set (including 1:4 as x,y,cellnr,time). Then create institution variable.
  gath_p25deg <- extracted[c(1:4, grep(".T", x=names(extracted)), grep(".Precip-PotEvap", x=names(extracted)))] %>% 
  mygather(key.col = "inst.var", val.col = name_of_values, gather.cols = -(1:4)) %>% 
  separate(col=inst.var, into = c("inst", "var"), sep = "[:punct:]", extra = "merge") # only split at the first '.' and not at the '-' in 'Precip-Potevap'
  
  # Add the EVI (averaged to 0.25 deg) to the dataset with an rbind (after timesteps beyond 2014-12-01 are removed, because these are not available for the models)
  evi_p25deg <- extract_var(filepath = "stack_p25deg.nc", variable = "evi", timeandloc = TRUE, FUN = ifelse(test = (name_of_values == "untr"), yes = evifun, no = FUN), add.args = list(name = "evi")) # The ifelse() is to make sure that evi is still scaled to its original range when no transformation takes place, switch is based on name_of_values.
  # For the rbind to work, the variable needs to be named, and it needs an institution entry: which is NA because it is not modeled but measured.
  evi_p25deg <- evi_p25deg %>% mutate(inst = as.character(NA), var = "Evi")
  evi_p25deg[[name_of_values]] <- evi_p25deg$evi # naming. Other values have been named in the mygather() above
  evi_p25deg <- evi_p25deg %>% select(-evi) %>% filter(time <= "2014-12-01")
  
  # Add the SILO data, and do the appropriate transformation according the supplied FUN and name_of_values
  silo_rain_p25deg <- extract_var(filepath = "silo_Precip_2000-2015.nc", variable = "rain", timeandloc = TRUE, FUN = NULL)
  silo_PotEvap_p25deg <- extract_var(filepath = "silo_PotEvap_2000-2015.nc", variable = "evap", timeandloc = TRUE, FUN = NULL)
  silo_transformed_p25deg <- full_join(x = silo_rain_p25deg, silo_PotEvap_p25deg) %>% 
    mutate(temp = rain + evap , inst = as.character(NA), var = "Precip-PotEvap") # Because of sign conventions. PotEvap was stored positive downward.
  silo_transformed_p25deg[[name_of_values]] <- FUN(frame = silo_transformed_p25deg, name = "temp")
  # make sure names match for rbind and temporal removal of anything outside 2000-02-01:2014-12-01
  silo_transformed_p25deg <- silo_transformed_p25deg %>% 
    filter(("2000-02-01" <= time) & (time <= "2014-12-01")) %>%
    select(-rain, -evap, -temp)
  
  
  gath_p25deg <- rbind(gath_p25deg, evi_p25deg, silo_transformed_p25deg)
  return(gath_p25deg)
}

# Create the (un)transformed results by calling the wrapped function:
Set2_p25deg_norm <- create_set2(FUN = normalisation, name_of_values = "norm")
Set2_p25deg_prob <- create_set2(FUN = probabilisation, name_of_values = "prob")
Set2_p25deg_orig <- create_set2(FUN = notransformation, name_of_values = "untr") # EVI still needs scaling even if there is no transformation, the function create_set2 therefore has a that switch depends on the name_of_values.
Set2_p25deg <- full_join(x = full_join(x = Set2_p25deg_orig, y = Set2_p25deg_norm), y = Set2_p25deg_prob)

# Save the set.
save(Set2_p25deg, file = "./Set2_p25deg.RData")

# ===========================
# Part 3 Check consistency of dataset 2 timeseries at a Single spot. 
# ===========================

#library(gridExtra)
#plot1 <- ggplot(Set2_p25deg %>% filter(cellnr==400), aes(x=time, y=untr)) + geom_line(aes(color=inst), show_guide = FALSE) + facet_wrap(~ var, scales = "free_y", ncol = 1) + labs(title = "cellnr400 original")
#plot2 <- ggplot(Set2_p25deg %>% filter(cellnr==400), aes(x=time, y=norm)) + geom_line(aes(color=inst), show_guide = FALSE) + facet_wrap(~ var, scales = "free_y", ncol = 1) + labs(title = "cellnr400 normalized")
#plot3 <- ggplot(Set2_p25deg %>% filter(cellnr==400), aes(x=time, y=prob)) + geom_line(aes(color=inst)) + facet_wrap(~ var, scales = "free_y", ncol = 1) + labs(title = "cellnr400 probabilized")
#grobs <- list(plot1, plot2, plot3)
#grid.arrange(grobs=grobs, layout_matrix = matrix(data = c(1,2,3,3), nrow = 1))
## Spatial mean response. Can imply an averaging of probabilities. BAD BAD PRACTICE!
#plot4 <- Set2_p25deg %>% group_by(time, inst, var) %>% summarise(spatmean = mean(normval, na.rm = TRUE)) %>% ggplot(aes(x=time, y=spatmean)) + geom_line(aes(color=inst)) + facet_wrap(~var, ncol = 1) + labs(title = "mean")

# ===========================
# Part 4 Write the DEM datasets.
# ===========================

# Cleanup because no functions are needed and the datasets get crazy in size: 1.5 Gb
rm(list = ls())
gc()

# The elevation value is an integer and in meters.
Set3_dem_highres <- as.tibble(as.data.frame(stack(x = "./merged_dem_highest_res.nc", varname = "Band1"), xy=TRUE))
names(Set3_dem_highres) <- c("x", "y", "Elev")
saveRDS(Set3_dem_highres, file = "./Set3_dem_highres.rds")
Set3_dem_p05deg <- as.tibble(as.data.frame(stack(x = "./merged_dem_p05deg.nc", varname = "Band1"), xy=TRUE))
names(Set3_dem_p05deg) <- c("x", "y", "Elev")
saveRDS(Set3_dem_p05deg, file = "./Set3_dem_p05deg.rds")
