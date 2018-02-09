#!/usr/bin/env Rscript

# Script for analysis 2. It uses the Set2 database at 0.25 degree. It has all spatio-temporal values of the modeled and observed series in it. These are both the original values, normalised values and cumulative probabilities. Each of these are used in a different part of the script for a different analysis of the drought response.
# The workflow is as follows: at the beginning of each part ancillary functions are defined (if neccesary for the analysis), then computations follow and figures or tables are exported. The observations to which we compare are either SILO Precip-PotEvap (with modeled Precip-PotEvap) or EVI (with modeled transpiration)
# Inherits the ~/Australia/RESULTS working directory from the bash script

require(tibble)
require(lubridate)
require(dplyr)
require(tidyr)
require(ggplot2)
require(hydroGOF)
require(RColorBrewer)
require(gridExtra)

# ===========================
# Part 0 Load the data
# ===========================

load("../DATA/Set2_p25deg.RData")

# Spatial classification by a left join: creates NA in the level variables when cells do not belong to a group.
class_p25 <- read.csv("../DATA/spatial_classification_p25deg.csv")
Set2_p25deg <- left_join(x = Set2_p25deg, y = class_p25)
# Order the spatial classification
Set2_p25deg$level2 <- factor(x = Set2_p25deg$level2, levels = c("westofdarling", "eastofdarling", "macquarie", "dry", "coopercreek", "warrego", "darling"))

# Add model labels
Set2_p25deg$inst <- factor(x = Set2_p25deg$inst, levels = c("ecmwf", "polytechfr", "uu", "anu"), labels = c("H-TESSEL", "ORCHIDEE", "PCR-GLOBWB", "W3RA"))

# ===========================
# Part 1 Original values
# ===========================

# Spatial and temporal pooling to compute the mean. The differences between them can inform us about whether one is wetter or drier in the absolute sense.
domainmeans <- Set2_p25deg %>% filter(var != "Evi") %>% group_by(var,inst) %>% summarise(meanval = mean(untr, na.rm=TRUE))
regionmeans <- Set2_p25deg %>% group_by(var,inst,level2) %>% summarise(meanregval = mean(untr, na.rm=TRUE), x = ymd("2012-10-01")) %>% filter(level2 %in% c("eastofdarling", "westofdarling", "macquarie"))
# Export a table that shows differences in the domain model means.
An2_ori_domainmean <- domainmeans %>% mutate(unit = ifelse(test = (var == "TotMoist"), yes = "kg m-2", no = "kg m-2 s-1" )) # the NA entry for the institution of Precip-PotEvap are the SILO observations.
write.csv(x = An2_ori_domainmean, file = "./An2_ori_domainmean.csv", row.names = FALSE)

# Graphs of the spatial mean value distribution through time for three interesting Regions. The Evi observarions and SILO observations are included. ORCHIDEE and PCR-GLOBWB seem to have the closest Precip-PotEvap estimate to SILO.
An2_ori_regset <- Set2_p25deg %>% group_by(var,inst,time, level2) %>% summarise(meanfieldval = mean(untr, na.rm=TRUE)) %>% filter(level2 %in% c("eastofdarling", "westofdarling", "macquarie"))

evi_ori_plot <- ggplot(data = An2_ori_regset %>% filter(var == "Evi"), aes(x=time, y=meanfieldval)) + 
  geom_line(col = "grey50") + 
  geom_text(data = regionmeans %>% filter(var == "Evi"), aes(label = round(x = meanregval, 2), x = x), y = 0.4, col = "grey50", size = 4) +
  facet_wrap(~level2, ncol = 3) +
  labs(x = "", y = "EVI [-]")

tveg_ori_plot <- ggplot(data = An2_ori_regset %>% filter(var == "TVeg"), aes(x=time, y=meanfieldval)) + 
  geom_line(aes(col = inst)) + 
  facet_wrap(~level2, ncol = 3) +
  geom_text(data = regionmeans %>% filter(var == "TVeg") %>% group_by(level2) %>% mutate(y = c(0.000045,0.00004,0.000035,0.00003)), aes(label = round(x = meanregval, 7), x = x, col = inst, y = y), size = 4) +
  labs(x = "", y = "TVeg [kg m−2 s−1]") + guides(color = FALSE)

pp_ori_plot <- ggplot(data = An2_ori_regset %>% filter(var == "Precip-PotEvap"), aes(x=time, y=meanfieldval)) + 
  geom_line(aes(col = inst)) + 
  facet_wrap(~level2, ncol = 3) +
  geom_text(data = regionmeans %>% filter(var == "Precip-PotEvap") %>% group_by(level2) %>% mutate(y = c(-0.0002,-0.00025,-0.0003,-0.00035,-0.0004)), aes(label = round(x = meanregval, 6), x = x, col = inst, y = y), size = 4) +
  labs(x = "", y = "Prec-PotEvap [kg m−2 s−1]") + guides(color = FALSE)

totmoist_ori_plot <- ggplot(data = An2_ori_regset %>% filter(var == "TotMoist"), aes(x=time, y=meanfieldval)) + 
  geom_line(aes(col = inst)) + 
  facet_wrap(~level2, ncol = 3) +
  geom_text(data = regionmeans %>% filter(var == "TotMoist") %>% group_by(level2) %>% mutate(y = c(900,800,700,600)), aes(label = round(x = meanregval, 0), x = x, col = inst, y = y), size = 4) +
  labs(x = "Time", y = "Tot Moist [kg m−2]") + guides(color = FALSE)

# This plot is adapted in InkScape, to align the whitespaces better and add a legend.
pdf(file = "./An2_ori_regseries.pdf", width = 11, height = 9)
grid.arrange(evi_ori_plot, tveg_ori_plot, pp_ori_plot, totmoist_ori_plot, ncol=1)
dev.off()

# # This distribution through time and space is discarded because it is less intuitive for this case than a timeserie depicting the drought.
# An2_ori_dist <- Set2_p25deg %>% ggplot(aes(x=untr)) + geom_freqpoly(aes(colour = inst, lty = is.na(inst), y = (..count..)/sum(..count..))) + facet_wrap(~var, scales = "free", ncol = 2) + labs(x = "untransformed values", y = "frequency", colour = "Models", title = "Distribution of measured (Evi) and modeled variables through time and space") + scale_linetype_manual(guide = FALSE, values = c("solid","twodash"))
# pdf(file = "./An2_ori_dist.pdf", width = 9, height = 6.5)
# An2_ori_dist
# dev.off()

# Kling-Gupta efficiency on the spatial and temporal pool. ORCHIDEE is the best. Further proof that H-TESSEL estimate deviates and that W3RA method of residual from energy balance is also different.
observed <- Set2_p25deg %>% filter(var == "Precip-PotEvap" & is.na(inst)) %>% pull(untr)
modelkge <- Set2_p25deg %>% filter(var == "Precip-PotEvap" & !is.na(inst)) %>% group_by(var,inst) %>% summarise(kge = KGE(sim = untr, obs = observed))
write.csv(x = modelkge, file = "./An2_ori_kge.csv", row.names = FALSE)

# ===========================
# Part 2 Standardized values, because seasonality and cell-climatology has to be removed for correlation.
# ===========================

# == 2.1 == Within variable autocorrelation

# Function for computation of Autocorrelation and confidence intervals. Needs to be supplied with a numerical vector that contains the timeserie
# Chatfield's Analysis of Time Series (1980) states that in the limit the variance in r is normally distributed and has a variance equal to 1/n_obs so stdev is (1/(n_obs))^0.5
# We make it into a two-tailed confidence interval. See also https://stats.stackexchange.com/questions/211628/
# It can report the full output, the maximum lag at which autocorrelation is significantly different from zero, and the last lag in which there is uninterrupted (subsequent) significant positive autocorrelation. 
acf_adapted <- function(timeserie, conf_level = 0.95, lag.max = 12, report = "full") {
  acf_output <- acf(timeserie, lag.max = lag.max, plot = FALSE)
  r_upper <- qnorm((1+conf_level)/2) * 1/sqrt(acf_output$n.used) # Build the confidence interval of a normal distribution around zero (outside these bounds the r differs significantly from zero)
  r_lower <- qnorm((1-conf_level)/2) * 1/sqrt(acf_output$n.used)
  output <- data.frame(r=acf_output$acf[,,,drop=TRUE], lag = acf_output$lag[,,,drop=TRUE], diff_from_zero = (acf_output$acf[,,,drop=TRUE] < r_lower | acf_output$acf[,,,drop=TRUE] > r_upper))
  if (report == "full") {
    return(output)
  } else if (report == "maxlag") {
    if (all(is.na(output$r))) { # short routine for timeseries that are all 0, can occur for certain cells in variable TVeg
      return(as.numeric(NA))
    } else {
      return(output$lag[max(which(output$diff_from_zero))]) # This happens for nearly all other timeseries (autocorrelation at lag 0 is always 1 so this is always significantly different from zero)
    }
  } else if (report == "sublag") {
    if (all(is.na(output$r))) { # short routine for timeseries that are all 0, can occur for certain cells in variable TVeg
      return(as.numeric(NA))
    } else {
      runlengths <- rle(output$diff_from_zero) # Provides a list of values and the length of its uninterrupted repetition in the boolean series.
      return(ifelse(runlengths$values[1], yes = output$lag[runlengths$lengths[1]], no = NA)) # We want the length of uninterrupted petition when the first series is TRUE. Then we select the lag belonging to that position in the Boolean series
    }
  } else {
    stop("provide correct way of reporting: full, maxlag or sublag")
  }
}

# Compute and plot spatial distributions (boxplot) of maximum lag with significant autocorrelation and the first sequence of subsequent lags with signigicant autocorrelation (numbers will be much smaller)
# There are cells where the original values contain 0 transpiration. They are not missing. E.g. Set2_p25deg %>% filter(var == "TVeg", inst == "ORCHIDEE", cellnr == 52). But they do result in warnings as the r becomes NaN, this is solved by putting also NA for maxlag or sublag
memory <- Set2_p25deg %>% group_by(var, inst, cellnr) %>% summarize(maxlag = ifelse(test = any(is.na(norm)), yes = as.numeric(NA), no = acf_adapted(timeserie = norm, lag.max = 180, report = "maxlag")),
  sublag = ifelse(test = any(is.na(norm)), yes = as.numeric(NA), no = acf_adapted(timeserie = norm, lag.max = 180, report = "sublag")))

# After summarizing we lose the Spatial Classication, so do a join again:
memory <- left_join(x = memory, y = class_p25)

max_box <- ggplot(data = memory, aes(x = var, y = maxlag, colour = inst)) + geom_boxplot() + labs(y = "maximum lag in months with significantly autocorrelated values", x ="Variable", colour = "Models") # The maximum lag was the previous way of doing things, but values are very high. W3RA seems to have only one type of memory in transpiration 
sub_box <- ggplot(data = memory, aes(x = var, y = sublag, colour = inst)) + geom_boxplot() + labs(y = "lag in months with uninterrupted significantly autocorrelated values", x ="Variable", colour = "Models") # Differences between max_box and sub_box interpretations are in currently in the log. Sub_box is the one currently exported for the report.

pdf(file = "./An2_submem_box.pdf", width = 7, height = 6)
sub_box
dev.off()

# Plot memory for the interesting regions.
sub_mem_set <- memory %>% filter(level2 %in% c("eastofdarling", "westofdarling", "macquarie"))
#plot_pal <- brewer.pal(n = length(unique(class_p25$level2)), name = "Set1")
#sub_class <- ggplot(data = sub_mem_set, aes(x = level2, y = sublag, colour = inst)) + geom_boxplot() + facet_wrap(~var)
#sub_class2 <- ggplot(data = sub_mem_set, aes(x = inst, y = sublag, fill = level2)) + geom_bar(stat = "summary", fun.y = "mean", position = "dodge") + facet_wrap(~var) + scale_color_manual(values = plot_pal, na.value = "grey50")
sub_class3 <- ggplot(data = sub_mem_set, aes(x = inst, y = sublag, color = inst)) + geom_boxplot() + facet_grid(level2~var)
#sub_spat <- ggplot(data = sub_mem_set, aes(x = var, y = sublag, colour = inst)) + geom_jitter(aes(shape = level2), stat = "summary", fun.y = "mean")
#sub_spat2 <- ggplot(data = memory %>% filter(!is.na(level2)), aes(x = var, y = sublag, fill = inst)) + geom_bar(stat = "summary", fun.y = "mean", position = "dodge") + facet_wrap(~level2)



# == 2.2 == Between variable (lagged) correlation: spatial pool

# Spearman rank correlation or kendalls tau between variables (an possibly for lagged combinations too)
# Function that prepares the data from the viewpoint of one variable, that is shifted with different lags by a different function (plus and minus, amount of months can be chosen) and compared to the other variables in their original form.
# It does a spatial pool (on timestep shift requires with the amount of cellnr's) with one named correlation matrix as output.
cor_norm_pool <- function(fullset, primvar = "TVeg", model = "W3RA", maxlag = 7) {
  
  stopifnot(any(unique(fullset$var) == primvar))
  stopifnot(any(unique(fullset$inst) == model))
  
  ncell <- length(unique(fullset$cellnr))
  primvars <- paste0(primvar, as.character((-maxlag):maxlag)) # Lag is applied in two ways: forward and delay
  othervars <- setdiff(unique(fullset$var), primvar)
  
  # Construct a dataframe with shifted series: the normal and lagged versions of the primary variable.
  primseries <- fullset$norm[which((fullset$var == primvar) & (fullset$inst == model))]
  stopifnot(ncell*length(unique(fullset$time)) == length(primseries))
  primshifted <- shiftseries(series = primseries, maxtimelag = maxlag, cellspertime = ncell)
  colnames(primshifted) <- primvars
  
  # Construct a dataframe with the other variables.Evi doesn't need to be filtered by model. Precip-PotEvap observations are not included, only modeled ones
  otherseries <- sapply(X = othervars, FUN = function(name) {
    if (name == "Evi") { 
      return(fullset$norm[which(fullset$var == name)])
    } else {
      return(fullset$norm[which((fullset$var == name) & (fullset$inst == model))])
    }
  })
  otherseries <- as.data.frame(otherseries)
  colnames(otherseries) <- othervars

  return(as.data.frame(cor(cbind(primshifted, otherseries), use = "complete.obs", method = "spearman")))
}

# Intermediate function that temporally lags a series. Because the series can be filled per cellnr, and we want to shift a single timestep each lag, it will fill 'fields' of size ncell with NA for each shift: 
shiftseries <- function(series, maxtimelag, cellspertime) {
  data <- as.data.frame(matrix(data = NA, nrow = length(series), ncol = length(-maxtimelag:maxtimelag)))
  appended <- c(rep(NA,maxtimelag*cellspertime), series, rep(NA,maxtimelag*cellspertime)) # first maxlag timesteps becomes NA
  for (i in 1:ncol(data)) {
    lag <- ((-maxtimelag):maxtimelag)[i]
    data[,i] <- appended[(1 + (maxtimelag-lag)*cellspertime):(length(series) + (maxtimelag-lag)*cellspertime)]
  }
  return(data)
}

# Vegetation is our main interest so primvar. Because we also have negative lags, we can see influence of other variables on TVeg and for positive shifst of TVeg on the other variables. Strictly spoken it is not influence but correlation or 'synchronicity in the timeseries'
# Then we also want a possibility to compute the spatial pooled correlation for the full field or for selected subregions, which is why the combination of correlation matrices to dataframes is wrapped in the following function:
compute_and_merge_to_frame <- function(poolset = Set2_p25deg, maximum_lag = 7, subregions = NULL) {
  if (is.null(subregions)) { # Meaning we want a full spatial pool
    labels <- "full"
  } else {
    labels <- subregions
  }
  
  # Loop of one when full spatial, of number of subregions when subregions were supplied
  corframes <- lapply(X = labels, FUN = function(label) { 
    if (label == "full") { # Either all regions are selected (even NA) when the set is full, otherwise only the subregion
      fullset <- poolset
    } else {
      fullset <- poolset %>% filter(level2 == label)
    } 
    # Compute the correlations
    w3ra_cor <- cor_norm_pool(fullset = fullset, primvar = "TVeg", model = "W3RA", maxlag = maximum_lag)
    htessel_cor <- cor_norm_pool(fullset = fullset, primvar = "TVeg", model = "H-TESSEL", maxlag = maximum_lag)
    pcr_cor <- cor_norm_pool(fullset = fullset, primvar = "TVeg", model = "PCR-GLOBWB", maxlag = maximum_lag) #  
    orch_cor <- cor_norm_pool(fullset = fullset, primvar = "TVeg", model = "ORCHIDEE", maxlag = maximum_lag) #
    
    # Put the four correlation matrices into a single frame.
    all_cor <- list("W3RA" = w3ra_cor, "H-TESSEL" = htessel_cor, "PCR-GLOBWB" = pcr_cor, "ORCHIDEE" = orch_cor)
    intermediate <- lapply(X= seq_along(all_cor), FUN = function(x) {
      var_interest <- all_cor[[x]][(1:(maximum_lag*2+1)),-(1:(maximum_lag*2+1))]
      var_interest$var_lag <- ((-maximum_lag):maximum_lag) # This is more understandable for plotting. -7 means that EVI_0 is compared to TVeg_7
      var_interest$Model <- names(all_cor)[x]
      return(gather(data = var_interest, key = "var", value = "spearmancor", -var_lag, -Model))
    })
    corframe <- do.call(what = rbind, args = intermediate)
    corframe$region <- label
    return(corframe)
  })
  return(do.call(what = rbind, args = corframes))
}

# Apply it for the complete field. 
full_cor <- compute_and_merge_to_frame(poolset = Set2_p25deg, maximum_lag = 7)

plot_full_cor <- ggplot(data = full_cor, aes(x = var_lag, y = spearmancor)) + geom_line(aes(color = Model)) + facet_wrap(~var) + labs(x = "Shift in variable timeseries [months]", y = "Correlation between the variable and transpiration")
pdf(file = "./An2_full_cor.pdf", width = 9.5, height = 5)
plot_full_cor
dev.off()

# Apply it for the three regions we are interested in. Currently used in the report.
spat_cor <- compute_and_merge_to_frame(poolset = Set2_p25deg, maximum_lag = 7, subregions = c("westofdarling", "eastofdarling","macquarie"))
plot_spat_cor <- ggplot(data = spat_cor, aes(x = var_lag, y = spearmancor)) + geom_line(aes(color = Model)) + facet_grid(region~var) + labs(x = "Shift in variable timeseries [months]", y = "Correlation between the variable and transpiration")
pdf(file = "./An2_spat_cor.pdf", width = 9.5, height = 9)
plot_spat_cor
dev.off()

# == 2.3 == Between variable (lagged) correlation: spatially distributed. Not used anymore, the spatial split is done above. Not informative to do it again per cell.

# For specific lags/variable combinations, the following function computes the correlation for each cell (899 vector combinations) and outputs a field of correlation values.
# Check: https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html for the way this function uses 'enquo' and !!:
cor_norm_cell <- function(fullset, model = "W3RA", var1 = "TVeg", lag1 = 0, var2 = "Evi") {
  
  ncell <- length(unique(fullset$cellnr))
  model <- enquo(model)
  var1 <- enquo(var1)
  var2 <- enquo(var2)
  
  if (lag1 >= 0) {
    var1set <- fullset %>% filter((inst == !!model & var == !!var1)) %>%
      mutate(var1lag = lag(x = .data$norm, n = ncell*lag1)) %>% # Lagging of var1 to calculate delayed effect of var1 on var2
      select(x, y, cellnr, time, var1lag)
  } else {
    var1set <- fullset %>% filter((inst == !!model & var == !!var1)) %>%
      mutate(var1lag = lead(x = .data$norm, n = ncell*abs(lag1))) %>% # Leading of var1 to calculate delayed effect of var2 on var1
      select(x, y, cellnr, time, var1lag)
  }
  
  if (quo_name(var2) == "Evi") { # Var 2 always stays within the same model, unless Measured Evi is chosen (has NA model entry)
    var2set <- fullset %>% filter(var == !!var2) %>% select(x, y, cellnr, time, norm)
  } else {
    var2set <- fullset %>% filter((inst == !!model & var == !!var2))  %>% select(x, y, cellnr, time, norm)
  }
  
  corrset <- var2set %>% mutate(var1lag = var1set$var1lag) %>% group_by(cellnr, x, y) %>% summarize(cellspear = cor(x = var1lag, y = norm, use = "complete.obs", method = "spearman"))
  return(corrset)
}

# Make a plot for interesting cases. This perhaps may be matched to maps of groundwater depths if this data in Analysis 1 gets better.
test <- cor_norm_cell(fullset = Set2_p25deg, model = "PCR-GLOBWB", var1 = "TVeg", lag1 = 0, var2 = "Evi")
ggplot(data = test, aes(x = x, y = y)) + geom_point(aes(color = cellspear), size = 5)

# == 2.3 == Between observations only: Correlation of Observed Precip-PotEvap (at different lags) with observed EVI 

cor_observations <- function(fullset, maxlag) {
  ppseries <- fullset %>% filter(var == "Precip-PotEvap" & is.na(inst)) %>% pull(norm) # gets to normalized SILO
  varnames <- paste0("Precip-PotEvap", as.character((-maxlag):maxlag)) 
  ncell <- length(unique(fullset$cellnr))
  
  ppshifted <- shiftseries(series = ppseries, maxtimelag = maxlag, cellspertime = ncell)
  colnames(ppshifted) <- varnames
  
  evi <- fullset %>% filter(var == "Evi") %>% select(norm) %>% rename(evi = norm) # gets to normalized EVI, rename for later use in the correlation matrix.
  
  return(as.data.frame(cor(cbind(ppshifted, evi), use = "complete.obs", method = "spearman")))
}

maximum_lag <- 7
cor_obs <- data.frame(spearmancor = cor_observations(fullset = Set2_p25deg, maxlag = maximum_lag)[(1:(maximum_lag*2+1)),-(1:(maximum_lag*2+1))], lag = (-maximum_lag):maximum_lag)
An2_ppevi_cor <- ggplot(data = cor_obs, aes(x = lag, y = spearmancor)) + geom_line() + labs(x = "shift in SILO Precip-PotEvap series in months", y = "Spearman rank correlation with measured EVI")
# These observational results are pretty counterintuitive. One would expect the observed meteorological aridity to precede the Observed EVI. Likely this follows from the way the potential evaporation is estimated.


# ===========================
# Part 3 Probability transformed values
# ===========================

# Drought defined as anything below the empirical 0.2 quantile? Then split into longer lasting droughts and shorter ones?
droughtset <- Set2_p25deg %>% mutate(dry = ifelse(prob < 0.2, yes = TRUE, no = FALSE)) %>% select(-prob, -untr, -norm)

# Timeseries of area under the influence of drought.
areaseries <- droughtset %>% group_by(time, inst, var) %>% summarise(areadry = sum(dry)/n(), sum = sum(dry), n = n())
plot_area <- areaseries %>% ggplot(aes(x = time, y = areadry)) + geom_line(aes(color = inst)) + 
  facet_wrap(~var, ncol = 1) + labs(x = "Time", y = "Dry fraction of the area", color = "Models")
pdf(file = "./An2_area_dry.pdf", width = 9.5, height = 5)
plot_area
dev.off()
