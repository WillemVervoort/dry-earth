#!/usr/bin/env Rscript

# Script for analysis 1
# The workflow is as follows: the Set1 dataset at p05deg contains EVI timeseries for each location. With several summary metrics (see methodology and table 1 of the report) we characterise a certain drought response and take out the temporal dimension.
# Set1 also already includes the aridity index from highres meteo sets. We add the groundwater observations to basically construct have a spatial collection of three continuous variables usable for regression: groundwater depth, climatological aridity.
# In the last part of the script these are explored, figures produced and summary statistics computed.

require(tibble)
require(lubridate)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

# ===========================
# Part 1 function definitions
# ===========================

# Function for a neirest neighbour assignment of the groundwater observations to the correct cellnr.
# Input is a dataframe with "cellnr", the "x" longitude and the "y" latitude of the cell midpoints, and the cellsize in decimal degrees.
# Input is also a frame with "x" lats and "y" lons of the point observations
# Output is a the same point observation frame with an extra column of the cellnr the observation falls in. Can be used for averaging later (e.g. when multiple are present in one cell.)
nn_assignment <- function(cellframe, obsframe, cellsize = 0.05) {
  stopifnot(all(has_name(cellframe, c("x", "y", "cellnr"))) & all(has_name(obsframe, c("x", "y"))))
  
  increment = cellsize/2
  # Create four extra colums with boundary longitudes and latitudes (ymax, ymin, xmax, xmin) by adding the increment in all directions.
  cellframe <- cellframe %>% mutate(ymax = y + increment, ymin = y - increment, xmax = x + increment, xmin = x - increment)
  result <- obsframe %>% rowwise() %>% # rowwise() practically is the same as creating a loop over x and y in the obsframe, this enables subsetting by boolean logic: 
    mutate(nncellnr = cellframe$cellnr[which((cellframe$xmax >= x & cellframe$xmin < x) & (cellframe$ymax >= y & cellframe$ymin < y))]) 
  return(result %>% ungroup())
}

# Function for conditional mutating.
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

# ===========================
# Part 2 Build set for analysis by summarizing EVI series and adding groundwater. 
# ===========================

load(file = "../DATA/Set1_p05deg.RData") # Loads the spatio-temporal EVI series

class_p05 <- read.csv("../DATA/spatial_classification_p05deg.csv", stringsAsFactors = FALSE) # spatial classification from pre-processing step 4

# Visual inspection: (can be uncommented if desired)
# Set1_p05deg %>% filter(cellnr %in% c(1,70,140)) %>% # rowwise increase of cellnr. Should be a transect, and show an aridity or EVI gradient
#   ggplot(aes(x=time, y = evi)) + geom_line(aes(color = factor(cellnr)))

# == 2.1 == Quality control

# There are 190629 low quality pixels. This is determined by the first two bits in the quality sequence, if "00" then the quality is good. See MODIS metadata document table 5. The low quality ones are either not produced or 'most probably cloudy'. Note that Analysis 2 currently does not have a quality control.
Set1_p05deg <- Set1_p05deg %>% mutate(goodqual = ifelse(test = (substr(x= quality, start = 1, stop = 2) == "00"), yes = TRUE, no = FALSE)) %>% 
  mutate(evi = ifelse(goodqual, evi, as.numeric(NA))) %>% #  We set these 190629 evi-values to NA
  select(-quality, -goodqual)

# == 2.2 == Convert the EVI timeseries to summary statistics.

th_nr_mo <- 4 # Threshold for how many monthly values should at least be present (within 12 months) to compute a yearly value 
th_nr_yr <- 6 # Threshold for how many yearly values should at least be present (within 17 full years and 10 years of drought) to compute a summary statistic

# Intra-annual variation: average of the within year Standard deviation in monthly EVI values
evi_intra <- Set1_p05deg %>% mutate(yr = year(time)) %>% group_by(cellnr, yr) %>% 
  summarise(stdev = ifelse(test = (sum(!is.na(evi)) >= th_nr_mo), yes = sd(evi, na.rm = TRUE), no = as.numeric(NA))) %>% 
  group_by(cellnr) %>% # Overrides previous grouping
  summarise(evi_intra = ifelse(test = (sum(!is.na(stdev)) >= th_nr_yr), yes = mean(stdev, na.rm = TRUE), no = as.numeric(NA))) # Computes the average of the Standard deviation. However, places NA when the number of present years is lower than the threshold.
# In the end this results in some negative values, because of negative mean evi's in certain years (and thus negative COV in that year) Leave uncorrected now STILL THE CASE WITH STDEV?

# na_set <- Set1_p05deg %>% mutate(yr = year(time)) %>% group_by(cellnr, yr) %>% summarise(nrna = sum(is.na(evi))) # there are quite some 

# Inter-annual variation: Standard deviation in annual average EVI values.
evi_inter <- Set1_p05deg %>% mutate(yr = year(time)) %>% group_by(cellnr, yr) %>% 
  summarise(ann_avg = ifelse(test = (sum(!is.na(evi)) >= th_nr_mo), yes = mean(evi, na.rm = TRUE), no = as.numeric(NA))) %>% 
  group_by(cellnr) %>%
  summarise(evi_inter = ifelse(test = (sum(!is.na(ann_avg)) >= th_nr_yr), yes = sd(ann_avg, na.rm = TRUE), no = as.numeric(NA))) # Automatically 

# Average of yearly means. Not equal to overall mean because group sizes differ
evi_mean <- Set1_p05deg %>% mutate(yr = year(time)) %>% group_by(cellnr, yr) %>% 
  summarise(ann_avg = ifelse(test = (sum(!is.na(evi)) >= th_nr_mo), yes = mean(evi, na.rm = TRUE), no = as.numeric(NA))) %>%
  group_by(cellnr) %>%
  summarise(evi_mean = ifelse(test = (sum(!is.na(ann_avg)) >= th_nr_yr), yes = mean(ann_avg, na.rm = TRUE), no = as.numeric(NA)))

# Minimum of yearly means.
evi_min <- Set1_p05deg %>% mutate(yr = year(time)) %>% group_by(cellnr, yr) %>% 
  summarise(ann_avg = ifelse(test = (sum(!is.na(evi)) >= th_nr_mo), yes = mean(evi, na.rm = TRUE), no = as.numeric(NA))) %>%
  group_by(cellnr) %>%
  summarise(evi_min = ifelse(test = (sum(!is.na(ann_avg)) >= th_nr_yr), yes = min(ann_avg, na.rm = TRUE), no = as.numeric(NA)))

# Maximum of yearly means.
evi_max <- Set1_p05deg %>% mutate(yr = year(time)) %>% group_by(cellnr, yr) %>% 
  summarise(ann_avg = ifelse(test = (sum(!is.na(evi)) >= th_nr_mo), yes = mean(evi, na.rm = TRUE), no = as.numeric(NA))) %>%
  group_by(cellnr) %>%
  summarise(evi_max = ifelse(test = (sum(!is.na(ann_avg)) >= th_nr_yr), yes = max(ann_avg, na.rm = TRUE), no = as.numeric(NA)))

# Change/recovery growth in mean evi from before 2010 (inside drought) to the period after 2011 (fully recovered)
evi_rec <- Set1_p05deg %>% mutate(before = time <= ymd("2009-12-31"), after = time >= ymd("2011-01-01")) %>%
  mutate(place = ifelse( test = (!before & !after), yes = "inter", no = ifelse( test = before, yes = "before", no = "after"))) %>%
  group_by(cellnr, place) %>% summarise(avg_temp = mean(evi, na.rm = TRUE), before = unique(before), after = unique(after)) %>%
  group_by(cellnr) %>%
  do((.[.$after,3] - .[.$before,3])) %>%  # for each cellnr: (new-old)/old * 100 percent #/.[.$before,3]*100
  rename(evi_rec = avg_temp) %>% ungroup()

# == 2.3 == Flatten temporal dimension of the aridity data and merge with summary statistics

# although the ariditity index has the samen value each timestep (temporal summary per cell), it has only been written inside the set until 2014-12-01. So we can just select the first timestep
cell_arid <- Set1_p05deg %>% select(x,y,cellnr,aridityindex,time) %>% 
  filter(time == unique(Set1_p05deg$time)[1]) %>%
  select(-time) # Removes all temporal dimension. Midpoints cellnr and aridity do not change

cell_arid <- full_join( x = full_join( x = full_join( x = full_join(x = full_join(x = full_join(x = cell_arid, y=evi_mean), y = evi_max), y = evi_min), y = evi_inter), y = evi_intra), y = evi_rec)

# Some cells lie in lakes, these generally at some point obtain negative values and due to salts show an opposite temporal trajectory.
# One of these large values (19.6) has a weird timeseries, also in the negative. Apparently it lies in Lake Numalla in Queensland:
#Set1_p05deg %>% filter(cellnr == 11557) %>% ggplot(aes(x = time, y = evi)) + geom_line() 

# Therefore we select those by a negative annual average and put NA:
cell_arid <- cell_arid %>% mutate_cond(condition = (is.na(evi_min) | evi_min < 0), evi_mean = as.numeric(NA), evi_max = as.numeric(NA), evi_min = as.numeric(NA), evi_inter = as.numeric(NA), evi_intra = as.numeric(NA), evi_rec = as.numeric(NA))

# == 2.4 == Join the spatial classification to the dataset
cell_arid <- left_join(x = cell_arid, y = class_p05)

# == 2.5 == Load groundwater data. Fill in missing Elevation values with the highres DEM. And add to the EVI metrics frame

all_gw_data <- readRDS(file = "../DATA/All_gw_data.rds")  %>% mutate(x = Lon, y=Lat)  # Loading and short renaming for the nn_assignment-function
# Although the precision of the groundwater observation coordinates has been improved, we still do not employ the highresolution DEM to look up elevations. We use the same resolution as the EVI and aridity dataset and therefore only have to do a single nn-assignment.
Set3_dem_p05deg <- readRDS(file = "../DATA/Set3_dem_p05deg.rds") %>% mutate(cellnr = 1:n()) # Loading of the DEM and assigning cellnrs in the zonal direction first.

# The groundwater data are also temporal summaries, and can be assigned by assigning the observations to the NN-gridpoint. However, there can be two records (SWL and SWL.longterm) and multiple observations per cell.
all_gw_nn <- nn_assignment(cellframe = cell_arid, obsframe = all_gw_data, cellsize = 0.05) # this results in 321 unique cellnrs. For the other observations we will have to average. The maximum of nobs in a single cell is apparently 33

# Give supplementary elevation as the variable 'Elev2' to replace 'Elev' where nothing is measured.
all_gw_nn$Elev2 <- vapply(X = all_gw_nn$nncellnr, FUN.VALUE = 1, FUN = function(x) Set3_dem_p05deg$Elev[which(Set3_dem_p05deg$cellnr == x)])
all_gw_nn <- all_gw_nn %>% mutate_cond(condition = is.na(Elev), Elev = Elev2)

gw_assigned_to_cell <- all_gw_nn %>% select(nncellnr, Elev, Depth.Compl, Date.Compl, Type, Period, N, SWL) %>% 
  group_by(nncellnr, Type, Period) %>% # This says: for all observations linked to a certain nncellnr (so multiple when two or more lie in the same cell): do the following
  summarise_all(funs(median(., na.rm = TRUE))) %>% # The median method is chosen because it makes more sense for the Dates.
  rename(cellnr = nncellnr) %>% ungroup() # Rename to variable that can be used for joining to the EVI summary statistics set.

# right_join means that all cellnrs in gw_assigned are taken from cell_arid
cell_complete <- full_join(x = cell_arid, y = gw_assigned_to_cell)
cell_filtered <- right_join(x = cell_arid, y = gw_assigned_to_cell)

# ===========================
# Part 3 Actual analysis
# ===========================

# == 3.0 == Maps for Willem
map_arid <- ggplot(cell_complete, aes(x = x, y = y)) + geom_point(aes(color = aridityindex), shape = 15) + scale_colour_gradientn(colours = brewer.pal(n = 11, name = "RdYlBu")) + geom_point(data = all_gw_nn, mapping = aes(shape = Type), color = "black") +
 scale_shape_manual(values = c(1,2)) + labs(shape = "groundwater data type", title = "Climatological AI (Prec/PotEvap) 2000-2014")

map_mean <- ggplot(cell_complete, aes(x = x, y = y)) + geom_point(aes(color = evi_mean), shape = 15) + geom_point(data = all_gw_nn, mapping = aes(shape = Type), color = "red") +
  scale_shape_manual(values = c(1,2)) + labs(shape = "groundwater data type", title = "Mean of annual average EVI 2000-2017")

map_max <- ggplot(cell_complete, aes(x = x, y = y)) + geom_point(aes(color = evi_max), shape = 15) + labs(title = "Maximum of annual average EVI 2000-2017")

map_min <- ggplot(cell_complete, aes(x = x, y = y)) + geom_point(aes(color = evi_min), shape = 15) + labs(title = "Minimum of annual average EVI 2000-2017")

map_inter <- ggplot(cell_complete, aes(x = x, y = y)) + geom_point(aes(color = evi_inter), shape = 15) + labs(title = "Stdev in annual average EVI values 2000-2017")

map_intra <- ggplot(cell_complete, aes(x = x, y = y)) + geom_point(aes(color = evi_intra), shape = 15) + labs(title = "Average of yearly Stdev in monthly EVI values 2000-2017")

map_rec <-  ggplot(cell_complete, aes(x = x, y = y)) + geom_point(aes(color = evi_rec), shape = 15) + labs(title = "Recovery in percent from mean EVI 2000-2009 to mean EVI 2011-2017") +
 scale_color_gradient2(limits = c(-60, 60))

map_class <- ggplot(cell_filtered, aes(x = x, y = y)) + geom_point(aes(color = level2, shape = level1)) + labs(title = "classification of the groundwater observations") + scale_shape_discrete(na.value = 15)

# == 3.1 == Exploratory plots

#plot_min <- ggplot(cell_complete %>% arrange(!is.na(SWL.longterm)), aes(x = aridityindex, y = evi_min)) + geom_point(aes(color = SWL.longterm, alpha = !is.na(SWL.longterm))) + scale_color_continuous(limits = c(0,40)) + labs(color = "SWL [m]", alpha = "SWL measured?", x = "aridity index", y = "Min annual EVI") 

plot_mean <- ggplot(cell_complete %>% arrange(!is.na(SWL)), aes(x = aridityindex, y = evi_mean)) + geom_point(aes(color = SWL, alpha = !is.na(SWL))) + scale_color_continuous(limits = c(0,50)) + labs(color = "SWL [m]", alpha = "SWL measured?", x = "aridity index", y = "Mean annual EVI") 
# There seems to be a cluster on top of the linear trend: Places where higher mean EVI's are reached because groundwater is shallow. The linear trend itself shows that drier areas are in the absolute sense on average less vegetated, throughout the whole timeseries.
pdf(file = "./An1_mean_cont.pdf", width = 7, height = 4.5)
plot_mean
dev.off()

plot_intra <- ggplot(cell_complete %>% arrange(!is.na(SWL)), aes(x = aridityindex, y = evi_intra)) + geom_point(aes(color = SWL, alpha = !is.na(SWL))) + scale_color_continuous(limits = c(0,50)) + labs(color = "SWL [m]", alpha = "SWL measured?", x = "aridity index", y = "Average intra-annual EVI stdev") 
pdf(file = "./An1_intra_cont.pdf", width = 7, height = 4.5)
plot_intra
dev.off()

plot_inter <- ggplot(cell_complete %>% arrange(!is.na(SWL)), aes(x = aridityindex, y = evi_inter)) + geom_point(aes(color = SWL, alpha = !is.na(SWL))) + scale_color_continuous(limits = c(0,50)) + labs(color = "SWL [m]", alpha = "SWL measured?", x = "aridity index", y = "Inter-annual EVI stdev")
# The interannual COV of annual averages were very much influenced by the fact that in this short timeseries drought is present for about halve the time. The extreme values of the drought influence the mean (reduce it) and are then suddenly less extreme from the Standard deviation point of view.
pdf(file = "./An1_inter_cont.pdf", width = 7, height = 4.5)
plot_inter
dev.off()

# Filtered version for the recovery plot (would mostly have to do with )
plot_rec <- ggplot(cell_complete %>% filter(Type == "series" & Period == "outside"), aes(x = aridityindex, y = evi_rec)) + geom_point(aes(color = SWL)) + scale_color_continuous(limits = c(0,50)) + labs(color = "SWL [m]", x = "aridity index", y = "EVI rec")

# == 3.2 == Spatial separation of the plots

# Full EVI series regionalized:
ggplot(cell_complete %>% arrange(!is.na(level2)), aes(x = aridityindex, y = evi_mean)) + geom_point(aes(color = level2)) # Large difference between east of darling and west of darling

# only for groundwater points
ggplot(cell_complete %>% arrange(!is.na(level2)), aes(x = aridityindex, y = evi_mean)) + geom_point(aes(color = level2, size = SWL)) + scale_size(trans = 'sqrt', limits = c(0,50), breaks = c(0,2,5,10,30,40)) + labs(title = "mean EVI in points with measured SWL, regions are colored.")

# In intra we really see the surface water supply stand out e.g. darling relative to east of and westof. and Macquarie
ggplot(cell_complete %>% arrange(!is.na(level2)), aes(x = aridityindex, y = evi_intra)) + geom_point(aes(color = level2, size = SWL)) + scale_size(trans = 'sqrt', limits = c(0,50), breaks = c(0,2,5,10,30,40)) + labs(title = "average of intra-annual stdev in EVI with measured SWL, regions are colored.")
ggplot(cell_complete %>% arrange(!is.na(level2)), aes(x = aridityindex, y = evi_intra)) + geom_point(aes(color = level2)) # Baselevel changes

# For inter we see most differences disappear. No difference east of darling west of darling despite differences in mean.
ggplot(cell_complete %>% arrange(!is.na(level2)), aes(x = aridityindex, y = evi_inter)) + geom_point(aes(color = level2))

# For recovery, we see that westofdarling recovers the most. And its groundwater measurements aren't that shallow. Actually a lot deeper than the others on average: cell_filtered %>% group_by(level2) %>% summarize(meanswl = mean(SWL))
ggplot(cell_complete %>% arrange(!is.na(level2)), aes(x = aridityindex, y = evi_rec)) + geom_point(aes(color = level2, size = SWL)) + scale_size(trans = 'sqrt', limits = c(0,50), breaks = c(0,2,5,10,30,40)) + labs(title = "recovery in EVI with measured SWL, regions are colored.") + ylim(0,0.06)

# == 3.3 == Fit different models to the recovery

# Only for evi_rec the pattern might be stratified enough. Here probably something can be said about the influence of groundwater on 
split_level <- 10
gw_shallow <- lm(evi_rec ~ aridityindex, data = cell_filtered %>% filter(SWL.longterm < split_level))
gw_deep <- lm(evi_rec ~ aridityindex, data = cell_filtered %>% filter(split_level  <= SWL.longterm))

pal1 <- brewer.pal(3, "Set1")

plot_rec <- ggplot(data = cell_filtered, aes(x = aridityindex, y = evi_rec)) + geom_point(aes(color = SWL.longterm >= split_level)) + scale_color_brewer(palette = "Set1") + geom_abline(intercept = gw_shallow$coefficients[1], slope = gw_shallow$coefficients[2], color = pal1[1]) + geom_abline(intercept = gw_deep$coefficients[1], slope = gw_deep$coefficients[2], color = pal1[2]) + labs(color = paste("SWL deeper than", split_level, "m", sep = " "), y = "Recovery of mean EVI after 2010 [percent]", x = "aridity index")
pdf(file = "./An1_rec_10m.pdf", width = 7, height = 4.5)
plot_rec
dev.off()

# A test for significant difference between deep and shallow models is not done with an interaction term (https://stats.stackexchange.com/questions/55501/test-a-significant-difference-between-two-slope-values) but by evaluating the sampling distrubution of difference in the intercept. This is a linear combination of two (assumed independent) sampling distributions with: sd_diff = sqrt( sd_1^2 + sd_2^2) and mu_diff = mu_1 - mu_2
# We will test whether deep has a higher intercept than shallow (and the positive difference differs significantly from zero.)
diff_mean <- gw_deep$coefficients[1] - gw_shallow$coefficients[1]
diff_sd <- sqrt( coef(summary(gw_deep))["(Intercept)", "Std. Error"]^2 + coef(summary(gw_shallow))["(Intercept)", "Std. Error"]^2)
qnorm(p = c(0.05, 0.1), mean = diff_mean, sd = diff_sd) # In the first case zero does not lie in the lower tail, in the second case it does. Therefore, at the 0.1 confidence level: a g.w. level below 10 m results in a larger recovery at the end of the drought than when it is shallower and the drought had less of an impact.

# However, probably the assumptions for a Least Squares Estimation are not valid. Expected value and variance of the residuals seem dependend on x. Non-linear model?  
plot(y = residuals.lm(gw_shallow), x = gw_shallow$model$aridityindex)
plot(y = residuals.lm(gw_deep), x = gw_deep$model$aridityindex)
