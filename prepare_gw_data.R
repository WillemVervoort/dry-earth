#!/usr/bin/env Rscript

# This script constructs the groundwater dataset from several raw databases, whose setup differs for Queensland and New South Wales.
require(tibble)
require(dplyr)
require(tidyr)
require(ggplot2)
require(lubridate)
setwd("~/Australia/Raw_gw_data/")

# ===========
# Part 0: Content of the dataset to be created
# ===========
# The starting point are unique "Work.No" for which we know that its borehole/well lies in our domain. Elevation and the completed depth are reported when available. The first type is a 'point' type. The Standing Water Level is therefore assiciated with the reported (bore completion) date (which we require to be in the study period, as SWL's from the early 1900's will have little meaning).  The second type is a 'series' from either the full study period, the period outside the drought or inside the drought and contians the average SWL. For each of these temporal aggregations the number of values on which the median was based is reported.

# Both Queensland data and new south wales data have some entries for the top values of the first waterbearing layers. We use this possible filter to exclude all observations from likely confined aquifers below this threshold in meters:
topthreshold <- 50

datanames <- c("Work.No", "Lon", "Lat", "Elev", "Depth.Compl", "Date.Compl", "Type", "Period", "N", "SWL")

# ===========
# Part 1: Function definitions
# ===========

# Function for conditional mutating.
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

# Filter of groundwater measurements. Simplest case is one point measurement, in which case only the observations inside study period are taken.
# The more difficult case is a timeseries. First the values inside the study period are taken, then if the maximum differences are above a certain value (often a saw-tooth pattern that looks like pumping tests) the outliers are removed, then the overall median is taken, and the median inside and outside the drought. For all the number of observations for the median is reported.
# Changes the measurementname to 'SWL' (negative values need to be removed beforehand)
# Is supplied a tibble, returns a filtered tibble
gw_filter <- function(dataset, groupname = Work.No, datename = Date, measurementname = Standing.Water.Level) {
  study <- list(start = "1990-01-01", end = "2017-09-01")
  drought <- list(start = "2000-01-01", end = "2010-12-31")
  maxdifth <- 20
  maxiqr <- 1.5
  
  groupname <- enquo(groupname)
  datename <- enquo(datename)
  measurementname <- enquo(measurementname)
  stopifnot(all((dataset %>% pull(!!measurementname)) >= 0, na.rm = TRUE)) # Pre-requisite that values are also positively oriented for clarity of the maxdif approach. 
  
  # If large groups then we are dealing with timeseries, and not with an accidental duplicate
  dataset <- dataset %>% group_by(!!groupname) %>% mutate(N = n(), SWL = !!measurementname)
  if (max(dataset$N) < 4) {
    print("point type detected")
    dataset <- dataset %>% filter((!!datename) >= study$start & (!!datename) <= study$end) %>% # Brackets around datename due to preference of !! operator
      select(-!!measurementname) %>%
      mutate(Type = "point", Period = as.character(NA))
    return(dataset %>% ungroup)
  } else {
    print("series type detected")
    dataset <- dataset %>% filter((!!datename) >= study$start & (!!datename) <= study$end) %>% 
      mutate(maxdif = max(SWL) - min(SWL), distance = (SWL - median(SWL))/IQR(SWL)) %>% # Based on maxdif (indicator of too extreme changes (only possibly by pumping or measurement errors)) we remove extreme values when they lie ouside 1.5 * inter-quartile-range
      filter(!((maxdif > maxdifth) & (distance < -maxiqr | distance > maxiqr))) %>%
      select(-!!measurementname, -maxdif, -distance)
    # Filtering is complete, now the temporal aggregation will be done
    print(nrow(dataset))
    agr_full <- dataset %>% summarise(SWL = median(SWL), N = n(), Period = "full", Type = "series")
    agr_inside <- dataset %>% mutate(inside = ((!!datename) >= drought$start & (!!datename) <= drought$end )) %>% 
      group_by(!!groupname, inside) %>% summarise(SWL = median(SWL), N = n(), Type = "series") %>%
      ungroup %>% mutate(Period = ifelse(test = inside, yes = "inside", no = "outside")) %>% select(-inside)
    return(full_join(x = agr_full, y = agr_inside))
  }
}

# Filter does not work well for temporal series of New South Wales Work.No %in% c("GW036865", "GW036964"), when it is implemented with maxdif > 20 and 2*IQR

# Function for conversion of degree, min, sec to decimal degree
dms_to_dd <- function(degrees = NA, minutes = NA, seconds = NA) {
  if (any(is.na(c(degrees, minutes, seconds)))) {
    return(NA)
  } else {
    return(degrees + minutes/60 + seconds/3600)
  }
}


# ===========
# Part 1: NSW
# ===========

# For NSW: BORE_GROUNDWATER has details about drilled depth, location, measured SWL. The TEMP_WATERLEVEL has measurements through time, sometimes overlappping.

nsw_raw_bore <- as.tibble(read.csv(file = "NSW/NSW_BORE_GROUNDWATER.csv", header = TRUE, stringsAsFactors = FALSE))

# Creation of a dataset that retains all work.no's in our domain. Including the other registered information.
nsw_domain_ref <- nsw_raw_bore %>% 
  mutate(Latitude = -Latitude) %>% # Correction to negative lat coordinates.
  filter((Latitude <= -25 & Latitude >= -32.25) & (Longitude >= 140 & Longitude <= 147.75)) %>% # In this set we have 4316/4339 unique Work.No's. The duplicates are exact duplicates: nsw_crop_bore %>% group_by(Work.No) %>% filter(n() > 1) %>% select(Work.No, Standing.Water.Level) %>% pull(Standing.Water.Level) therefore: 
  distinct(Work.No, Date.Completed.X, Completed.Depth, Standing.Water.Level, Elevation, Longitude, Latitude) %>% 
  mutate(Date.Completed.X = dmy(Date.Completed.X)) %>%
  select(Work.No, Longitude, Latitude, Elevation, Completed.Depth, Date.Completed.X)

# Filter the point standing water level measurements of this Bore dataset. 
nsw_clean_bore <- nsw_raw_bore %>% 
  filter(Work.No %in% nsw_domain_ref$Work.No) %>% # Implicit spatial crop.
  mutate(Date.Completed.X = dmy(Date.Completed.X)) %>% # Date conversion
  select(Work.No, Standing.Water.Level, Date.Completed.X) %>% 
  filter(Standing.Water.Level >= 0) %>% # Filter the negative values and do not retain those without measurements.
  gw_filter(groupname = Work.No, datename = Date.Completed.X, measurementname = Standing.Water.Level) %>% # Does the temporal selection
  distinct()

# Load a second dataset that includes timeseries, and can be linked to the lat/lon information through Work numbers.
nsw_raw_temp <- as.tibble(read.csv(file = "NSW/NSW_TEMP_WATERLEVEL.csv", header = TRUE, stringsAsFactors = FALSE)) # Here the duplicate Work.No are registrations at different times. Some work.numbers seem to be random comments of a very large length. There also are empty character and of one character worknumbers. So we have to take only the ones with a reasonable length

nsw_cleaner_temp <- nsw_raw_temp %>% mutate(lworkno = nchar(Work.No)) %>% filter( 6 < lworkno & lworkno <= 8) %>% 
  select(-Pipe.No, -From.Interval.No, -To.Interval.No, -Measurement.Src.Description, -Comments, -lworkno) %>% # SWL.From.Ground is a character vector, which upon conversion has more NA's than Standing.Water.Level. When both values are present they do not differ much. However, for some entries SWL.From.Ground is available and the regular one not. So in the next step the regular NA will replaced, and both are merged to one variable (preferably not from Swl.From.Ground)
  mutate(Swl.From.Ground = as.numeric(Swl.From.Ground)) %>% 
  mutate_cond(condition = (is.na(Standing.Water.Level) & !is.na(Swl.From.Ground)), Standing.Water.Level = Swl.From.Ground) %>%
  mutate(Standing.Water.Level = replace(Standing.Water.Level, Standing.Water.Level>=999, NA)) %>% # This is to correct missing values that are not noted as NA but as 9999.00 or 99999.00
  select(-Swl.From.Ground) %>%
  mutate(Measurement.Date.X = dmy(Measurement.Date.X)) %>%  
  filter((Work.No %in% nsw_domain_ref$Work.No) & Standing.Water.Level >= 0) # Implicit spatial cropping to the domain, and filtering of zeros.

nsw_clean_temp <- gw_filter(dataset = nsw_cleaner_temp, groupname = Work.No, datename = Measurement.Date.X, measurementname = Standing.Water.Level) # Does the temporal selection, plus filtering for outliers, plus computation of medians

# #--------------
# # Short investigation into the type of timeseries. They have a large temporal extent and there are 20 that have over 100 measurements.
# nsw_cleaner_temp %>% group_by(Work.No) %>% mutate(n = n()) %>% filter(n >= 100)
# nsw_cleaner_temp %>% group_by(Work.No) %>% mutate(n = n(), med = median(Standing.Water.Level)) %>% filter(n >= 100) %>% ggplot(aes(x = Measurement.Date.X, y = Standing.Water.Level)) + geom_line(aes(color = Work.No)) + geom_point(aes(color = Work.No, y=med), stat = "identity")
# nsw_cleaner_temp %>% group_by(Work.No) %>% mutate(n = n()) %>% ggplot(aes(x = Measurement.Date.X)) + geom_histogram() # Distribution of observations through time.
# nsw_cleaner_temp %>% group_by(Work.No) %>% mutate(n = n()) %>% ggplot(aes(x = n)) + geom_histogram() # Distribution of observations through time.
# 
# # For some wells there seem to be a lot of variability (drops of more than 40 m) And there seem to be trends of decline. And seasonality
# nsw_cleaner_temp %>% group_by(Work.No) %>% mutate(n = n(), maxdif = max(Standing.Water.Level) - min(Standing.Water.Level)) %>% filter(n >= 1000 & maxdif < 50 & Measurement.Date.X > "1990-01-01") %>% ggplot(aes(x = Measurement.Date.X, y = Standing.Water.Level)) + geom_line(aes(color = Work.No))
# 
# # Perhaps we want to filter only on the maximum difference in the period 1990+ 
# nsw_cleaner_temp %>% filter(Measurement.Date.X > "1990-01-01") %>% mutate(predrought = (Measurement.Date.X < "2000-01-01")) %>% 
#   group_by(Work.No) %>% mutate(maxdif = max(Standing.Water.Level) - min(Standing.Water.Level)) %>%
#   filter(maxdif < 50) %>% #And construct a pre-drought average and post-drought average?
#   group_by(Work.No, predrought) %>% summarise(meanSWL = mean(Standing.Water.Level), n = n(), medDate = median(Measurement.Date.X)) %>% # Now keep only the averages based on 10+ values. These are 7.000 values
#   filter(n >= 10) %>% 
#   ggplot(aes(x = predrought, y = meanSWL)) + geom_boxplot()
# 
# #-----------------

# As mentioned, the NSW dataset has no Top-of-aquifer variable. But there exists a separate table that contains some water bearing layers for NSW, we run a cross-comparison to filter out those Work.No where the first water bearing layer starts very deep (> 50 m).
nsw_raw_bear <- as.tibble(read.csv(file ="NSW/NSW_WATER_BEARING_ZONES.csv", header = TRUE, stringsAsFactors = FALSE))
nsw_workno_too_deep <- nsw_raw_bear %>% group_by(Work.No) %>% 
  mutate(Top.Depth = min(From.Depth, To.Depth)) %>% select(Work.No, From.Depth, To.Depth, Top.Depth) %>% ungroup() %>%
  filter(Top.Depth > topthreshold) %>% pull(Work.No) %>% unique() # We want to cross compare (not merge) so only keep the Work.No with very deep Top.Depths. A lot of values, because no spatial crop was done.

# First the temporally filtered point measurements are joined to the nsw_domain_ref. 
# Then the aggregated timeseries are added. 
# Then some extra filter oprations are done.
nsw_merged <- rbind( right_join(x=nsw_domain_ref, y=nsw_clean_bore), right_join(x = nsw_domain_ref, y = nsw_clean_temp)) %>% # The first right join adds about 315 SWL measurements, The second right join adds 572 (but these are possibly full, inside and outside). No points without gw registration are added.
  filter(!(Work.No %in% nsw_workno_too_deep)) %>% # removes 887-760 = 127 observations where water bearing layers start too deep
  select(Work.No, Longitude, Latitude, Elevation, Completed.Depth, Date.Completed.X, Type, Period, N, SWL) # Order for renaming.

names(nsw_merged) <- datanames
# To see how point based relates to full SWL aggregate
# nsw_merged %>% group_by(Work.No) %>% filter(n() > 1) %>% arrange(Work.No) %>% filter(all(c("series", "point") %in% Type) & (Period %in% c(as.character(NA), "full"))) %>% ggplot(aes(x = Type, y = SWL)) + geom_boxplot()

# ===========
# Part 2: QLD
# ===========

# For QLD: four tables are relevant: REGISTRATION, AQUIFER, ELEVATIONS, WATER_LEVELS.TXT

qld_raw_reg <- as.tibble(read.csv(file = "QLD/REGISTRATION.csv", header = TRUE, stringsAsFactors = FALSE))

# Creation of a dataset that retains all work.no's in our domain. Including the other registered information, like latitude, drilled date. The decimal degree latitude and longitudes are derived from the deg,min,sec strings that are separated by -, 10363 strings are empty, therefore NA's are introduced.
qld_domain_ref <- qld_raw_reg %>% select(RN, DRILLED_DATE, LAT, LNG) %>% group_by(RN) %>%
  do(GIS_LAT = do.call(what = dms_to_dd, args = as.list(as.integer(unlist(strsplit(x = .$LAT, split = "-"))))), 
     GIS_LON = do.call(what = dms_to_dd, args = as.list(as.integer(unlist(strsplit(x = .$LNG, split = "-"))))),
     DRILLED_DATE = .$DRILLED_DATE) %>%
  unnest() %>% mutate(GIS_LAT = -GIS_LAT) %>% # Correction for the fact that deg,min,sec is stored positively
  filter((GIS_LAT <= -25 & GIS_LAT >= -32.25) & (GIS_LON >= 140 & GIS_LON <= 147.75)) %>% # There don't seem to be any duplicates here.
  mutate(DRILLED_DATE = dmy(DRILLED_DATE))

qld_raw_elev <- as.tibble(read.csv(file = "QLD/ELEVATIONS.csv", header = TRUE, stringsAsFactors = FALSE)) # Here there seem to be multiple measurements. PIPE == "X" seems to be the surface measurment. 
qld_clean_elev <- qld_raw_elev %>% filter(PIPE == "X") %>% # Pipe = X is the ground measurement
  mutate(ELEVATION = replace(ELEVATION, ELEVATION <= -100, NA)) %>% # This deals with the -9999 as NA and a -114 value that cannot be true.
  mutate(RDATE = dmy(RDATE)) %>% # For selection of the latest measurement
  group_by(RN) %>% mutate(maxDATE = max(RDATE)) %>% # Create dummy var. for filtering when two or more datums are available. Not enough info is available to distinguish the effect of the different datums
  filter(RDATE == maxDATE) %>% select(RN, ELEVATION) %>% ungroup() %>%
  filter(RN %in% qld_domain_ref$RN) # implicit spatial crop.

# This aquifer database can provides us with 347 SWL measurements in the domain, that are measured on the "DRILLING DATE" in the qld_domain_ref set. So first the numbers of are extracted for which the first waterbearing layer is too deep. Separately, the SWL's from the set are extracted and put through the point-measurement temporal filter.
qld_raw_aq <- as.tibble(read.csv(file = "QLD/AQUIFER.csv", header = TRUE, stringsAsFactors = FALSE)) # Duplicates are present. qld_raw_aq %>% group_by(RN) %>% summarise(nocc = n()). 
qld_rn_too_deep <- qld_raw_aq %>% filter(REC==1) %>% # Registration of the shallowest TOP for each location. This is the first registration with REC=1. We are only interested in the possibly phreatic aquifers that influence the rootzone.
  mutate(TOP = replace(TOP, TOP > 2000, NA)) %>% # Do deal with the possibly wrong registrations and nodatas of 9000+
  filter(TOP > topthreshold) %>% pull(RN) %>% unique()
qld_clean_aq <- qld_raw_aq %>% filter(REC==1) %>%
  mutate(SWL = -SWL) %>% # Correction for the fact that Queensland reports them negative downwards, while we want positive downwards.
  filter(SWL >= 0) %>% # Filter the negative values and do not retain those without measurements, similar like NSW
  mutate(RDATE = dmy(RDATE)) %>% rename(Standing.Water.Level = SWL) %>% # Rename to avoid conflict in the gw_filter function
  filter(RN %in% qld_domain_ref$RN) %>% # Implicit spatial crop.
  gw_filter(groupname = RN, datename = RDATE, measurementname = Standing.Water.Level) %>% # Currently leaves 287 observations (but possibly still to deep, filtering happens later)
  select(RN, Type, Period, N, SWL) # RDATE can be discarded because it is the same as Drilled_date
  
qld_raw_temp <- as.tibble(read.table(file = "QLD/WATER_LEVELS.txt", sep = "|", header = TRUE, stringsAsFactors = FALSE)) # Pretty large file. 300 Mb

# # short investigation into temporal range. Looks like a lot of observations suited for time-serie comparison. 
# #-------------
# qld_clean_temp <- qld_raw_temp %>% select(RN, MEASUREMENT, RDATE) %>%
#   mutate(RDATE = dmy(RDATE)) %>%
#   filter(RN %in% qld_clean_reg$RN) %>% # implicit spatial cropping
#   group_by(RN) %>% mutate(n = n())
# ggplot(data = qld_clean_water, aes(x = RDATE)) + geom_histogram() # Many observations in the range 1950-1980
# ggplot(data = qld_clean_water, aes(x = n)) + geom_histogram() # Only a few 100+ timeseries
# 
# qld_clean_water %>% filter(RDATE > "1990-01-01") %>% ungroup() %>% mutate(RN = factor(RN)) %>% filter(n > 10) %>% ggplot(aes(x = RDATE, y = MEASUREMENT)) + geom_line(aes(group = RN, color = RN))
# 
# #-------------

qld_cleaner_temp <- qld_raw_temp %>% select(RN, MEASUREMENT, RDATE) %>%
  mutate(RDATE = dmy(RDATE), MEASUREMENT = -MEASUREMENT) %>% # Correction for the fact that Queensland reports them negative downwards, while we want positive downwards.
  filter((RN %in% qld_domain_ref$RN) & MEASUREMENT >= 0) # implicit spatial crop. Most timeseries seem to have around 24 length. 

qld_clean_temp <- gw_filter(dataset = qld_cleaner_temp, groupname = RN, datename = RDATE, measurementname = MEASUREMENT) # Does the temporal selection, plus filtering for outliers, plus computation of medians

# First the elevations are added to the registration data within the domain. Then this set is joined to SWL from either the points (clean_aq) or to the aggregated timeseries (clean_temp).
# Afterwards there is a a filtering of RN's that are too deep.
qld_domain_reg <- left_join(x = qld_domain_ref, y = qld_clean_elev)
qld_merged <- rbind( right_join(x = qld_domain_reg , y = qld_clean_aq), right_join(x = qld_domain_reg , y = qld_clean_temp)) %>%
  filter(!(RN %in% qld_rn_too_deep)) %>%
  mutate(RN = as.character(RN), Depth.Compl = as.numeric(NA)) %>% # To make it match to the content of NSW
  select(RN, GIS_LON, GIS_LAT, ELEVATION, Depth.Compl, DRILLED_DATE, Type, Period, N, SWL) # Order for renaming.

names(qld_merged) <- datanames

# ===========
# Part 3: Complete set
# ===========

all_gw_data <- rbind(qld_merged, nsw_merged)
saveRDS(all_gw_data, file = "../All_gw_data.rds")
plotgw <- ggplot(all_gw_data, aes(x=Lon, y=Lat)) + geom_point(aes(color = SWL, shape = Type)) + scale_color_continuous(limits = c(0,50))
