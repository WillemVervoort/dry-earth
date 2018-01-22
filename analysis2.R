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

# ===========================
# Part 0 Load the data
# ===========================

load("../DATA/Set2_p25deg.RData")

# Spatial classification by a left join: creates NA in the level variables when cells do not belong to a group.
class_p25 <- read.csv("../DATA/spatial_classification_p25deg.csv")
Set2_p25deg <- left_join(x = Set2_p25deg, y = class_p25)

# Add model labels
Set2_p25deg$inst <- factor(x = Set2_p25deg$inst, levels = c("ecmwf", "polytechfr", "uu", "anu"), labels = c("H-TESSEL", "ORCHIDEE", "PCR-GLOBWB", "W3RA"))

# ===========================
# Part 1 Original values
# ===========================

# Spatial and temporal pooling to compute the mean. The differences between them can inform us about whether one is wetter or drier in the absolute sense.
modelmeans <- Set2_p25deg %>% filter(var != "Evi") %>% group_by(var,inst) %>% summarise(meanval = mean(untr, na.rm=TRUE), meannorm = mean(norm, na.rm= TRUE), meanprob = mean(prob, na.rm=TRUE))
# Export a table that shows differences in the model means
An2_ori_mean <- modelmeans %>% select(-meannorm, -meanprob) %>% mutate(unit = ifelse(test = (var == "TotMoist"), yes = "kg m-2", no = "kg m-2 s-1" )) # the NA entry for the institution of Precip-PotEvap are the SILO observations.
write.csv(x = An2_ori_mean, file = "./An2_ori_mean.csv", row.names = FALSE)

#ggplot(data = modelmeans, aes(x=inst, y=meanval)) + geom_bar(stat = "identity") + facet_wrap(~var, scales = "free_y") # Simple histogram of the table values

# Graphs of the spatial mean value distribution through time. ORCHIDEE and PCR-GLOBWB seem to have the closest Precip-PotEvap estimate to SILO.
An2_ori_mean <- Set2_p25deg %>% group_by(var,inst,time) %>% summarise(meanfieldval = mean(untr, na.rm=TRUE)) %>% ggplot(aes(x=time, y=meanfieldval)) + geom_line(aes(color = inst)) + facet_wrap(~var, scales = "free_y", ncol = 1)
# The same goes for the graph of the original value as a distribution through time and space
An2_ori_dist <- Set2_p25deg %>% ggplot(aes(x=untr)) + geom_freqpoly(aes(colour = inst, lty = is.na(inst), y = (..count..)/sum(..count..))) + facet_wrap(~var, scales = "free", ncol = 2) + labs(x = "untransformed values", y = "frequency", colour = "Models", title = "Distribution of measured (Evi) and modeled variables through time and space") + scale_linetype_manual(guide = FALSE, values = c("solid","twodash"))
pdf(file = "./An2_ori_dist.pdf", width = 9, height = 6.5)
An2_ori_dist
dev.off()

# Kling-Gupta efficiency on the spatial and temporal pool. ORCHIDEE is the best.
observed <- Set2_p25deg %>% filter(var == "Precip-PotEvap" & is.na(inst)) %>% pull(untr)
modelkge <- Set2_p25deg %>% filter(var == "Precip-PotEvap" & !is.na(inst)) %>% group_by(var,inst) %>% summarise(kge = KGE(sim = untr, obs = observed))

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

max_box <- ggplot(data = memory, aes(x = var, y = maxlag, colour = inst)) + geom_boxplot() + labs(y = "maximum lag in months with significantly autocorrelated values", x ="Variable", colour = "Models") # W3RA seems to have only one type of memory in transpiration 
sub_box <- ggplot(data = memory, aes(x = var, y = sublag, colour = inst)) + geom_boxplot() + labs(y = "lag in months with uninterrupted significantly autocorrelated values", x ="Variable", colour = "Models")
sub_class <- ggplot(data = memory %>% filter(!is.na(level2)), aes(x = level1, y = sublag, colour = inst)) + geom_boxplot() + facet_wrap(~var)

pdf(file = "./An2_submem_box.pdf", width = 7, height = 6)
sub_box
dev.off()

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

# Vegetation is our main interest so primvar. Because we also have negative lags, we can see influence of other variables on TVeg and of TVeg on the other variables. Strictly spoken it is not influence but correlation or 'synchronicity in the timeseries'
maximum_lag <- 7
w3ra_cor <- cor_norm_pool(fullset = Set2_p25deg, primvar = "TVeg", model = "W3RA", maxlag = maximum_lag)
htessel_cor <- cor_norm_pool(fullset = Set2_p25deg, primvar = "TVeg", model = "H-TESSEL", maxlag = maximum_lag)
pcr_cor <- cor_norm_pool(fullset = Set2_p25deg, primvar = "TVeg", model = "PCR-GLOBWB", maxlag = maximum_lag) #  
orch_cor <- cor_norm_pool(fullset = Set2_p25deg, primvar = "TVeg", model = "ORCHIDEE", maxlag = maximum_lag) # 

# Prepare a plot of correlations with TVeg. X-axis are the different lags.
all_cor <- list("W3RA" = w3ra_cor, "H-TESSEL" = htessel_cor, "PCR-GLOBWB" = pcr_cor, "ORCHIDEE" = orch_cor)
intermediate <- lapply(X= seq_along(all_cor), FUN = function(x) {
  var_interest <- all_cor[[x]][(1:(maximum_lag*2+1)),-(1:(maximum_lag*2+1))]
  var_interest$TVeg_lag <- ((-maximum_lag):maximum_lag)
  var_interest$Model <- names(all_cor)[x]
  return(gather(data = var_interest, key = "var", value = "spearmancor", -TVeg_lag, -Model))
})

plot_cor <- ggplot(data = do.call(what = rbind, args = intermediate), aes(x = TVeg_lag, y = spearmancor)) + geom_line(aes(color = Model)) + facet_wrap(~var) + labs(x = "Shift in modeled transpiration timeseries [months]", y = "Correlation between transpiration and other variable ")
pdf(file = "./An2_TVeglag_cor.pdf", width = 9.5, height = 5)
plot_cor
dev.off()

# == 2.3 == Between variable (lagged) correlation: spatially distributed

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
