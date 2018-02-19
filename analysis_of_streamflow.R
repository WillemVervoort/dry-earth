#!/usr/bin/env Rscript

# inherits the RESULTS working directory from the bash script

require(tibble)
require(lubridate)
require(dplyr)
require(tidyr)
require(ggplot2)
require(maps)
require(RColorBrewer)
#require(gridExtra)

# ===========================
# Part 0 Load the data
# ===========================

csv_files <- list.files(path = "../DATA", pattern = "str_[[:alpha:]]+.csv", full.names = TRUE)
# Daily mean discharges in m3 s-1

store <- lapply(X = seq_along(csv_files), FUN = function(i) {
  data_f <- as.tibble(read.csv(csv_files[i], skip=9, stringsAsFactors = FALSE)) # Measurements
  
  # Some extra info
  station_name <- as.character(read.csv(csv_files[i], header = FALSE, nrows = 1, stringsAsFactors = FALSE)[[2]])
  region <- substr(x = csv_files[i], start = 13, stop = nchar(csv_files[i])-4) # Matches spatial classification
  lat <- read.csv(csv_files[i], header = FALSE, skip = 2, nrows = 1, stringsAsFactors = FALSE)[[2]]
  lon <- read.csv(csv_files[i], header = FALSE, skip = 3, nrows = 1, stringsAsFactors = FALSE)[[2]]
  
  cutoff <- quantile(data_f$Value,0.95, na.rm=TRUE) # For calculation of the flood frequency curve, discard the zero values

  data_f <- data_f %>% 
    #filter(Value > cutoff) %>%
    mutate(region = region, station_name = station_name, long = lon, lat = lat) %>%
    mutate(reference = t(as.data.frame(strsplit(X.Timestamp, split = "[+]")))[,1], add = t(as.data.frame(strsplit(X.Timestamp, split = "[+]")))[,2]) %>% # The plus sign does not work in regular date parsing, so we split it to the reference and the hours that need to be added
    mutate(Time = ymd_hms(reference) + hm(add)) %>% 
    select(-Interpolation.Type, -Quality.Code, -X.Timestamp, -reference, -add) #
})

plot_df <- do.call(rbind,store)

# ===========================
# Part 1 Streamflow plots
# ===========================

# normalisation of values for each station. Might be useful.
plot_df <- plot_df %>% group_by(region) %>% mutate(normalized = scale(Value))

plot_df$station_name <- factor(x = plot_df$station_name, levels = unique(plot_df$station_name), labels = c("Cooper Creek at Cullyamurra Water Hole", "Darling river at Louth", "Macquarie river at Pillicawarrina", "Warrego river at Cunnamulla Weir"))

plot_pal <- brewer.pal(n = 7, name = "Set1")[c(1,2,5,6)] # To match the colors used for these regions 

# Unmodified hydrographs
pl <- ggplot(plot_df,aes(x = Time, y = Value)) + 
  geom_line(aes(colour=region)) +
  xlab("Date") + ylab("Flow m3 s-1") +
  scale_color_manual(values = plot_pal)
pl

# Plot of normalized values. Not very useful, because where is the zero discharge?
flow_duration_norm <- ggplot(data = plot_df, aes(x = normalized)) + geom_step(aes(y=-..y..+1, color = region),stat="ecdf") + coord_flip() + labs(x = "Discharge [m3 s-1]", y = "Probability of exceedence") + ylim(0,0.3) + scale_color_manual(values = plot_pal)

# Plot of unmodified values but on a logarithmic scale. This nicely fits with the interpretation that warrego is the least under influence of surface water supply. And that macquarie has the most steady stream, even though the peaks are not extremely large.
flow_duration_log <- ggplot(data = plot_df, aes(x = Value)) + geom_step(aes(y=-..y..+1, color = region),stat="ecdf") + labs(x = "Discharge [m3 s-1]", y = "Probability of exceedence") + scale_x_log10() + scale_color_manual(values = plot_pal) + coord_flip() # the -y+1 is to transform non-exceedence (standard in cdf) to exceedence (standard in flow duration curves)
pdf(file = "./An1_flow_duration_log.pdf", width = 5, height = 3)
flow_duration_log
dev.off()

# ===========================
# Part 2 Map making
# ===========================

australia <- map_data(map = "world", region = "Australia")
studyarea <- data.frame(long = c(140,147.75,140,147.75), lat = c(-25,-25,-32.25,-32.25), order = c(1,2,4,3))
mapaus <-ggplot(data = australia, aes(x = long, y = lat)) + geom_polygon(aes(group = group), fill = NA, color = "black") + coord_fixed(xlim = c(110,155), ylim = c(-8,-45)) + geom_polygon(data = studyarea, aes(order = order), fill = NA, color = "grey50") + geom_point(data = plot_df %>% group_by(station_name) %>% summarise(long = long[1], lat = lat[1]), aes(color = station_name)) + scale_color_manual(values = plot_pal) + labs(x = "longitude", y = "latitude", color = "Discharge station") 
pdf(file = "./studyarea_and_stations.pdf", width = 7.5, height = 3.5)
mapaus
dev.off()

