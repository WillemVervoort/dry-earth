#!/usr/bin/env Rscript

# Runs in the DATA directory and creates a dataframe that for special cellnrs lists its classification values. Multiple classification levels are possible. For example cellnr's in the cooper creek area are part of the northern transect. Two tables are created because the cellnrs obviously differ with resolution (899 at 0.25 degree and 22475 at 0.05 degree). These can later be loaded and joined to the analysis datasets by cellnr to group the cellnr's according this classification.

require(dplyr)

# Bounding boxes in which the groundwater observations should fall
st_mac <- list(west = 147.3, north = -31, east = 147.75, south = -31.4)
st_edar <- list(west = 144.6, north = -31, east = 146, south = -31.75)
st_dar <- list(west = 143.5, north = -31, east = 144.25, south = -31.5)
st_wdar <- list(west = 141, north = -31, east = 142, south = -31.75)

nt_war <- list(west = 145.8, north = -26, east = 146.5, south = -27)
nt_dry <- list(west = 143.3, north = -26.5, east = 144.3, south = -27.25)
nt_cc <- list(west = 141.7, north = -26, east = 142.2, south = -26.4) # With this small box we select only areas in the ephemeral system creating a more unified hydrological selection than when also dry neighbouring cells are selected. The downside is that there are no groundwater observations.

# Function that creates the cellnrs lists with the bounding boxes above, then the hierarchy of lists is transformed to the table format. Calls upon two other functions below.
table_at_resolution <- function(resolutionset) {
  
  # Gather cellnrs at 0.05 degree for each of the transects. The southern one is dominated by summer rainfall
  indif_rain <- list(macquarie = give_cellnr(set = resolutionset, boxlist = st_mac),
                     eastofdarling= give_cellnr(set = resolutionset, boxlist = st_edar),
                     darling= give_cellnr(set = resolutionset, boxlist = st_dar),
                     westofdarling= give_cellnr(set = resolutionset, boxlist = st_wdar))
  
  sum_rain <- list(warrego = give_cellnr(set = resolutionset, boxlist = nt_war),
                   dry = give_cellnr(set = resolutionset, boxlist = nt_dry),
                   coopercreek = give_cellnr(set = resolutionset, boxlist = nt_cc))
  
  hierarchy <- list(indif_rain = indif_rain, sum_rain = sum_rain)
  
  # Go to the table format with function below and return that
  return(create_references(list_lists = hierarchy))
}

# Function that extracts the respective cellnrs, for a single dataset and a single bounding box
give_cellnr <- function(set, boxlist) {
  set %>% filter( x >= boxlist$west & x <= boxlist$east & y <= boxlist$north & y >= boxlist$south) %>% pull(cellnr) %>% unique() # Y coordinates are negative values.
}

# This function can have as input a list of lists of lists etc.. So multiple named levels of classification are possible
create_references <- function(list_lists) {
  #print(str(list_lists))
  stopifnot(all(!duplicated(unlist(list_lists, recursive = TRUE, use.names = FALSE)))) # One cell cannot belong to multiple classifications.
  
  # Recursive function that will get the names from each level and paste them together and outputs a vector of strings.
  getnames <- function(l, name = "") {
    if (!is.list(l)) {
      return(rep(name, length.out = length(l))) # The repetition is to make sure that all cellnumbers found at the deepest level recieve the same set of names (they are classified similarly)
    } else {
      namelist <- lapply(X = 1:length(l), FUN = function(x) getnames(l = l[[x]], name = paste(name, names(l)[x], sep = "."))) # This is the difficult part. Note the index x for looping over the entries of each list, extracting a name to supply as an argument, and then diving into the (possible) list beneath by doing l[[x]] 
      return(unlist(namelist, use.names = FALSE))
    }
  }
  
  gathernames <- strsplit(x = getnames(list_lists), split = ".", fixed = TRUE) # Calls the function and later converts the pasted strings to multiple vectors of separate ones. First string entry is empty, because of name = ""
  maxlength <- max(sapply(X = gathernames, FUN = length)) # Maxdepth is -1 because of empty
  namecols <- as.data.frame(t(sapply(X = gathernames, FUN = function(x) x[2:maxlength]))) # Here we convert the list output to columns. Discarding the first empty string, and automatically putting NA when level is not present.
  references <- cbind(unlist(list_lists, recursive = TRUE, use.names = FALSE), namecols) # The recursive unlist() gets the cellnrs, and we join them to the name columns
  names(references) <- c("cellnr", paste0("level", 1:ncol(namecols)))
  return(references)
}

# Load the 0.05 and 0.25 degree datasets created in the previous pre-processing step, call the functions and write the tables.
load("./Set1_p05deg.RData")
spatial_classification_p05deg <- table_at_resolution(resolutionset = Set1_p05deg) # Check the group-sizes: 
# as.tbl(spatial_classification_p05deg) %>% group_by(level1, level2) %>% summarise(n = n())
write.csv(x = spatial_classification_p05deg, file = "./spatial_classification_p05deg.csv", row.names = FALSE)

load("./Set2_p25deg.RData")
spatial_classification_p25deg <- table_at_resolution(resolutionset = Set2_p25deg)
write.csv(x = spatial_classification_p25deg, file = "./spatial_classification_p25deg.csv", row.names = FALSE)



