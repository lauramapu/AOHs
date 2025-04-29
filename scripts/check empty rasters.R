
# from the process of more than 3000 species for the year 2010 we got several ones with no data
# so we should check if they are empty because there is actually no suitable habitat for those species in that year
# and while following the translation by lumbierres
# or if we are having any kind of issue with the code
# so

library(terra)
library(dplyr)
library(mapview)

# set path
dir_path <- 'Spatial_Data/AOHs'  # Replace with your actual directory path

# we found that empty rasters are all those that weight 1989kb or less
# that is, ordering by ascending size, from Congosorex polli

# get file information
files_info <- file.info(list.files(dir_path, recursive=TRUE, full.names = TRUE))

# filter files smaller than 1989 KB
empty_files <- files_info[files_info$size <= 1989, ]

# load a random one
random <- rast(rownames(empty_files)[21])
plot(random) # empty indeed
species_name <- random@pntr$names

# so we need to load the original distribution
mammals <- vect('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp')
mammal <- mammals[mammals$sci_name==species_name,]
mapview(mammal)

hab_pref <- readRDS('Habitats/mammal_habitat_preferences.rds')[[species_name]]
elev_range <- readRDS('Habitats/mammal_elevation_ranges.rds')[[species_name]]

base_files <- list.files('Spatial_Data/AOHs/baselayers', 
                         pattern='.tif',
                         full.names=T)
base <- rast(base_files[3])

# we can plot this with all pixels because it takes too much so we're just cropping and extracting pixel ranges
dist <- base %>%
  crop(mammal) %>%
  mask(mammal)
dist <- classify(dist, cbind(NA, 0)) # NAs complicate computations
unique(values(dist))
# here we can see several pixel values with the suitable habitat according to IUCN
# here it seems our problem is because artificial/terrestrial habitat is coded with '14' 
# and our codes are 141 and 144
# so we should check all habitat codes in the iucn search so that we have no discrepancies

hab_pref <- readRDS('Habitats/mammal_habitat_preferences.rds')
elev_range <- readRDS('Habitats/mammal_elevation_ranges.rds')

# Combine all data frames in the list
combined_df <- do.call(rbind, hab_pref)

# Get unique combinations
unique_habitats <- unique(combined_df[c("Habitat_Code", "Description")])

# Sort by Habitat_Code
unique_habitats <- unique_habitats[order(unique_habitats$Habitat_Code), ]
#View(unique_habitats)
# we have both 14_1 and 14 so we need to fix this somehow when doing the aoh
# check our habitat codes
habitats <- read.csv('Habitats/translation_by_lumbierres.csv')
habitats$code_lumbierres
# this would be the only case in which we need to do something because it is the only habitat for which we have subgroups

# redo the workflow we had adding an ifelse

# get habitat codes and elevation range from the current species
habitat_codes <- hab_pref[[mammal$sci_name]]$Habitat_Code
elevation_range <- elev_range[[mammal$sci_name]]

# compute elevation thresholds (/10)
lo_e <- elevation_range$Lower_Elevation_Limit / 10
hi_e <- elevation_range$Upper_Elevation_Limit / 10

# if any is NA replace with zero (no elevation range so we assume all values are suitable)
if (any(is.na(c(lo_e, hi_e)))) {
  lo_e <- 0
  hi_e <- 999
}

# initialize logical raster (all FALSE)
cond <- dist
cond[] <- FALSE

for (code_raw in habitat_codes) {
  
  # Handle special cases for codes starting with '14'
  if (startsWith(as.character(code_raw), '14')) {
    
    # Convert to proper numeric code based on suffix
    if (code_raw == '14') {
      code_vec <- 140:149  # Expand 14 to 140–149
    } else if (code_raw == '14_1' || code_raw == '14_2') {
      code_vec <- 141
    } else if (code_raw == '14_3' || code_raw == '14_6') {
      code_vec <- 143
    } else if (code_raw == '14_4' || code_raw == '14_5') {
      code_vec <- 144
    } else {
      warning(paste("Unknown 14_x code:", code_raw))
      next
    }
    
  } else {
    # For all other codes: just erase subgroup and convert to numeric
    code_vec <- as.numeric(gsub("_.*", "", code_raw))
  }
  
  # Apply condition logic for one or multiple habitat codes
  for (code in code_vec) {
    minv <- code * 1000 + lo_e
    maxv <- code * 1000 + hi_e
    cond <- cond | (dist >= minv) & (dist <= maxv)
  }
}

# convert trues to ones
binary_mask <- ifel(cond, 1, NA)
names(binary_mask) <- mammal$sci_name

# FINALLY DONNEEEEEEEEE
# we had three problems here
# 1. we have both '14' and '14_x' so we were not extracting any of these because of the subgroups
# 2. also because we had 14*1000 instead of 140*1000
# 3. and then because when no elevation is provided we consider lo_e and hi_e as zero that returned same value for
# minv and maxv so we did not have a range of elevation values so no values were taken

# all these changes already implemented in scripts 02