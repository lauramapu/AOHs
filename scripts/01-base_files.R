
############
# workflow to filter IUCN distribution maps per elevation range and habitat correspondence for each species
############################

# there's a difference with the original work by Lumbierres et al (2022) because we need to use ESA-CCI
# to have all years in our study range while she uses CGLS-LC100
# also habitat translation is sometimes confusing

# if an error occurs when crop/mask (terra) telling something like too many values to write read this:
# https://github.com/rspatial/terra/issues/1686
# and this to do the clean installation:
# https://stackoverflow.com/questions/70925962/r-packages-raster-fail-to-upload-while-searching-for-terra-last-version

# in this script, we are generating all files needed for every year

### 1. crop/mask mammal ranges with tropical forest mask and extract wanted distributions
# great part of this code extracted from 02_IUCN_Ranges by Iago)

rm(list=ls())

library(terra)
library(dplyr)
library(data.table)
library(readxl)
library(rredlist)
library(iucnredlist) # devtools::install_github('IUCN-UK/iucnredlist')
library(stringr)
library(tictoc)

#Import pantropical forest zone shp
forest <- terra::vect('Spatial_Data/Tropical_Forest/TropicalForest.shp')
# dissolve to fix bad geometries and save
forest$id <- 1
forest <- aggregate(forest, by = 'id', dissolve = TRUE)
writeVector(forest, 'Spatial_Data/Tropical_Forest/tropicalmask.shp', overwrite=T)

#Import terrestrial mammals distribution ranges from IUCN.

# download file if needed from https://www.iucnredlist.org/resources/files/7864a4a2-f7ed-422f-8a43-7fc42e3883c9
# (you must fill form and place under the folder 'Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals)
mammal_ranges <- terra::vect('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/MAMMALS_TERRESTRIAL_ONLY.shp')

# here I select the ones used in Lumbierres et al 2022:

# "we selected range polygons with extant and probably extant presence; native, reintroduced,
# and assisted colonization origin; and resident seasonality for non-migratory species
# (all mammals and the 8,979 non-migratory birds) [...].
# For 18 mammal and 22 bird species categorized as Critical Endangered, there were no presence
# polygons coded as extant or probably extant. To assist the conservation of these species,
# we produced AOH maps using the possibly extinct polygon for these taxa."

# codes of attributes are available in Spatial_Data/IUCN_Standard_attributes_for_spatial_data_v1.20_2024.xlsx
# our interest codes are the following:

# $presence -> extant == 1; Probably extant == 2; Possibly Extinct == 4
# $origin -> native == 1; reintroduced == 2, Assisted colonisation == 6
# $seasonal -> resident == 1

# spatvector to df to operate better
mammal_df <- as.data.frame(mammal_ranges, geom = "WKT")

# extract CR species
critical_species <- unique(mammal_df$sci_name[mammal_df$category == "CR"])

# find CR species with no polygons $presence == 1 or 2
cr_no_presence_12 <- mammal_df %>%
  filter(sci_name %in% critical_species) %>%
  group_by(sci_name) %>%
  summarise(has_presence_12 = any(presence %in% c(1, 2))) %>%
  filter(!has_presence_12) %>%
  pull(sci_name)

# set conditions
# first for the original criteria
condition_original <- 
  (mammal_ranges$presence %in% c(1, 2)) &
  (mammal_ranges$origin %in% c(1, 2, 6)) &
  (mammal_ranges$seasonal == 1)

# add presence == 4 polygons for CR species lacking of presence == 1 or 2
condition_cr_presence4 <- 
  (mammal_ranges$sci_name %in% cr_no_presence_12) & 
  (mammal_ranges$presence == 4)

# combine conditions
final_condition <- condition_original | condition_cr_presence4

# subset
mammal_ranges <- mammal_ranges[final_condition, ]

# Combine multiple polygons of each species into one.
mammal_ranges <- mammal_ranges[, c('id_no','sci_name')]
mammal_ranges <- aggregate(mammal_ranges, by = 'id_no', dissolve = TRUE)

#Intersect mammal distribution range with tropical forest extent
forest <- project(forest, mammal_ranges)
tropical_ranges <- terra::intersect(mammal_ranges, forest)
rm(mammal_ranges)

# save single shp
terra::writeVector(tropical_ranges, 
                   'Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp',
                   overwrite = TRUE)

# #Save each species as single shp
# for (i in seq_along(unique(tropical_ranges$sci_name))) {
#   
#   sci_name <- unique(tropical_ranges$sci_name)[i]
#   
#   terra::writeVector(tropical_ranges[which(tropical_ranges$sci_name == sci_name), ], 
#                      paste0('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/', sci_name, '.shp'),
#                      overwrite = TRUE)
# }

rm(tropical_ranges)

### 2. extract habitat and elevation ranges information on each species with the IUCN API
# (extracted from script 04_Habitat_Preferences done by Iago, just for mammals)

# initialize lists to store habitats and elevation ranges
mammal_habitat_preferences <- list()
mammal_elevation_ranges <- list()

# initialize token and api
token <- 'D98jKWcys1KfrHVAngXsno85KWjNw6s2qyyt'  # better if you get your own token
api <- init_api(token) 

# import tropical terrestrial mammals created before
mammals <- vect('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp')
# our field of interest is 'id_no' which is an identifier for each species
species <- data.frame(id = unique(mammals$id_no), sci_name = unique(mammals$sci_name))

# this takes a lot of time so you can just unsilence this and load the files
mammal_habitat_preferences <- readRDS('Habitats/mammal_habitat_preferences.rds')
mammal_elevation_ranges <- readRDS('Habitats/mammal_elevation_ranges.rds')

tic()
for (i in 1:nrow(species)) {
  print(i)  
  
  # get id and name
  species_id <- species[i, 'id']
  species_name <- as.character(species[i, 'sci_name'])
  
  # verify if species is already processed
  if (species_name %in% names(mammal_habitat_preferences)) {next}
  
  # Extract most recent assessment 
  assessment_raw <- iucnredlist::assessments_by_sis_id(api, species_id)
  assessment_raw <- assessment_raw %>% 
    dplyr::filter(year_published == max(year_published)) %>%
    dplyr::filter(latest == TRUE)
  
  if (is.null(assessment_raw$assessment_id)) {
    print(paste('No data for species:', species_name, 'Skipping...'))
    next  
  }
  
  # Obtain assessments from species
  a_data <- assessment_data_many(api, assessment_raw$assessment_id, wait_time = 0.5)
  habitats <- extract_element(a_data, 'habitats')
  
  #filter suitable habitat of the species
  if (!is.null(habitats) && nrow(habitats) > 0) {
    habitats <- habitats %>% 
      filter(suitability == 'Suitable') %>% 
      select('Description' = description, 'Habitat_Code' = code) %>%  
      distinct() 
    
    # save results in list
    mammal_habitat_preferences[[species_name]] <- habitats
  }
  
  # extract elevation info if it exists
  elevation <- extract_element(a_data, 'supplementary_info')
  
  # just verify if there is info about elevation limits
  if (!is.null(elevation) && all(c('upper_elevation_limit', 'lower_elevation_limit') %in% colnames(elevation))) {
    elevation_data <- elevation %>% 
      select('Upper_Elevation_Limit' = upper_elevation_limit, 
             'Lower_Elevation_Limit' = lower_elevation_limit) %>% 
      mutate('IUCN_Species' = species_name) #add species name
  } else {
    # If there is no data, create the df with NA values
    elevation_data <- data.frame(
      'Upper_Elevation_Limit' = NA,
      'Lower_Elevation_Limit' = NA,
      'IUCN_Species' = species_name
    )
  }
  
  mammal_elevation_ranges[[species_name]] <- elevation_data
}
toc()  

# Save results 
saveRDS(mammal_habitat_preferences, 'Habitats/mammal_habitat_preferences.rds')
saveRDS(mammal_elevation_ranges, 'Habitats/mammal_elevation_ranges.rds')

### 3. process SRTM files
# SRTM was extracted and resampled using GEE (3 files)
# to get these files, run the script 'gee_srtm_global' in GEE and place them in the 'Spatial_Data/SRTM90' folder

srtm_files <- list.files('Spatial_Data/SRTM90', pattern='.tif$', full.names=TRUE)
srtm_list <- lapply(srtm_files, rast)
srtm_merged <- do.call(merge, srtm_list)

forest <- project(forest, srtm_merged)

srtm <- srtm_merged %>%
  crop(forest) %>%
  mask(forest)

writeRaster(srtm, 'Spatial_Data/SRTM90/strm_300m_trop.tif', overwrite=T)

### 4. translate habitats as described in IUCN to ESA-CCI land uses basing on Lumbierres et al 2021 method

# Lumbierres generates correlations between habitats and land uses and she provides thresholds to convert 

# first load correspondences between esa and iucn
# table from Lumbierres et al 2021 (Figure 3)
# this is basically a correlation table between ESA and IUCN classes
translation <- read_xlsx('Habitats/esa-habitats.xlsx') %>%
  mutate(across(3:13, ~ {
    x <- ifelse(.x == '-', NA, .x) # replace '-' with NA
    as.numeric(x) # convert to numeric
  }))

# ESA-CCI legend (PDF converted to excel)
cci_legend <- read_xlsx('Habitats/CCI-LC_Maps_Legend.xlsx') %>%
  # new column to match land use groups used in Lumbierres et al 2021
  mutate(lumbierres = c(
    NA,
    'Cropland, rainfed',
    'Cropland, rainfed: herbaceous cover',
    'Cropland, rainfed: tree or shrub cover',
    'Cropland irrigated or post-flooding',
    'Mosaic cropland (>50%) /natural vegetation (tree, shrub, herbaceous cover)(<50%)',
    'Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%)/cropland (<50%)',
    'Tree cover, broadleaved, evergreen, closed to open (>15%)', 
    'Tree cover, broadleaved, deciduous, closed to open (>15%)', # grouped
    'Tree cover, broadleaved, deciduous, closed to open (>15%)',
    'Tree cover, broadleaved, deciduous, closed to open (>15%)',
    'Tree cover, needleleaved, evergreen, closed to open (>15%)', # grouped
    'Tree cover, needleleaved, evergreen, closed to open (>15%)',
    'Tree cover, needleleaved, evergreen, closed to open (>15%)',
    'Tree cover, needleleaved, deciduous, closed to open (>15%)', # grouped
    'Tree cover, needleleaved, deciduous, closed to open (>15%)',
    'Tree cover, needleleaved, deciduous, closed to open (>15%)',
    'Tree cover, mixed leaf type (broadleaved and needleleaved)',
    'Mosaic tree and shrub (>50%) / herbaceous cover (<50%)',
    'Mosaic herbaceous cover (>50%) / tree and shrub (<50%)',
    'Shrubland', # grouped
    'Shrubland',
    'Shrubland',
    'Grassland',
    'Lichens and mosses',
    'Sparse vegetation (tree, shrub,herbaceous cover)(<15%)', # grouped
    'Sparse vegetation (tree, shrub,herbaceous cover)(<15%)',
    'Sparse vegetation (tree, shrub,herbaceous cover)(<15%)',
    'Sparse vegetation (tree, shrub,herbaceous cover)(<15%)',
    'Tree cover, flooded, fresh or brakish water',
    'Tree cover, flooded, saline water',
    'Shrub or herbaceous cover, flooded, fresh/ saline/brakish water',
    'Urban area',
    'Bare areas', # grouped
    'Bare areas',
    'Bare areas',
    'Water bodies', # grouped
    'Water bodies'
  ))

# iucn habitats with code and description (extracted from sup materials in Lumbierres et al 2021)
iucn <- tibble::tibble(
  code = c(
    'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H8', 
    'H14.1', 'H14.2', 'H14.3', 'H14.6', 'H14.4', 'H14.5', 'H15'
  ),
  name = c(
    'Forest',
    'Savanna',
    'Shrubland',
    'Grassland',
    'Wetlands',
    'Rocky Areas',
    'Desert',
    'Artificial arable and pasture lands: Arable Land', # 1
    'Artificial arable and pasture lands: Pastureland',# 2
    'Artificial degraded forest and plantation: Plantations', # 3
    'Artificial degraded forest and plantation: Degraded Forest', # 6
    'Artificial urban and rural gardens: Rural Gardens', # 4
    'Artificial urban and rural gardens: Urban Areas', # 5
    'Artificial Aquatic'
  ),
  description = c(
    'Forest consists of a continuous stand of trees and includes both forested areas (generally with a closed canopy) and wooded areas.',
    'Savannas are transitional between grasslands and forests. They are ecosystems dominated by a grass ground cover with an overstorey of widely spaced trees and shrubs.',
    'Also referred to as scrub, bushland and thicket.',
    'Native grasslands are comprised of grasses and broadleaved herbaceous plants, and are either without woody plants, or the latter are very sparsely distributed.',
    'Areas of marsh, fen, peatland or water, whether natural or artificial, permanent or temporary, with water that is static or flowing, fresh, brackish or salt, including areas of marine water the depth of which at low tide does not exceed six meters.',
    'Cliffs, mountain peaks, talus, feldmark.',
    'Consists of arid landscapes with a sparse plant cover, except in depressions where water accumulates. The sandy, stony, or rocky substrate contributes more to the appearance of the landscape than does the vegetation.',
    'Includes cereal fields, rice paddies, perennial crops, orchards and groves.',
    'Includes fertilized or re-seeded permanent grasslands, sometimes treated with selective herbicides, with very impoverished flora and fauna. Also includes secondary grasslands and wooded farmland.',
    'A plantation is an intentional planting of a crop, on a larger scale, usually for uses other than cereal production or pasture. The term is currently most often used for plantings of trees and shrubs. The term tends also to be used for plantings maintained on economic bases other than that of subsistence farming.',
    'Former subtropical or tropical forest that has been extensively cleared or impacted by human activities. Often there is some degree of regeneration or there are small fragments of forest remaining.',
    'Rural gardens are located in a rural setting, serving families whose main income comes from wage labor (rural or urban)... [description continues, see source]',
    'Usually metropolitan and commercial areas dominated by asphalt, concrete and roof. Includes buildings, lawns and parks.',
    'These are human-made wetland habitats.'
  ),
  habitat = c(
    'Forest',
    'Savanna',
    'Shrubland',
    'Grassland',
    'Wetlands',
    'Rocky areas',
    'Desert',
    'Artificial arable and pasture lands',
    'Artificial arable and pasture lands',
    'Artificial degraded forest and plantation',
    'Artificial degraded forest and plantation',
    'Artificial urban areas and rural gardens',
    'Artificial urban areas and rural gardens',
    'Artificial aquatic'
  )
)

# tertiles:
# low: 1.138:1.351
# mid: 1.362:1.712
# high: 1.743-13.720

# so for each tertile we are creating two columns (one for name and other for code)

# names
iucn$thr_low <- NA_character_
iucn$thr_mid <- NA_character_
iucn$thr_high <- NA_character_

# codes
iucn$thr_low_code <- NA_character_
iucn$thr_mid_code <- NA_character_
iucn$thr_high_code <- NA_character_

# define thresholds
thresholds <- list(
  low = c(1.138, 1.351),
  mid = c(1.362, 1.712),
  high = c(1.743, 13.720)
)

# iterate through each habitat in iucn
for (i in 1:nrow(iucn)) {
  habitat_name <- iucn$habitat[i]
  
  if (!habitat_name %in% colnames(translation)) next
  
  # initialize lists for names and codes
  low_classes <- mid_classes <- high_classes <- character()
  low_codes <- mid_codes <- high_codes <- character()
  
  for (j in 1:nrow(translation)) {
    thr_class <- translation$`IUCN habitat class Land-cover class`[j]
    value <- translation[[habitat_name]][j]
    
    if (!is.na(value)) {
      # get corresponding codes from cci_legend (may be multiple)
      codes <- cci_legend$Value[cci_legend$lumbierres == thr_class]
      codes_str <- paste(codes, collapse = "; ")  # combine multiple codes
      
      if (value >= thresholds$low[1] && value <= thresholds$low[2]) {
        low_classes <- c(low_classes, thr_class)
        low_codes <- c(low_codes, codes_str)
      } else if (value >= thresholds$mid[1] && value <= thresholds$mid[2]) {
        mid_classes <- c(mid_classes, thr_class)
        mid_codes <- c(mid_codes, codes_str)
      } else if (value >= thresholds$high[1] && value <= thresholds$high[2]) {
        high_classes <- c(high_classes, thr_class)
        high_codes <- c(high_codes, codes_str)
      }
    }
  }
  
  # assign names and codes
  if (length(low_classes) > 0) {
    iucn$thr_low[i] <- paste(low_classes, collapse = "; ")
    iucn$thr_low_code[i] <- paste(low_codes, collapse = "; ")
  }
  if (length(mid_classes) > 0) {
    iucn$thr_mid[i] <- paste(mid_classes, collapse = "; ")
    iucn$thr_mid_code[i] <- paste(mid_codes, collapse = "; ")
  }
  if (length(high_classes) > 0) {
    iucn$thr_high[i] <- paste(high_classes, collapse = "; ")
    iucn$thr_high_code[i] <- paste(high_codes, collapse = "; ")
  }
}

# somehow we get many NAs in the middle of the strings, so we need to clean

# function to remove NA from strings
clean_na_strings <- function(x) {
  if (is.na(x)) {
    return(NA_character_)  # keep original NA values
  }
  # split string, remove "NA", and recombine
  parts <- strsplit(x, "; ")[[1]]
  parts_clean <- parts[parts != "NA"]
  if (length(parts_clean) == 0) {
    return(NA_character_)  # return NA if all parts were "NA"
  } else {
    return(paste(parts_clean, collapse = "; "))
  }
}

# apply cleaning to all code columns
iucn$thr_low_code <- sapply(iucn$thr_low_code, clean_na_strings)
iucn$thr_mid_code <- sapply(iucn$thr_mid_code, clean_na_strings)
iucn$thr_high_code <- sapply(iucn$thr_high_code, clean_na_strings)

write.csv(iucn, 'Habitats/translation_by_lumbierres.csv', row.names=F)
