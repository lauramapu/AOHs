
# script 01 and 02 merged just for birds
# we can skip several computations already done

rm(list=ls())

library(terra)
library(dplyr)
library(data.table)
library(readxl)
library(rredlist)
library(iucnredlist) # devtools::install_github('IUCN-UK/iucnredlist')
library(stringr)
library(tictoc)
library(purrr)

# save log to file with same name as script
sink(paste0('log-', rstudioapi::getSourceEditorContext()$path %>% basename(), '.txt'), split=T, append=F)

store_slurm <- '/mnt/netapp1/Store_CSIC/home/csic/byc/abl/lamapu/AOHs/' # replace

forest <- vect('Spatial_Data/Tropical_Forest/tropicalmask.shp')

### 1. crop BirdLife distributions with tropical forest mask and save
# this step extracted from script 02 by iago

#Import bird distribution ranges from IUCN/Bird_Life
bird_ranges <- terra::vect(paste0(store_slurm, 'Spatial_Data/BirdLife_Range_Maps_Birds/BOTW.gdb')) # All_Species layer
# import second layer to extract categories
iucn_cats <- as.data.frame(terra::vect(paste0(store_slurm, 'Spatial_Data/BirdLife_Range_Maps_Birds/BOTW.gdb'),
                                       layer = 'BirdLifeTaxonomicChecklist_2022')) %>%
  rename(sisid=SISID)

# spatvector to df to operate better
bird_df <- as.data.frame(bird_ranges, geom = 'WKT')
bird_df <- left_join(bird_df, iucn_cats[, c('sisid', 'F2022_IUCN_RedList_Category')], by = 'sisid')
rm(iucn_cats)

# Project forest mask to ranges projection
forest <- terra::project(forest, crs(bird_ranges))

# here I select the ones used in Lumbierres et al 2022:

# 'we selected range polygons with extant and probably extant presence; native, reintroduced,
# and assisted colonization origin; and resident seasonality for non-migratory species
# (all birds and the 8,979 non-migratory birds). For migratory birds (1,816 species), we kept
# separated the ranges for breeding (1,446 species), non-breeding (1,550 species) and a
# combination of resident and uncertain (1,290 species) seasonality. [...]
# For 18 bird and 22 bird species categorized as Critical Endangered, there were no presence
# polygons coded as extant or probably extant. To assist the conservation of these species,
# we produced AOH maps using the possibly extinct polygon for these taxa.'

# codes of attributes are available in Spatial_Data/IUCN_Standard_attributes_for_spatial_data_v1.20_2024.xlsx
# our interest codes are the following:

# $presence -> extant == 1; Probably extant == 2; Possibly Extinct == 4
# $origin -> native == 1; reintroduced == 2, Assisted colonisation == 6
# $seasonal -> resident == 1; 2 == breeding; 3 == non-breeding; passant == 4; uncertain == 5

# first we need to separate migratory from non-migratory birds

# we assume all species that contain polygons with breeding, non-breeding or passant seasonality are migratory
# if they have no polygons with such categories (only resident or uncertain) we consider non-migratory
non_migratory <- bird_df %>%
  group_by(sci_name) %>%
  summarise(is_migratory = any(seasonal %in% c(2,3,4))) %>%
  filter(!is_migratory) %>%
  pull(sci_name)

migratory <- bird_df %>%
  group_by(sci_name) %>%
  summarise(is_migratory = any(seasonal %in% c(2,3,4))) %>%
  filter(is_migratory) %>%
  pull(sci_name)

# extract CR species
critical_species <- unique(bird_df$sci_name[bird_df$F2022_IUCN_RedList_Category == 'CR'])

# find CR species with no polygons $presence == 1 or 2
cr_no_presence_12 <- bird_df %>%
  filter(sci_name %in% critical_species) %>% # filter critical
  group_by(sci_name) %>%
  summarise(has_presence_12 = any(presence %in% c(1, 2))) %>%
  filter(!has_presence_12) %>%
  pull(sci_name)

rm(bird_df); gc()

# non-migratory

# set conditions
# first for the original criteria
condition_original <- 
  (bird_ranges$sci_name %in% non_migratory) &
  (bird_ranges$presence %in% c(1, 2)) &
  (bird_ranges$origin %in% c(1, 2, 6)) &
  (bird_ranges$seasonal == 1)

# add presence == 4 polygons for CR species lacking of presence == 1 or 2
condition_cr_presence4 <- 
  (bird_ranges$sci_name %in% non_migratory) &
  (bird_ranges$sci_name %in% cr_no_presence_12) & 
  (bird_ranges$presence == 4)

# combine conditions
final_condition <- condition_original | condition_cr_presence4

# subset
non_migratory <- bird_ranges[final_condition, ]

# function to aggregate by species, intersect with tropical mask and write

final_processing <- function(ranges, range_name) {
  ranges <- aggregate(ranges, by = 'sisid', dissolve = TRUE)
  ranges <- ranges[, c('sisid','sci_name')]
  names(ranges) <- c('IUCN_ID','IUCN_Species')
  tropical_ranges <- terra::intersect(ranges, forest)
  tropical_ranges <-terra::simplifyGeom(tropical_ranges, tolerance = 0.01) #keep tolerance super low
  terra::writeVector(tropical_ranges,
                     # no format so it writes a single file, we read the same way as any spatial file
                     paste0(store_slurm, 'Spatial_Data/BirdLife_Range_Maps_Birds/', range_name),
                     overwrite = TRUE)
}

final_processing(non_migratory, 'nonmigratory')
rm(non_migratory, condition_original, condition_cr_presence4, final_condition); gc()

# migratory

# breeding polygons
# the original criteria
condition_breeding <- 
  (bird_ranges$sci_name %in% migratory) &
  (bird_ranges$presence %in% c(1, 2)) &
  (bird_ranges$origin %in% c(1, 2, 6)) &
  (bird_ranges$seasonal == 2)
# add presence == 4 polygons for CR species lacking of presence == 1 or 2 
condition_cr_presence4 <- 
  (bird_ranges$sci_name %in% migratory) &
  (bird_ranges$sci_name %in% cr_no_presence_12) & 
  (bird_ranges$presence == 4) &
  (bird_ranges$seasonal == 2)
# combine conditions
breeding <- condition_breeding | condition_cr_presence4
# subset and write
breeding <- bird_ranges[breeding, ]
rm(condition_breeding, condition_cr_presence4); gc()
final_processing(breeding, 'breeding')
rm(breeding); gc()

# non-breeding polygons
condition_nonbreeding <- 
  (bird_ranges$sci_name %in% migratory) &
  (bird_ranges$presence %in% c(1, 2)) &
  (bird_ranges$origin %in% c(1, 2, 6)) &
  (bird_ranges$seasonal == 3)
condition_cr_presence4 <- 
  (bird_ranges$sci_name %in% migratory) &
  (bird_ranges$sci_name %in% cr_no_presence_12) & 
  (bird_ranges$presence == 4) &
  (bird_ranges$seasonal == 3)
nonbreeding <- condition_nonbreeding | condition_cr_presence4
nonbreeding <- bird_ranges[nonbreeding, ]
rm(condition_nonbreeding, condition_cr_presence4); gc()
final_processing(nonbreeding, 'nonbreeding')
rm(nonbreeding); gc()

# resident or uncertain polygons
condition_resident <- 
  (bird_ranges$sci_name %in% migratory) &
  (bird_ranges$presence %in% c(1, 2)) &
  (bird_ranges$origin %in% c(1, 2, 6)) &
  (bird_ranges$seasonal %in% c(1, 5))
condition_cr_presence4 <- 
  (bird_ranges$sci_name %in% migratory) &
  (bird_ranges$sci_name %in% cr_no_presence_12) & 
  (bird_ranges$presence == 4) &
  (bird_ranges$seasonal %in% c(1, 5))
resident_uncertain <- condition_resident | condition_cr_presence4
resident_uncertain <- bird_ranges[resident_uncertain, ]
rm(condition_resident, condition_cr_presence4); gc()
final_processing(resident_uncertain, 'resident_uncertain')
rm(resident_uncertain); gc()

rm(list=ls()); gc()

### 2. extract IUCN information on habitat and elevation range for each species

allbirds <- list.files(paste0(store_slurm, 'Spatial_Data/BirdLife_Range_Maps_Birds'),
                       pattern='^[^.]*$',
                       full.names=T) %>%
  lapply(vect) %>%
  lapply(as.data.frame)

allspecies <- data.frame()
for (i in seq_along(allbirds)) {
  ranges <- allbirds[[i]]
  species <- ranges$IUCN_Species
  ids <- ranges$IUCN_ID
  allspecies <- rbind(allspecies, data.frame(sci_name=species, sisid=ids))
}
allspecies <- distinct(allspecies) # erase duplicates

# Initialize habitat preferences list
bird_habitat_preferences <- list()
bird_elevation_ranges <- list()

# initialize token and api
token <- 'D98jKWcys1KfrHVAngXsno85KWjNw6s2qyyt'  # better if you get your own token
api <- init_api(token) 

tic()
for (i in 1:nrow(allspecies)) {
  
  cat(i, '/', nrow(allspecies),'\n')
  
  species_id <- allspecies[i, 'sisid']
  species_name <- as.character(allspecies[i, 'sci_name'])
  
  # **Check if species is already in the habitat preferences list**
  if (species_name %in% names(bird_habitat_preferences)) {
    next  # Skip to the next species
  }
  
  # Obtain species assessments
  assessment_raw <- iucnredlist::assessments_by_sis_id(api, species_id)
  assessment_raw <- assessment_raw %>% 
    dplyr::filter(year_published == max(year_published)) %>%
    dplyr::filter(latest == TRUE)
  
  # If no data, assume all habitats are suitable (Gallego-Zamorano assumption)
  if (is.null(assessment_raw$assessment_id)) next
  
  # Retrieve habitat data
  a_data <- assessment_data_many(api, assessment_raw$assessment_id, wait_time = 0.5)
  habitats <- extract_element(a_data, 'habitats')
  
  # Filter suitable habitats
  if (!is.null(habitats) && nrow(habitats) > 0) {
    habitats <- habitats %>%
      filter(suitability == 'Suitable') %>%
      select('Description' = description, 'Habitat_Code' = code, 'Season' = season) %>%
      distinct()
    
    # Store result with species name
    bird_habitat_preferences[[species_name]] <- habitats
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
  
  bird_elevation_ranges[[species_name]] <- elevation_data
}
toc() # 40425.73 sec elapsed

saveRDS(bird_habitat_preferences, 'Habitats/bird_habitat_preferences.rds')
saveRDS(bird_elevation_ranges, 'Habitats/bird_elevation_ranges.rds')

### 3. finally, generate AoH for each species using the habitat-elevation tif

# this process can't be parallelized in local because not enough RAM (some iterations take more than 50gb)
# if you want to send this to cesga and parallelize, you must wrap/unwrap spatraster/spatvector objects
# outside and inside the loop because they are non-exportable
# (see https://future.futureverse.org/articles/future-4-non-exportable-objects.html for more info)

hab_pref <- readRDS('Habitats/bird_habitat_preferences.rds')
elev_range <- readRDS('Habitats/bird_elevation_ranges.rds')

base_files <- list.files(paste0(store_slurm, 'Spatial_Data/AOHs/baselayers'), 
                         pattern='.tif',
                         full.names=T)

birds_files <- list.files(paste0(store_slurm, 'Spatial_Data/BirdLife_Range_Maps_Birds'),
                          pattern='^[^.]*$', # match when do not contain dots
                          full.names=T)

translation <- read.csv('Habitats/translation_by_lumbierres.csv')

# function to convert habitat codes (e.g., '1_9', '14_4') to translation$code format (e.g., 'H1', 'H14.1')
convert_habitat_code <- function(code_raw) {
  # split code into main and subcomponents (e.g., '1_9' → c('1', '9'))
  parts <- strsplit(as.character(code_raw), '_')[[1]]
  main <- parts[1]
  
  # handle special case for '14' (has subgroups in translation)
  if (main == '14') {
    if (length(parts) == 1) {
      # ff code is just '14', return all 14.x codes
      return(paste0('H14.', 1:6))
    } else {
      # if code is '14_x', map to specific subgroup
      subgroup <- parts[2]
      if (subgroup %in% c('1', '2')) return('H14.1')
      if (subgroup %in% c('3', '6')) return('H14.3')
      if (subgroup %in% c('4', '5')) return('H14.4')
      warning('Unknown 14 subgroup: ', code_raw)
      return(NA)
    }
  } else {
    # for all other codes, use main part (e.g., '1_9' → 'H1')
    return(paste0('H', main))
  }
}

# nested loop in which first we select a year for the base layer
# then a group of birds (non-migratory, migratory breeding, migratory non-breeding, migratory resident-uncertain)
# and then extract distributions for that year and group

for (i in seq_along(base_files)) {
  
  # get year
  year <- regmatches(base_files[i], regexpr('\\d{4}', base_files[i]))
  
  # read and project base layer
  base <- rast(base_files[i])
  
  for (j in seq_along(birds_files)) {
    
    # get type of bird
    type <- basename(birds_files[j])
    
    birds <- vect(birds_files[j])
    
    # generate output path per year and type of bird
    output_dir <- paste0(store_slurm, 'Spatial_Data/AOHs/birds/', type, '/', year, '/')
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    for (k in seq_along(birds)) {
      
      # get species
      bird <- birds[k, ]
      output_file <- paste0(output_dir, bird$IUCN_Species, '.tif')
      
      # skip if the species is already processed
      if (file.exists(output_file)) {
        message('Skipping ', bird$sci_name, ' (already processed)')
        next
      }
      
      # crop and mask base with the current distribution
      r <- base %>%
        crop(bird) %>%
        mask(bird)
      
      # get habitat codes and elevation range from the current species
      habitats <- hab_pref[[bird$IUCN_Species]]
      # with birds we must select the proper season
      if (type == 'breeding') {
        habitat_codes <- habitats[habitats$Season=='Breeding Season', ]$Habitat_Code
      } else if (type == 'nonbreeding') {
        habitat_codes <- habitats[habitats$Season=='Non-Breeding Season', ]$Habitat_Code
      } else {
        habitat_codes <- habitats[habitats$Season=='Resident', ]$Habitat_Code
      }
      if (length(habitat_codes) == 0 || is.na(habitat_codes) || habitat_codes == '') {
        cat('Species', bird$IUCN_Species, 'and', type, 'season', 'skipped because no suitable habitat found in IUCN API.\n')
        next
      }
      
      # # if there are no suitable habitats, we assume all habitats are suitable
      # if (is.null(habitat_codes)){
      #   habitat_codes <- c(1:8, '14_1', '14-2', '14_3', '14_4', '14_5', '14_6', 15)
      # }
      elevation_range <- elev_range[[bird$IUCN_Species]]
      
      # compute elevation thresholds (/10)
      lo_e <- elevation_range$Lower_Elevation_Limit / 10
      hi_e <- elevation_range$Upper_Elevation_Limit / 10
      
      # if any is NA replace with zero (no elevation range so we assume all values are suitable)
      if (any(is.na(c(lo_e, hi_e)))) {
        lo_e <- 0
        hi_e <- 999
      }
      
      # initialize logical raster (all FALSE)
      cond <- r
      cond[] <- FALSE
      
      for (code_raw in habitat_codes) {
        
        # convert raw code to translation$code format
        code_conv <- convert_habitat_code(code_raw)
        
        # get landuse codes for the hightest tertile only
        # (this can be modified if lower tertiles are needed)
        landuse_codes <- translation[translation$code==code_conv,'thr_high_code'] 
        if (length(landuse_codes) == 0 || is.na(landuse_codes) || landuse_codes == '') {
          cat('Species', bird$IUCN_Species, '+' , type, 'season + and habitat code', code_conv, 'skipped because habitat is not terrestrial.\n')
          next
        }
        # convert to numeric
        code_vec <- as.numeric(strsplit(landuse_codes, ';', fixed = TRUE)[[1]])
        
        # apply logic condition to each land use code
        # landuse * 1000 and between min-max elevation
        # operator | ensures this process is additive
        for (code in code_vec) {
          minv <- code * 1000 + lo_e
          maxv <- code * 1000 + hi_e
          cond <- cond | (r >= minv) & (r <= maxv)
        }
      }
      
      # convert trues to ones
      binary_mask <- ifel(cond, 1, NA)
      names(binary_mask) <- bird$sci_name
      
      # write output
      writeRaster(binary_mask, output_file, overwrite = TRUE)
      message('Processed and saved: ', bird$sci_name)
    }
  }
}

on.exit(sink())
