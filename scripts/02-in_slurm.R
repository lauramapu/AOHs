
### 1. process ESA CCI and translate
# crop/mask with tropical forest
# from this step the process is per year (ESA-CCI)

rm(list=ls())

library(terra)
library(dplyr)
library(data.table)
library(readxl)
library(purrr)

store_slurm <- '/mnt/netapp1/Store_CSIC/home/csic/byc/abl/lamapu/AOHs/' # replace

forest <- vect(paste0(store_slurm, 'Spatial_Data/Tropical_Forest/tropicalmask.shp'))

# load translation csv 
habitats <- read.csv(paste0(store_slurm, 'Habitats/translation_by_lumbierres.csv'))
is <- habitats$Value
become <- habitats$code_lumbierres

# load srtm
srtm <- rast(paste0(store_slurm, 'Spatial_Data/SRTM90/strm_300m_trop.tif'))

# download all needed ESA files
# in the following vector just write down the years you want
years <- c(1995, 2000, 2005, 2010, 2015)
# set download directory
download_dir <- paste0(store_slurm, 'Spatial_Data/ESA-LC')
if (!dir.exists(download_dir)) dir.create(download_dir, recursive = TRUE)
# iterate though these years
for (i in years) {
  year <- years[i]
  ftp_url <- paste0('ftp://geo10.elie.ucl.ac.be/CCI/LandCover/byYear/ESACCI-LC-L4-LCCS-Map-300m-P1Y-', i, '-v2.0.7.tif')
  dest_file <- paste0(download_dir, 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-', i, '-v2.0.7.tif')
  download.file(ftp_url, destfile=dest_file, mode='wb')
}

# load all esa files (even if it is one or multiple)
esa_files <- list.files(
  path = paste0(store_slurm, 'Spatial_Data/ESA-LC/'),
  pattern = 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-.*-v2\\.0\\.7\\.tif',
  full.names = TRUE
) %>%
  lapply(rast)

# reproyect mask to match esa files
forest <- project(forest, esa_files[[1]]) 

# output_dir <- 'Spatial_Data/ESA-LC/processed/'
# if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_final <- paste0(store_slurm, 'Spatial_Data/AOHs/baselayers/')
if (!dir.exists(output_final)) dir.create(output_final, recursive = TRUE)

# process each esa file
imap(esa_files, function(esa, i) { # use imap instead of lapply to get index
  
  # extract year from filename
  year <- regmatches(basename(sources(esa)), regexpr('\\d{4}', basename(sources(esa))))
  
  # Check if output file already exists
  out_name <- paste0(output_final, 'baselayer_', year, '.tif')
  if (file.exists(out_name)) {
    cat("Skipping year", year, "- output already exists:", out_name, "\n")
    return(NULL)  # Skip this iteration
  }
  
  # crop and mask
  x <- esa %>% 
    crop(forest) %>% 
    mask(forest)
  
  # delete original
  file.remove(sources(esa))
  
  # in first iteration resample SRTM to ESA and round
  if (i==1) { 
    srtm <<- srtm %>% # <<- to keep the change in the global env
      resample(x, method='bilinear') %>%
      app(fun = function(x) {round(x/10)}) # round value to /10
  }
  
  # reclass land uses to habitats 
  y <- subst(x, is, become)
  
  # multiply per 1000
  y <- app(y, fun = function(i) {i*1000})
  
  # and sum to srtm
  habitat_elev <- srtm + y
  
  writeRaster(habitat_elev, out_name, overwrite=TRUE)
  cat(out_name, ' done\n')
})

rm(srtm, esa_files)

### 2. finally, generate AoH for each species using the habitat-elevation tif

# this process can't be parallelized in local because not enough RAM (some iterations take more than 50gb)
# if you want to send this to cesga and parallelize, you must wrap/unwrap spatraster/spatvector objects
# outside and inside the loop because they are non-exportable
# (see https://future.futureverse.org/articles/future-4-non-exportable-objects.html for more info)

hab_pref <- readRDS('Habitats/mammal_habitat_preferences.rds')
elev_range <- readRDS('Habitats/mammal_elevation_ranges.rds')

mammals <- vect(paste0(store_slurm, 'Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp'))

base_files <- list.files(paste0(store_slurm, 'Spatial_Data/AOHs/baselayers'), 
                         pattern='.tif',
                         full.names=T)

# nested loop in which first we select a year for the base layer
# and then extract distributions for that year
for (i in seq_along(base_files)) {
  
  # get year
  year <- regmatches(base_files[i], regexpr('\\d{4}', base_files[i]))
  
  # generate output path per year
  output_dir <- paste0(store_slurm, 'Spatial_Data/AOHs/', year, '/')
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # read and project base layer
  base <- base_files[i]
  base <- rast(base) %>%
    project(mammals, method='near')
  
  for (j in seq_along(mammals)) {
    
    # get species
    mammal <- mammals[j, ]
    output_file <- paste0(store_slurm, 'Spatial_Data/AOHs/', year, '/', mammal$sci_name, '.tif')
    
    # skip if the species is already processed
    if (file.exists(output_file)) {
      message("Skipping ", mammal$sci_name, " (already processed)")
      next
    }
    
    # crop and mask base with the current distribution
    r <- base %>%
      crop(mammal) %>%
      mask(mammal)
    
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
    
    # write output
    writeRaster(binary_mask, output_file, overwrite = TRUE)
    message("Processed and saved: ", mammal$sci_name)
  }
}
