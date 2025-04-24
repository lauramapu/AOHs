
### 1. process ESA CCI and translate
# crop/mask with tropical forest
# from this step the process is per year (ESA-CCI)

rm(list=ls())

library(terra)
library(dplyr)
library(data.table)
library(readxl)
library(purrr)

forest <- vect('Spatial_Data/Tropical_Forest/tropicalmask.shp')

# load translation csv 
habitats <- read.csv('Habitats/translation_by_lumbierres.csv')
is <- habitats$Value
become <- habitats$code_lumbierres

# load srtm
srtm <- rast('Spatial_Data/SRTM90/strm_300m_trop.tif')

# download all needed ESA files
# in the following vector just write down the years you want
years <- c(1995, 2000, 2005, 2010, 2015)
# set download directory
download_dir <- 'Spatial_Data/ESA-LC'
# iterate though these years
for (i in years) {
  year <- years[i]
  ftp_url <- paste0('ftp://geo10.elie.ucl.ac.be/CCI/LandCover/byYear/ESACCI-LC-L4-LCCS-Map-300m-P1Y-', i, '-v2.0.7.tif')
  dest_file <- paste0(download_dir, 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-', i, '-v2.0.7.tif')
  download.file(ftp_url, destfile=dest_file, mode='wb')
}

# load all esa files (even if it is one or multiple)
esa_files <- list.files(
  path = 'Spatial_Data/ESA-LC/',
  pattern = 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-.*-v2\\.0\\.7\\.tif',
  full.names = TRUE
) %>%
  lapply(rast)

# reproyect mask to match esa files
forest <- project(forest, esa_files[[1]]) 

output_dir <- 'Spatial_Data/ESA-LC/processed/'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_final <- 'Spatial_Data/AOHs/baselayers/'
if (!dir.exists(output_final)) dir.create(output_final, recursive = TRUE)

# process each esa file
imap(esa_files, function(esa, i) { # use imap instead of lapply to get index
  
  # extract year from filename
  year <- regmatches(basename(sources(esa)), regexpr('\\d{4}', basename(sources(esa))))
  
  # crop and mask
  x <- esa %>% 
    crop(forest) %>% 
    mask(forest)
  
  # save
  out_name <- paste0(output_dir, 'esa_', year, '_trop.tif')
  writeRaster(x, out_name, overwrite = TRUE)
  
  # delete original
  file.remove(sources(esa))
  
  # in first iteration resample SRTM to ESA and round
  if (i==1) { 
    srtm <<- srtm %>% # <<- to keep the change in the global env
      resample(x, method='bilinear') %>%
      app(fun = function(x) {round(x/10)}) # round value to /10
  }
  else {}
  
  # reclass land uses to habitats 
  y <- subst(x, is, become)
  # y <- classify(ESA, cbind(is, become)) # another solution
  
  # multiply per 1000
  y <- app(y, fun = function(i) {i*1000})
  
  out_name <- paste0(output_dir, 'esa_', year, '_reclass.tif')
  writeRaster(y, out_name, overwrite=T)
  
  # and sum to srtm
  habitat_elev <- srtm + y

  # the resulting tif has pixels coded so that first 3 numbers are habitat (from 1 to 150)
  # and last 3 numbers are elev/10 (should range from 0 to 700 or something like that because top elevation is ~7000m)
  
  out_name <- paste0(output_final, 'baselayer_', year, '.tif')
  writeRaster(habitat_elev, out_name, overwrite=T)
  
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

mammals <- vect('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp')

base_files <- list.files('Spatial_Data/AOHs/baselayers', 
                         pattern='.tif',
                         full.names=T)

# nested loop in which first we select a year for the base layer
# and then extract distributions for that year
for (year in base_files) {
  
  base <- base_files[year]
  base <- base %>%
    rast() %>%
    project(mammals, method='near')
  
  for (i in seq_along(mammals)) {
    mammal <- mammals[i, ]
    output_file <- paste0('Spatial_Data/AOHs/', mammal$sci_name, '.tif')
    
    # Skip if the file already exists
    if (file.exists(output_file)) {
      message("Skipping ", mammal$sci_name, " (already processed)")
      next  # Jump to the next iteration
    }
    
    # Proceed with processing if file doesn't exist
    r <- base %>%
      crop(mammal) %>%
      mask(mammal)
    
    habitat_codes <- as.numeric(sub('_.*', '', hab_pref[[mammal$sci_name]]$Habitat_Code))
    elevation_range <- elev_range[[mammal$sci_name]]
    
    # 1) Compute elevation thresholds (in raster units)
    lo_e <- elevation_range$Lower_Elevation_Limit / 10
    hi_e <- elevation_range$Upper_Elevation_Limit / 10
    
    # if any is NA replace with zero (no elevation range so we assume all values are suitable)
    if (any(is.na(c(lo_e, hi_e)))) {
      lo_e <- 0
      hi_e <- 0
    }
    
    # 2) Initialize logical raster (all FALSE)
    cond <- r
    cond[] <- FALSE
    
    # 3) For each habitat code, set TRUE where r ∈ [code*1000 + lo_e, code*1000 + hi_e]
    for (code in habitat_codes) {
      minv <- code * 1000 + lo_e
      maxv <- code * 1000 + hi_e
      cond <- cond | ((r >= minv) & (r <= maxv))
    }
    
    # 4) Convert logical → binary (1/NA)
    binary_mask <- ifel(cond, 1, NA)
    names(binary_mask) <- mammal$sci_name
    
    # 5) Write output
    writeRaster(binary_mask, output_file, overwrite = TRUE)
    message("Processed and saved: ", mammal$sci_name)
  }
}
# 575 already done, starting at 17:05, stopping at 8am the next day ~1700
# # starting again at 8.40, 