
### 1. process ESA CCI and translate
# crop/mask with tropical forest

rm(list=ls())

library(terra)
library(dplyr)
library(data.table)
library(readxl)
library(purrr)

forest <- vect('Spatial_Data/Tropical_Forest/tropicalmask.shp')

# load srtm
srtm <- rast('Spatial_Data/SRTM90/strm_300m_trop.tif')

# download all needed ESA files
# in the following vector just write down the years you want
years <- c(1995, 2000, 2005, 2010, 2015)

# set download directory
download_dir <- 'Spatial_Data/ESA-LC/'
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
  path = 'Spatial_Data/ESA-LC/',
  pattern = 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-.*-v2\\.0\\.7\\.tif',
  full.names = TRUE
) %>%
  lapply(rast)

# reproyect mask to match esa files
forest <- project(forest, esa_files[[1]]) 

# output_dir <- 'Spatial_Data/ESA-LC/processed/'
# if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_final <- 'Spatial_Data/AOHs/baselayers/'
if (!dir.exists(output_final)) dir.create(output_final, recursive = TRUE)

# process each esa file
imap(esa_files, function(esa, i) { # use imap instead of lapply to get index
  
  # extract year from filename
  year <- regmatches(basename(sources(esa)), regexpr('\\d{4}', basename(sources(esa))))
  
  # crop, mask and multiply per 1000
  x <- esa %>% 
    crop(forest) %>% 
    mask(forest) %>%
    app(fun = function(x) {x*1000})
  
  # # save
  # out_name <- paste0(output_dir, 'esa_', year, '_trop.tif')
  # writeRaster(x, out_name, overwrite = TRUE)
  
  # delete original
  file.remove(sources(esa))
  
  # in first iteration resample SRTM to ESA and round
  if (i==1) { 
    srtm <<- srtm %>% # <<- to keep the change in the global env
      resample(x, method='bilinear') %>%
      app(fun = function(x) {round(x/10)}) # round value to /10
  }
  else {}
  
  # and sum to srtm
  habitat_elev <- srtm + x

  # the resulting tif has pixels coded so that first 3 numbers are habitat (from 1 to 150)
  # and last 3 numbers are elev/10 (should range from 0 to 700 or something like that because top elevation is ~7000m)
  
  out_name <- paste0(output_final, 'baselayer_', year, '.tif')
  writeRaster(habitat_elev, out_name, overwrite=T)
  
  cat(out_name, ' done\n')
  
})

rm(srtm, esa_files)

### 2. finally, generate AoH for each species using the landuse-elevation tif

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

# check that this works
convert_habitat_code('1_9')    # Returns 'H1'
convert_habitat_code('14_4')   # Returns 'H14.4'
convert_habitat_code('14')     # Returns c('H14.1', 'H14.2', ..., 'H14.6')
convert_habitat_code('5')      # Returns 'H5'

# nested loop in which first we select a year for the base layer
# and then extract distributions for that year
for (i in seq_along(base_files)) {
  
  # get year
  year <- regmatches(base_files[i], regexpr('\\d{4}', base_files[i]))
  
  # generate output path per year
  output_dir <- paste0('Spatial_Data/AOHs/mammals/', year, '/')
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # read base layer
  base <- rast(base_files[i])

  for (j in seq_along(mammals)) {
    
    # get species
    mammal <- mammals[j, ]
    output_file <- paste0(output_dir, mammal$sci_name, '.tif')
    
    # skip if the species is already processed
    if (file.exists(output_file)) {
      message('Skipping ', mammal$sci_name, ' (already processed)')
      next
    }
    
    # crop and mask base with the current distribution
    r <- base %>%
      crop(mammal) %>%
      mask(mammal)
    
    # get habitat codes and elevation range from the current species
    habitat_codes <- hab_pref[[mammal$sci_name]]$Habitat_Code
    if (length(habitat_codes) == 0 || is.na(habitat_codes) || habitat_codes == '') {
      cat('Species', mammal$sci_name, 'skipped because no suitable habitat found in IUCN API.\n')
      next
    }
    # # if there are no suitable habitats, we assume all habitats are suitable
    # if (is.null(habitat_codes)){
    #   habitat_codes <- c(1:8, '14_1', '14-2', '14_3', '14_4', '14_5', '14_6', 15)
    # }
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
    cond <- r
    cond[] <- FALSE
    
    for (code_raw in habitat_codes) {
      
      # convert raw code to translation$code format
      code_conv <- convert_habitat_code(code_raw)
      
      # get landuse codes for the highest tertile only
      # (this can be modified if lower tertiles are needed)
      landuse_codes <- translation[translation$code==code_conv,'thr_high_code'] 
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
    names(binary_mask) <- mammal$sci_name
    
    # write output
    writeRaster(binary_mask, output_file, overwrite = TRUE)
    message('Processed and saved: ', mammal$sci_name)
  }
}
