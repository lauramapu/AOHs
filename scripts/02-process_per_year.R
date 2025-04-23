
### 1. process ESA CCI and translate
# crop/mask with tropical forest
# from this step the process is per year (ESA-CCI)
# trial done with with 2010

rm(list=ls())

library(terra)
library(dplyr)
library(data.table)
library(readxl)

forest <- vect('Spatial_Data/Tropical_Forest/tropicalmask.shp')

# load translation csv 
habitats <- read.csv('Habitats/translation_by_lumbierres.csv')
is <- habitats$Value
become <- habitats$code_lumbierres

# load and resample SRTM to match ESA
srtm <- rast('Spatial_Data/SRTM90/strm_300m_trop.tif') %>%
  resample(x, method='bilinear')
# now we round srtm value to /10
srtm <- app(srtm, fun = function(x) {round(x/10)})

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
lapply(esa_files, function(esa) {
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
  
  # reclass land uses to habitats 
  y <- subst(x, is, become)
  # y <- classify(ESA, cbind(is, become)) # another solution
  
  out_name <- paste0(output_dir, 'esa_', year, '_reclass.tif')
  writeRaster(y, out_name, overwrite=T)
  
  # multiply per 1000
  y <- app(y, fun = function(i) {i*1000})
  # and sum to srtm
  habitat_elev <- srtm + y

  # the resulting tif has pixels coded so that first 3 numbers are habitat (from 1 to 150)
  # and last 3 numbers are elev/10 (should range from 0 to 700 or something like that because top elevation is ~7000m)
  
  out_name <- paste0(output_final, 'baselayer_', year, '.tif')
  writeRaster(habitat_elev, out_name, overwrite=T)
  
  cat(out_name, ' done/n')
  
})

rm(srtm, esa_files, x, y, habitat_elev)

### 2. finally, generate AoH for each species using the habitat-elevation tif

# this process can't be parallelized in local because not enough RAM (some iterations take more than 50gb)
# if you want to send this to cesga and parallelize, you must wrap/unwrap spatraster/spatvector objects
# outside and inside the loop because they are non-exportable
# (see https://future.futureverse.org/articles/future-4-non-exportable-objects.html for more info)

hab_pref <- readRDS('Habitats/mammal_habitat_preferences.rds')
elev_range <- readRDS('Habitats/mammal_elevation_ranges.rds')

mammals <- vect('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp')

base <- rast('Spatial_Data/AOHs/aoh_2010.tif') %>%
  project(mammals, method='near')

for (i in seq_along(mammals)) {
  
  # aoh <- unwrap(aoh_w)
  # mammals <- unwrap(mammals_w)
  
  mammal <- mammals[i,]
  
  r <- aoh %>%
    crop(mammal) %>%
    mask(mammal)
  
  habitat_codes <- as.numeric(sub('_.*', '',hab_pref[[mammal$sci_name]]$Habitat_Code))
  elevation_range <- elev_range[[mammal$sci_name]]
  
  # 1) compute the elevation thresholds (in the raster's units)
  lo_e <- elevation_range$Lower_Elevation_Limit / 10
  hi_e <- elevation_range$Upper_Elevation_Limit / 10
  
  # 2) initialize a logical raster (all FALSE) matching 'r'
  cond <- r
  cond[] <- FALSE
  
  # 3) for each habitat code, set TRUE where r ∈ [code*1000 + lo_e, code*1000 + hi_e]
  for(code in habitat_codes){
    
    #if (lo_e == NA | hi_e == NA) {lo_e<-0; hi_e<-0} # set elev to zero if there is no data
    
    minv <- code * 1000 + lo_e
    maxv <- code * 1000 + hi_e
    
    cond <- cond |
      ((r >= minv) & (r <= maxv))
  }
  
  # 5) convert logical → binary (1/0)
  binary_mask <- ifel(cond, 1, NA)
  names(binary_mask) <- mammal$sci_name
  
  # 6) write out if you like
  writeRaster(binary_mask, paste0('Spatial_Data/AOHs/',mammal$sci_name,'.tif'), overwrite=TRUE)
}

registerDoSEQ()  
