
############
# full workflow to filter IUCN distribution maps per elevation range and habitat correspondence for each species
############################

# there's a difference with the original work by Lumbierres et al (2022) because we need to use ESA-CCI
# to have all years in our study range while she uses CGLS-LC100
# also habitat translation could be better for some of them (croplands classified as wetlands :/)

# if an error occurs when crop/mask (terra) telling something like too many values to write read this:
# https://github.com/rspatial/terra/issues/1686
# and this to do the clean installation:
# https://stackoverflow.com/questions/70925962/r-packages-raster-fail-to-upload-while-searching-for-terra-last-version


### 1. crop/mask mammal ranges with tropical forest mask and extract wanted distributions
# great part of this code extracted from 02_IUCN_Ranges by Iago)

rm(list=ls())

library(terra)
library(dplyr)
library(data.table)
library(readxl)

#Import pantropical forest zone shp
forest <- terra::vect("Spatial_Data/Tropical_Forest/TropicalForest.shp")
# dissolve to fix bad geometries and save
forest$id <- 1
forest <- aggregate(forest, by = "id", dissolve = TRUE)
writeVector(forest, 'Spatial_Data/Tropical_Forest/tropicalmask.shp', overwrite=T)

#Import mammal distribution ranges from IUCN.
mammal_ranges <- terra::vect("Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/MAMMALS_TERRESTRIAL_ONLY.shp")

# here I select the ones used in Lumbierres et al 2022:
# we selected range polygons with extant and probably extant presence; native, reintroduced,
# and assisted colonization origin; and resident seasonality for non-migratory species
# (all mammals and the 8,979 non-migratory birds)

legend <- unique(mammal_ranges$legend)

extant <- legend[grep('^Extant ', legend)]
extant <- extant[c(1,4,5)]
pextant <- legend[grep('^Probably Extant ', legend)]
pextant <- pextant[1]
wanted_ranges <- c(extant, pextant)
# revise this because I dont have it clear

# Filter out wanted polygons
mammal_ranges <- mammal_ranges[mammal_ranges$legend %in% wanted_ranges, ]

# Combine multiple polygons of each species into one.
mammal_ranges <- mammal_ranges[, c("id_no",'sci_name')]
mammal_ranges <- aggregate(mammal_ranges, by = "id_no", dissolve = TRUE)

#Intersect mammal distribution range with tropical forest extent
forest <- project(forest, mammal_ranges)
tropical_ranges <- terra::intersect(mammal_ranges, forest)
rm(mammal_ranges)

# save single shp
terra::writeVector(tropical_ranges, 
                   paste0('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp'),
                   overwrite = TRUE)

#Save each species as single shp
for (i in seq_along(unique(tropical_ranges$sci_name))) {
  
  sci_name <- unique(tropical_ranges$sci_name)[i]
  
  terra::writeVector(tropical_ranges[which(tropical_ranges$sci_name == sci_name), ], 
                     paste0('Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/', sci_name, '.shp'),
                     overwrite = TRUE)
}
  
rm(tropical_ranges)

### 2. extract habitat and elevation ranges information on each species with the IUCN API
# (great part extracted from script 04_Habitat_Preferences done by Iago, just for mammals)

library(rredlist)
library(iucnredlist) # devtools::install_github("IUCN-UK/iucnredlist")
library(stringr)
library(tictoc)

# initialize lists to store habitats and elevation ranges
mammal_habitat_preferences <- list()
mammal_elevation_ranges <- list()

# initialize token and api
token <- "D98jKWcys1KfrHVAngXsno85KWjNw6s2qyyt"  
api <- init_api(token) 

# import tropical terrestrial mammals created before
mammals <- vect("Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp")
# our field of interest is 'id_no' which is an identifier for each species
species <- data.frame(id = unique(mammals$id_no), sci_name = unique(mammals$sci_name))

# Temporizador de ejecución
tic()  # Iniciar el temporizador

# Iteración sobre cada fila del dataframe (cada especie)
for (i in 1:nrow(species)) {
  print(i)  # Imprimir el número de iteración (especie actual)
  
  # Obtain ID and 
  species_id <- species[i, "id"]
  species_name <- as.character(species[i, "sci_name"])
  
  # Verify if species was already processes
  if (species_name %in% names(mammal_habitat_preferences)) {next}
  
  # Extract most recent assessment 
  assessment_raw <- iucnredlist::assessments_by_sis_id(api, species_id)
  assessment_raw <- assessment_raw %>% 
    dplyr::filter(year_published == max(year_published)) %>%
    dplyr::filter(latest == TRUE)
  
  if (is.null(assessment_raw$assessment_id)) {
    print(paste("No data for species:", species_name, "Skipping..."))
    next  
  }
  
  # Obtain assessments from species
  a_data <- assessment_data_many(api, assessment_raw$assessment_id, wait_time = 0.5)
  habitats <- extract_element(a_data, "habitats")
  
  #filter suitable habitat of the species
  if (!is.null(habitats) && nrow(habitats) > 0) {
    habitats <- habitats %>% 
      filter(suitability == "Suitable") %>% 
      select("Description" = description, "Habitat_Code" = code) %>%  
      distinct() 
    
    # save results in list
    mammal_habitat_preferences[[species_name]] <- habitats
  }
  
  # extract elevation info if it exists
  elevation <- extract_element(a_data, "supplementary_info")
  
  # just verify if there is info about elevation limits
  if (!is.null(elevation) && all(c("upper_elevation_limit", "lower_elevation_limit") %in% colnames(elevation))) {
    elevation_data <- elevation %>% 
      select("Upper_Elevation_Limit" = upper_elevation_limit, 
             "Lower_Elevation_Limit" = lower_elevation_limit) %>% 
      mutate("IUCN_Species" = species_name) #add species name
  } else {
    # If there is no data, create the df with NA values
    elevation_data <- data.frame(
      "Upper_Elevation_Limit" = NA,
      "Lower_Elevation_Limit" = NA,
      "IUCN_Species" = species_name
    )
  }
  
  mammal_elevation_ranges[[species_name]] <- elevation_data
}

toc()  

# Save results 
saveRDS(mammal_habitat_preferences, "Habitats/mammal_habitat_preferences.rds")
saveRDS(mammal_elevation_ranges, "Habitats/mammal_elevation_ranges.rds")

### 3. process ESA CCI + SRTM and translate
# crop/mask with tropical forest
# from this step the process is per year (ESA-CCI)
# trial with 2010

forest <- vect('Spatial_Data/Tropical_Forest/tropicalmask.shp')
ESA <- rast("Spatial_Data/ESA-LC/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2010-v2.0.7.tif")

forest <- project(forest, ESA)

ESA <- ESA %>% 
  crop(forest) %>% 
  mask(forest)

writeRaster(ESA, 'Spatial_Data/ESA-LC/esa_2010_trop.tif', overwrite=T)

# do the same with SRTM but also merging because we have 3 rasters for the whole dataset
# SRTM was extracted and resampled using GEE
srtm_files <- list.files('Spatial_Data/SRTM90', pattern='.tif$', full.names=TRUE)
srtm_list <- lapply(srtm_files, rast)
srtm_merged <- do.call(merge, srtm_list)

forest <- project(forest, srtm_merged)

srtm <- srtm_merged %>%
  crop(forest) %>%
  mask(forest)
writeRaster(srtm, 'Spatial_Data/SRTM90/strm_300m_trop.tif', overwrite=T)

# translate ESA-CCI land uses to habitats as described in IUCN basing on Lumbierres et al 2021

# first load correspondences between esa and iucn and esa legend
translation <- read_xlsx('Habitats/esa-habitats.xlsx') %>%
  mutate(across(3:13, ~ {
    x <- ifelse(.x == "-", NA, .x) # replace '-' with NA
    as.numeric(x) # convert to numeric
  }))
cci_legend <- read_xlsx('Habitats/CCI-LC_Maps_Legend.xlsx') %>%
  # new column to match land use groups used in Lumbierres et al 2021
  mutate(lumbierres = c(
    NA,
    "Cropland, rainfed",
    "Cropland, rainfed: herbaceous cover",
    "Cropland, rainfed: tree or shrub cover",
    "Cropland irrigated or post-flooding",
    "Mosaic cropland (>50%) /natural vegetation (tree, shrub, herbaceous cover)(<50%)",
    "Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%)/cropland (<50%)",
    "Tree cover, broadleaved, evergreen, closed to open (>15%)",
    "Tree cover, broadleaved, deciduous, closed to open (>15%)",
    "Tree cover, broadleaved, deciduous, closed to open (>15%)",
    "Tree cover, broadleaved, deciduous, closed to open (>15%)",
    "Tree cover, needleleaved, evergreen, closed to open (>15%)",
    "Tree cover, needleleaved, evergreen, closed to open (>15%)",
    "Tree cover, needleleaved, evergreen, closed to open (>15%)",
    "Tree cover, needleleaved, deciduous, closed to open (>15%)",
    "Tree cover, needleleaved, deciduous, closed to open (>15%)",
    "Tree cover, needleleaved, deciduous, closed to open (>15%)",
    "Tree cover, mixed leaf type (broadleaved and needleleaved)",
    "Mosaic tree and shrub (>50%) / herbaceous cover (<50%)",
    "Mosaic herbaceous cover (>50%) / tree and shrub (<50%)",
    "Shrubland",
    "Shrubland",
    "Shrubland",
    "Grassland",
    "Lichens and mosses",
    "Sparse vegetation (tree, shrub,herbaceous cover)(<15%)",
    "Sparse vegetation (tree, shrub,herbaceous cover)(<15%)",
    "Sparse vegetation (tree, shrub,herbaceous cover)(<15%)",
    "Sparse vegetation (tree, shrub,herbaceous cover)(<15%)",
    "Tree cover, flooded, fresh or brakish water",
    "Tree cover, flooded, saline water",
    "Shrub or herbaceous cover, flooded, fresh/ saline/brakish water",
    "Urban area",
    "Bare areas",
    "Bare areas",
    "Bare areas",
    "Water bodies",
    "Water bodies"
  ))

# match habitat class with highest positive value to each land use
numeric_cols <- translation %>% select(where(is.numeric)) %>% names()
translation <- translation %>%
  rowwise() %>%
  mutate(
    max_value = max(c_across(all_of(numeric_cols)), na.rm = TRUE),
    habitat = numeric_cols[which.max(c_across(all_of(numeric_cols)))]
  ) %>%
  ungroup()

# iucn habitats with code and description (extracted from sup materials in Lumbierres et al 2021)
iucn <- tibble::tibble(
  code = c(
    "H1", "H2", "H3", "H4", "H5", "H6", "H8", 
    "H14.1", "H14.2", "H14.3", "H14.6", "H14.4", "H14.5", "H15"
  ),
  name = c(
    "Forest",
    "Savanna",
    "Shrubland",
    "Grassland",
    "Wetlands",
    "Rocky Areas",
    "Desert",
    "Artificial arable and pasture lands: Arable Land",
    "Artificial arable and pasture lands: Pastureland",
    "Artificial degraded forest and plantation: Plantations",
    "Artificial degraded forest and plantation: Degraded Forest",
    "Artificial urban and rural gardens: Rural Gardens",
    "Artificial urban and rural gardens: Urban Areas",
    "Artificial Aquatic"
  ),
  description = c(
    "Forest consists of a continuous stand of trees and includes both forested areas (generally with a closed canopy) and wooded areas.",
    "Savannas are transitional between grasslands and forests. They are ecosystems dominated by a grass ground cover with an overstorey of widely spaced trees and shrubs.",
    "Also referred to as scrub, bushland and thicket.",
    "Native grasslands are comprised of grasses and broadleaved herbaceous plants, and are either without woody plants, or the latter are very sparsely distributed.",
    "Areas of marsh, fen, peatland or water, whether natural or artificial, permanent or temporary, with water that is static or flowing, fresh, brackish or salt, including areas of marine water the depth of which at low tide does not exceed six meters.",
    "Cliffs, mountain peaks, talus, feldmark.",
    "Consists of arid landscapes with a sparse plant cover, except in depressions where water accumulates. The sandy, stony, or rocky substrate contributes more to the appearance of the landscape than does the vegetation.",
    "Includes cereal fields, rice paddies, perennial crops, orchards and groves.",
    "Includes fertilized or re-seeded permanent grasslands, sometimes treated with selective herbicides, with very impoverished flora and fauna. Also includes secondary grasslands and wooded farmland.",
    "A plantation is an intentional planting of a crop, on a larger scale, usually for uses other than cereal production or pasture. The term is currently most often used for plantings of trees and shrubs. The term tends also to be used for plantings maintained on economic bases other than that of subsistence farming.",
    "Former subtropical or tropical forest that has been extensively cleared or impacted by human activities. Often there is some degree of regeneration or there are small fragments of forest remaining.",
    "Rural gardens are located in a rural setting, serving families whose main income comes from wage labor (rural or urban)... [description continues, see source]",
    "Usually metropolitan and commercial areas dominated by asphalt, concrete and roof. Includes buildings, lawns and parks.",
    "These are human-made wetland habitats."
  )
)
iucn <- iucn %>%
  dplyr::mutate(code_numeric = c(
    1, 2, 3, 4, 5, 6, 8, 141, 142, 143, 146, 144, 145, 150
  )) %>%
  dplyr::mutate(habitat = c(
    "Forest",
    "Savanna",
    "Shrubland",
    "Grassland",
    "Wetlands",
    "Rocky Areas",
    "Desert",
    "Artificial arable and pasture lands",
    "Artificial arable and pasture lands",
    "Artificial degraded forest and plantation",
    "Artificial degraded forest and plantation",
    "Artificial urban areas and rural gardens",
    "Artificial urban areas and rural gardens",
    "Artificial Aquatic"
  ))
code_numeric <- iucn$code_numeric
iucn$code_lumbierres <- c(1, 2, 3, 4, 5, 6, 8, 141, 141, 143, 143, 144, 144, 150)

# new df matching things
habitats <- cci_legend %>%
  select(Value, Label, lumbierres)

colnames(translation)[1] <- 'lumbierres'
habitats <- left_join(habitats, translation, by = 'lumbierres')

habitats <- left_join(habitats, iucn, by = 'habitat')
habitats <- habitats %>%
  select(Value, Label, lumbierres, habitat, code_lumbierres)
habitats$code_lumbierres[1] <- 0

# load cropped esa tif
ESA <- rast('Spatial_Data/ESA-LC/esa_2010_trop.tif')
# reclass land uses to habitats in esa tif
is <- habitats$Value
become <- habitats$code_lumbierres
x <- subst(ESA, is, become)
writeRaster(x, 'Spatial_Data/ESA-LC/esa_2010_trop_reclass.tif', overwrite=T)
# y <- classify(ESA, cbind(is, become)) # another solution

# now join ESA reclassed and SRTM so that we have it all in a single tif
# first we need to resample SRTM to match ESA
srtm <- rast('Spatial_Data/SRTM90/strm_300m_trop.tif') %>%
  resample(x, method='bilinear')
# now we round srtm value to /10
srtm <- app(srtm, fun = function(x) {round(x/10)})
# also multiply ESA per 1000
x <- app(x, fun = function(x) {x*1000})
# and sum both
habitat_elev <- srtm + x
# the resulting tif has pixels coded so that first 3 numbers are habitat (from 1 to 150)
# and last 3 numbers are elev/10 (should range from 0 to 700 or something like that because top elevation is ~7000m)

dir.create('Spatial_Data/AOHs')
writeRaster(habitat_elev, 'Spatial_Data/AOHs/aoh_2010.tif', overwrite=T)

# 4. finally, generate a tif file or AOH for each species using the habitat-elevation tif
# erase everything and load all needed files
rm(list=ls())

hab_pref <- readRDS('Habitats/mammal_habitat_preferences.rds')
elev_range <- readRDS('Habitats/mammal_elevation_ranges.rds')

mammals <- vect("Spatial_Data/IUCN_Range_Maps_Terrestrial_Mammals/Terrestrial_Mammals_TropicalRanges.shp")

aoh <- rast('Spatial_Data/AOHs/aoh_2010.tif') %>%
              project(mammals, method='near')
writeRaster(aoh, 'Spatial_Data/AOHs/aoh_2010_reproj.tif')

aoh <- rast('Spatial_Data/AOHs/aoh_2010_reproj.tif')

for (i in seq_along(mammals)) {
  
          # aoh <- unwrap(aoh_w)
          # mammals <- unwrap(mammals_w)
          
          mammal <- mammals[i,]
          
          r <- aoh %>%
            crop(mammal) %>%
            mask(mammal)
          
          habitat_codes <- as.numeric(sub("_.*", "",hab_pref[[mammal$sci_name]]$Habitat_Code))
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
          writeRaster(binary_mask, paste0('Spatial_Data/AOHs/',mammal$sci_name,".tif"), overwrite=TRUE)
}

registerDoSEQ()  
