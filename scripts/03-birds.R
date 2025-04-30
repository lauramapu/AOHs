
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

forest <- vect('Spatial_Data/Tropical_Forest/tropicalmask.shp')

### 1. crop BirdLife distributions with tropical forest mask and save
# this step extracted from script 02 by iago

# skip if file already exists
out_name <- "Spatial_Data/BirdLife_Range_Maps_Birds/Birds_TropicalRanges.shp"

if (file.exists(out_name)) {
  cat("Skipping BirdLife tropical intersect - output already exists:", out_name, "\n")
}

else {

  #Import bird distribution ranges from IUCN/Bird_Life
  bird_ranges <- terra::vect("Spatial_Data/BirdLife_Range_Maps_Birds/BOTW.gdb")
  #head(bird_ranges)
  
  # here I select the ones used in Lumbierres et al 2022:
  
  # "we selected range polygons with extant and probably extant presence; native, reintroduced,
  # and assisted colonization origin; and resident seasonality for non-migratory species
  # (all mammals and the 8,979 non-migratory birds). For migratory birds (1,816 species), we kept
  # separated the ranges for breeding (1,446 species), non-breeding (1,550 species) and a
  # combination of resident and uncertain (1,290 species) seasonality. [...]
  # For 18 mammal and 22 bird species categorized as Critical Endangered, there were no presence
  # polygons coded as extant or probably extant. To assist the conservation of these species,
  # we produced AOH maps using the possibly extinct polygon for these taxa."
  
  # codes of attributes are available in Spatial_Data/IUCN_Standard_attributes_for_spatial_data_v1.20_2024.xlsx
  # our interest codes are the following:
  
  # $presence -> extant == 1; Probably extant == 2; Possibly Extinct == 4
  # $origin -> native == 1; reintroduced == 2, Assisted colonisation == 6
  # $seasonal -> resident == 1; 2 == breeding; 3 == non-breeding; passant == 4; uncertain == 5
  
  # first for the non-migratory
  
  # spatvector to df to operate better
  bird_df <- as.data.frame(bird_ranges, geom = "WKT")
  
  # extract CR species
  critical_species <- unique(bird_df$sci_name[bird_df$category == "CR"])
  
  # find CR species with no polygons $presence == 1 or 2
  cr_no_presence_12 <- bird_df %>%
    filter(sci_name %in% critical_species) %>%
    group_by(sci_name) %>%
    summarise(has_presence_12 = any(presence %in% c(1, 2))) %>%
    filter(!has_presence_12) %>%
    pull(sci_name)
  
  # set conditions
  # first for the original criteria
  condition_original <- 
    (bird_ranges$presence %in% c(1, 2)) &
    (bird_ranges$origin %in% c(1, 2, 6)) &
    (bird_ranges$seasonal == 1)
  
  # add presence == 4 polygons for CR species lacking of presence == 1 or 2
  condition_cr_presence4 <- 
    (bird_ranges$sci_name %in% cr_no_presence_12) & 
    (bird_ranges$presence == 4)
  
  # combine conditions
  final_condition <- condition_original | condition_cr_presence4
  
  # subset
  non_migratory <- bird_ranges[final_condition, ]
  
  # migratory
  
  # set conditions
  # first for the original criteria
  condition_breeding <- 
    (bird_ranges$presence %in% c(1, 2)) &
    (bird_ranges$origin %in% c(1, 2, 6)) &
    (bird_ranges$seasonal == 2)
  condition_nonbreeding <- 
    (bird_ranges$presence %in% c(1, 2)) &
    (bird_ranges$origin %in% c(1, 2, 6)) &
    (bird_ranges$seasonal == 3)
  condition_resident <- 
    (bird_ranges$presence %in% c(1, 2)) &
    (bird_ranges$origin %in% c(1, 2, 6)) &
    (bird_ranges$seasonal %in% c(1, 5))
  
  # add presence == 4 polygons for CR species lacking of presence == 1 or 2
  condition_cr_presence4 <- 
    (bird_ranges$sci_name %in% cr_no_presence_12) & 
    (bird_ranges$presence == 4)
  
  # combine conditions
  breeding <- condition_breeding | condition_cr_presence4
  nonbreeding <- condition_nonbreeding | condition_cr_presence4
  resident_uncertain <- condition_resident | condition_cr_presence4
  
  # Combine multiple polygons of each species into one.
  bird_ranges <- aggregate(bird_ranges, by = "sisid", dissolve = TRUE)
  bird_ranges <- bird_ranges[, c("sisid","sci_name")]
  names(bird_ranges) <- c("IUCN_ID","IUCN_Species")
  
  #Project mammal distribution range to tropical forest shapefile
  bird_ranges <- terra::project(bird_ranges, crs(forest))
  
  tropical_ranges <- terra::intersect(bird_ranges, forest)
  tropical_ranges <-terra::simplifyGeom(tropical_ranges, tolerance = 0.01) #keep tolerance super low
  #not simplified ranges can not be saved due to size of file
  rm(bird_ranges)
  
  terra::writeVector(tropical_ranges, out_name, overwrite = TRUE)

}

### 2. extract IUCN information on habitat and elevation range for each species


