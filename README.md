# Areas of Habitat for tropical terrestrial mammals

Based on the workflow conducted in Lumbierres et al 2022 to calculate Areas of Habitat for all tropical terrestrial mammals.

`/scripts`: contains the needed scripts to conduct the workflow

`/Spatial_Data`: contains some of the needed spatial data to conduct the workflow, some files must be generated elsewhere, downloaded manually or downloaded through code

`/Habitats`: contains the needed files to translate ESA CCI land uses to IUCN habitats with the method proposed in Lumbierres et al 2021

### Some general considerations about the translation CCI-IUCN

1. IUCN habitats are described in very general terms, so we have a single category for all kind of forests, grasslands, etc.
2. Since this is based on correlations, some of them at first seem incoherent such as croplands being classified as wetlands, but this concrete case might be explained by rice fields and other inundated croplands which are highly abundant in the tropical range.
3. Some habitats are lost due to low correlations with every land use, such as savannahs. Maybe some cases should be reconsidered.
4. In general, I found this translation a little bit confussing, so maybe it is good that some of you check that I did it right.
