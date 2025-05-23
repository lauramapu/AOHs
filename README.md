# Areas of Habitat for tropical terrestrial mammals and birds

Based on the workflow conducted in Lumbierres et al 2022 to calculate Areas of Habitat for all tropical terrestrial mammals and birds.

`/scripts`: contains the needed scripts to conduct the workflow

`/Spatial_Data`: contains some of the needed spatial data to conduct the workflow, some files must be generated elsewhere, downloaded manually or downloaded through code

`/Habitats`: contains the needed files to translate ESA CCI land uses to IUCN habitats with the method proposed in Lumbierres et al 2021

### Some general considerations about the translation CCI-IUCN

1. IUCN habitats are described in very general terms, so we have a single category for all kind of forests, grasslands, etc.
2. Some habitats get no land use correlations in certain tertiles, such as rocky areas and artificial aquatic in low tertile, and artificial degraded forest and platation in high tertile, so might be lost in translation (look for them in tokyo :P). Same happens with some land uses that has no high correlations with any habitat, such as mosaic classes, so pixels classified with those land uses are not suitable for any species due to this translation.
3. I did this only with the highest tertile of the positive correlations but it can be easily modified (or you can just repeat the AOH loop replacing 'thr_high_code' with 'thr_mid/low_code').

### Other general considerations

1. When filtering habitats per species with the IUCN package, some species get no habitat codes and in those cases I assumed no habitats are suitable so no AOH can be generated for the species. This happens for some birds (I think) and for some mammals that has some habitats coded as Unkwonn instead of Suitable or Marginal. There is a silenced line in the AOH generation in which you can set all habitats as suitable or you can also add Unknown habitats when extracting data through the IUCN API.
2. I suppose due to data availability issues, some species have their elevation range set to the same altitude for both min and max (see Dendromus vernayi). This limits a lot the generated AOH.
3. I am almost sure that there are no issues at all with the mammals workflow because I deeply checked results, but I cannot say the same for birds because I did not had enough time.
4. Birds with category Critical Endangered and with code 4 (possibly extinct) in attribute 'presence' are only 2 while in the original paper Lumbierres says they are 22, I guess either the category or code for those species changed but I did not had time to check that.
5. I compared my results for some species with the results obtained by Lumbierres and some notable differences are found. However, results are not fully comparable since we used different land uses and thus different translation, and also the range polygons used might be also different since there may have been updates in some species.
