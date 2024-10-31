#calculate range size for all species in the community 

library(rgbif)
library(dplyr)
library(red)

#download data from GBIF--------------------------

#get taxon key 
name_backbone("Rumex_densiflorus")

#download occurrence data from GBIF directly
occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"),
  pred_gte("year",1990),
  pred("taxonKey",2888957),
  format = "SIMPLE_CSV",
  user="leah.veldhuisen", 
  pwd="Columbia2305", 
  email="leah.veldhuisen@gmail.com"
)

d <- occ_download_get('0009396-241024112534372') %>%
  occ_download_import()

##clean up data list to only have species name, lat and long##############
df1 <- d %>%
  select(species, decimalLatitude, decimalLongitude,coordinateUncertaintyInMeters)
df1 <- df1[df1$coordinateUncertaintyInMeters < 1000, ]
df1 <- na.omit(df1)
df1 <- df1 %>%
  select(species, decimalLatitude, decimalLongitude)
df1$species <- sub(" ",".", df1$species)

#extract species occurrences
species_occ <- df1

#get occurrence points only
points_species <- data.frame(species_occ[,c("decimalLongitude", "decimalLatitude")])

#calculate AOO and EOO with red package----------------------------
#AOO
aoo(points_species)

#EOO
eoo(points_species)
