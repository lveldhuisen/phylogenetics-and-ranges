library(tidyverse)
library(dplyr)
library(picante)
library(geiger)
library(ape)
library(vegan)
library(forcats)
library(broom)
library(janitor)
library(patchwork)
library(car)

#data prep-----
plots <- read.csv("comm_phylo_analyses/RMBL_abundance_2021-2022.csv")

#match species names with phylogeny----
#bring in S&B taxa list
SB_list <- read.csv("SmithBrown18_taxaalist.csv")

#show missing species
missing_species <- plots %>% filter(!Species %in% SB_list$x)

#replace species names with phylogeny replacements, updated names, remove typos, etc
plots$Species[plots$Species == 'Agoseris_glauca'] <- 'Agoseris_glauca_var._dasycephala'
plots$Species[plots$Species == 'Arctostaphylos_uva-ursi'] <- 'Arctostaphylos_uvaursi'
plots$Species[plots$Species == 'Bromelica_spectabilis'] <- 'Melica_spectabilis'
plots$Species[plots$Species == 'Carex_spp.'] <- 'Carex_albonigra'
plots$Species[plots$Species == 'Chamerion _angustifolium'] <- 'Chamerion_angustifolium'
plots$Species[plots$Species == 'Chrysothamnus viscidiflorus'] <-'Chrysothamnus_viscidiflorus'
plots$Species[plots$Species == 'Gayophytum_diffusum'] <- 'Gayophytum_diffusum_subsp._diffusum'
plots$Species[plots$Species == 'Gayophytum_ramosissimum '] <-'Gayophytum_ramosissimum'
plots$Species[plots$Species == 'Erigeron_elatior'] <- 'Erigeron_grandiflorus'
plots$Species[plots$Species == 'Helianthella_quinquenervis'] <- 'Helianthella_uniflora'
plots$Species[plots$Species == 'Hydrophyllum_capitatum'] <- 'Hydrophyllum_capitatum_var._capitatum'
plots$Species[plots$Species == 'Penstemon_crandallii'] <- 'Penstemon_crandallii_subsp._crandallii'
plots$Species[plots$Species == 'Potentilla_gracilis'] <- 'Potentilla_gracilis_var._flabelliformis'
plots$Species[plots$Species == 'Rhodiola_intregrifolia'] <- 'Rhodiola_integrifolia'
plots$Species[plots$Species == 'Veratrum_tenuipetalum'] <- 'Veratrum_virginicum'
plots$Species[plots$Species == 'Elymus_triticoides'] <- 'Leymus_triticoides'
plots$Species[plots$Species == 'Fleum_pretense'] <- 'Phleum_pratense'
plots$Species[plots$Species == 'Erigeron_eatoni'] <- 'Erigeron_eatonii'
  
#remove unidentified species 
to_delete <- c("Achinaterum_spp.","Agoseris_spp. ", "Arnica_spp.", 
               "Asteraceae_spp.","Bromopsis_spp.","Bromus_spp.","Cirsium_spp.",
               "Epilobium_spp.","Erigeron_glacialis","Fragaria_spp.",
               "Galium_spp.","Unknown","Poa_spp.","Poa_pratensis","Potentilla_spp.",
                "Senecio_bigelovii","Senecio_crassulus","Solidago_spp.",
               "Veratrum_spp.", "Bromus_inermis",
               "Capsella_bursa-pastoris", "Rumex_cripus","")

plots_clean <- plots %>% filter(!Species %in% to_delete)

#check to see all species are in phylogeny
missing_species <- plots_clean %>% filter(!Species %in% SB_list$x)
