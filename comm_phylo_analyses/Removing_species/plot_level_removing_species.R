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
library(forstringr)
library(reshape2)

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
  
#remove unidentified and non-native species 
to_delete <- c("Achinaterum_spp.","Agoseris_spp. ", "Arnica_spp.", 
               "Asteraceae_spp.","Bromopsis_spp.","Bromus_spp.","Cirsium_spp.",
               "Epilobium_spp.","Erigeron_glacialis","Fragaria_spp.",
               "Galium_spp.","Unknown","Poa_spp.","Poa_pratensis","Potentilla_spp.",
                "Senecio_bigelovii","Senecio_crassulus","Solidago_spp.",
               "Veratrum_spp.", "Bromus_inermis",
               "Capsella_bursa-pastoris", "Rumex_cripus","Phleum_pratense",
               "Elymus_glaucus","Leymus_triticoides","Festuca_saximontana",
               "Heracleum_sphondylium","Juncus_drummondii","Melica_spectabilis",
               "Unknown","Taraxacum_officinale","Thlaspi_arvense","Tragopogon_dubius",
               "Bromus_carinatus", "Bromus_ciliatus","Achnatherum_lettermanii")

plots_clean <- plots %>% filter(!Species %in% to_delete)

#check to see all species are in phylogeny
missing_species2 <- plots_clean %>% filter(!Species %in% SB_list$x)

#sum abundance across years
plots_clean_summed <- plots_clean %>% 
  group_by(Site, Species) %>% 
  summarise_if(
    is.numeric,
    sum,
    na.rm = TRUE
  )

#remove summed year column 
plots_clean_summed <- subset(plots_clean_summed, select = -Total)

#fill plots column 
test <- plots_clean_summed %>%
  pivot_longer(cols = c("Plot.1","Plot.2","Plot.3","Plot.4","Plot.5"), 
               names_to = "Plot", values_to = "Total")

#change plot number to not include "plot"
test <- test %>% mutate(across(Plot, gsub, pattern="Plot.", replacement=""))

#rename
plots_list <- test

#create new column to merge site and plot
plots_list$PlotID <- paste(plots_list$Site,"_",plots_list$Plot)

#get rid of spaces 
plots_list <- str_rm_whitespace_df(plots_list)
plots_list$PlotID <- gsub('\\s+', '', plots_list$PlotID)

#make matrix to prune phylogeny----
matrix_for_pruning <- subset(plots_list, select = -c(Total, Site, Year, Plot, PlotID))

#have only one replicate of each species value
matrix_for_pruning <- matrix_for_pruning[!duplicated(matrix_for_pruning), ]
sp_vec <- as.vector(matrix_for_pruning)
matrix_for_pruning$Total <- 1

#check with whole site species list------
edi_df <- read.csv("results/Raw_results_EDI.csv")
edi_df$Species <- gsub('\\s+', '_', edi_df$Species)

check <- edi_df %>% filter(!Species %in% matrix_for_pruning$Species)

#transpose to be matrix and format 
matrix <- dcast(matrix_for_pruning, Total ~ Species)
matrix <- matrix %>% remove_rownames %>% column_to_rownames(var='Total')
row.names(matrix)[row.names(matrix) == "1"] <- "Total"
names(matrix)[names(matrix) == 'Arctostaphylos_uva-ursi'] <- 'Arctostaphylos_uvaursi'

#prune phylogeny-----
#bring in S&B tree
SBtree <- read.tree("ALLMB.tre")

#prune tree
pruned.tree <- treedata(SBtree, unlist(matrix[1, matrix[1,]>0]), 
                        warnings = F)$phy

#check tree
plot(pruned.tree)
is.rooted(pruned.tree)
is.binary(pruned.tree)
Ntip(pruned.tree)

#make species lists for each plot -----
Road1 <- plots_list %>% filter(PlotID %in% "Road_1")
Road2 <- plots_list %>% filter(PlotID %in% "Road_2")
Road3 <- plots_list %>% filter(PlotID %in% "Road_3")
Road4 <- plots_list %>% filter(PlotID %in% "Road_4")
Road5 <- plots_list %>% filter(PlotID %in% "Road_5")
Pfeiler1 <- plots_list %>% filter(PlotID %in% "Pfeiler_1")
Pfeiler2 <- plots_list %>% filter(PlotID %in% "Pfeiler_2")
Pfeiler3 <- plots_list %>% filter(PlotID %in% "Pfeiler_3")
Pfeiler4 <- plots_list %>% filter(PlotID %in% "Pfeiler_4")
Pfeiler5 <- plots_list %>% filter(PlotID %in% "Pfeiler_5")
PBM1 <- plots_list %>% filter(PlotID %in% "PBM_1")
PBM2 <- plots_list %>% filter(PlotID %in% "PBM_2")
PBM3 <- plots_list %>% filter(PlotID %in% "PBM_3")
PBM4 <- plots_list %>% filter(PlotID %in% "PBM_4")
PBM5 <- plots_list %>% filter(PlotID %in% "PBM_5")

#we made individual community matrics for remove in excel becuase it is
#easier to order according to abundance and range size 