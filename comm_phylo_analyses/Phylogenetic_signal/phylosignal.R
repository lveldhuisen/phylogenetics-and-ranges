#calculate phylogenetic diversity of abundance and range size
library(tidyverse)
library(dplyr)
library(picante)
library(geiger)
library(ape)

#import S&B phylogeny--------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

##make community data matrix#### 
allsitesmatrix <- read.table("comm_phylo_analyses/Phylogenetic_signal/species_matrix_fortree.txt", 
                                   sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_abundance_matrix[32,
                                                            PBM_abundance_matrix[32,]>0]),
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

specieslist <- pruned.tree$tip.label

##Blomberg's K#####
###abundance#####
###make trait dataframe####
abundance_df <- read.csv("comm_phylo_analyses/Phylogenetic_signal/abundance_trait_data.csv")

abundance_df = subset(abundance_df, select = -c(X, X.1, X.2) )

abundance_df <- abundance_df %>%
  group_by(Species) %>%
  summarise(across(c(Mean_abundance), sum))

abundance_df <- abundance_df[ order(match(abundance_df$Species, specieslist$specieslist)), ]

###calculate signal#####
phylosignal(abundance_df, pruned.tree, reps = 5000, checkdata = TRUE)

###range size######