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
PBM_abundance_matrix <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/PBM_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_abundance_matrix[32,PBM_abundance_matrix[32,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

##Blomberg's K#####
###abundance#####
phylosignal(x, pruned.tree, reps = 5000, checkdata = TRUE)

###range size######