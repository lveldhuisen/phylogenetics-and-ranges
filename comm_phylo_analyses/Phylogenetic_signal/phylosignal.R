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
pruned.tree <- treedata(SBtree, unlist(allsitesmatrix[2,allsitesmatrix[2,]>0]),
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

specieslist <- pruned.tree$tip.label
specieslist <- as.data.frame(specieslist)

##Blomberg's K#####
###abundance#####
###make trait dataframe####
abundance_df <- read.csv("comm_phylo_analyses/Phylogenetic_signal/abundance_trait_data.csv")

abundance_df = subset(abundance_df, select = -c(X, X.1, X.2) )

abundance_df <- abundance_df %>%
  group_by(Species) %>%
  summarise(across(c(Mean_abundance), sum))

abundance_df <- abundance_df[ order(match(abundance_df$Species, specieslist$specieslist)), ]

abundance_df <- abundance_df %>% remove_rownames %>% column_to_rownames(var="Species")
abundance_df <- df2vec(abundance_df, colID=1)

###calculate signal#####
phylosignal(abundance_df, pruned.tree, reps = 5000, checkdata = TRUE)

###range size######
##make community data matrix#### 
allsitesmatrix_range <- read.table("comm_phylo_analyses/Phylogenetic_signal/rangesize_comm_matrix.txt", 
                             sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(allsitesmatrix_range[1,allsitesmatrix_range[1,]>0]),
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

specieslist <- pruned.tree$tip.label
specieslist <- as.data.frame(specieslist)

#get range size data
rangesize_df <- read.csv("comm_phylo_analyses/Phylogenetic_signal/range_size_results.csv")

rangesize_df = subset(rangesize_df, select = -c(EOO) )

rangesize_df <- rangesize_df[-c(67:89), ]

rangesize_df <- rangesize_df %>%
  group_by(Species) %>%
  summarise(across(c(AOO), sum))

rangesize_df <- rangesize_df[ order(match(rangesize_df$Species, specieslist$specieslist)), ]

rangesize_df <- rangesize_df %>% remove_rownames %>% column_to_rownames(var="Species")
rangesize_df <- df2vec(rangesize_df, colID=1)

###calculate signal#####
phylosignal(rangesize_df, pruned.tree, reps = 5000, checkdata = TRUE)
