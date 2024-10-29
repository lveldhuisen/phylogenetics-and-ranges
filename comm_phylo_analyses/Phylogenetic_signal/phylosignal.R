#calculate phylogenetic signal of abundance and range size, also 
#contains code for figure 2 A-B

library(tidyverse)
library(dplyr)
library(picante)
library(geiger)
library(ape)
library(phytools)
library(viridis)
library(viridisLite)
library(patchwork)

#import S&B phylogeny--------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

##make community data matrix#### 
allsitesmatrix <- read.table("comm_phylo_analyses/Phylogenetic_signal/species_matrix_fortree.txt", sep = "\t", header = T, row.names = 1)

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
abundance_df <- read.csv("comm_phylo_analyses/Phylogenetic_signal/abundance_trait_data_new.csv")

colnames(abundance_df)[colnames(abundance_df) == 'X'] <- 'Species'
colnames(abundance_df)[colnames(abundance_df) == 'x'] <- 'Mean_abundance'

abundance_df <- abundance_df %>%
  group_by(Species) %>%
  summarise(across(c(Mean_abundance), sum))

abundance_df <- abundance_df[ order(match(abundance_df$Species, 
                                          specieslist$specieslist)), ]

abundance_df <- abundance_df %>% remove_rownames %>% column_to_rownames(var="Species")
abundance_df <- df2vec(abundance_df, colID=1)

write.csv(abundance_df, file = "comm_phylo_analyses/Phylogenetic_signal/abundance_trait_data_new.csv")

###calculate signal#####
phylosignal(abundance_df, pruned.tree, reps = 5000, checkdata = TRUE) #with picante
phylosig(pruned.tree, abundance_df, method="K", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=10) #with phytools, gives same answer

###range size######
##make community data matrix#### 
allsitesmatrix_range <- read.table("comm_phylo_analyses/Phylogenetic_signal/rangesize_matrix.txt", sep = "\t", header = T, row.names = 1)

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
colnames(rangesize_df)[colnames(rangesize_df) == 'AOO..km2.'] <- 'AOO'
rangesize_df <- rangesize_df %>%
  group_by(Species) %>%
  summarise(across(c(AOO), sum))

rangesize_df <- rangesize_df[ order(match(rangesize_df$Species, specieslist$specieslist)), ]

rangesize_df <- rangesize_df %>% remove_rownames %>% column_to_rownames(var="Species")
rangesize_df <- df2vec(rangesize_df, colID=1)

###calculate signal#####
phylosignal(rangesize_df, pruned.tree, reps = 5000, checkdata = TRUE)

##pagels lambda####
###abundance######
phylosig(pruned.tree, abundance_df, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=10)

fitContinuous(pruned.tree, abundance_df, SE = 0, model = c("lambda"), bounds= list(), 
              control = list(method = c("subplex","L-BFGS-B"), niter = 5000, 
                             FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)

###range size#####
phylosig(pruned.tree.range, rangesize_df, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=10)

fitContinuous(pruned.tree.range, rangesize_df, SE = 0, model = c("lambda"), 
              bounds= list(), control = list(method = c("subplex","L-BFGS-B"), 
                                             niter = 500000, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)

#make phylogeny with traits mapped on----------
#abundance-----
##make community data matrix#### 
allsitesmatrix <- read.table("comm_phylo_analyses/Phylogenetic_signal/comm_matrix_forfig.txt", 
                             sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree.forfig <- treedata(SBtree, unlist(allsitesmatrix[1,allsitesmatrix[1,]>0]),
                        warnings = F)$phy
plot(pruned.tree.forfig)
is.rooted(pruned.tree)

specieslist <- pruned.tree$tip.label
specieslist <- as.data.frame(specieslist)

#trait data
abundance <- read.csv("comm_phylo_analyses/Phylogenetic_signal/abundance_trait_data_new.csv")
abundance <- abundance %>% 
  rename(species = X,
         trait_value = x) 
abundance <- abundance %>% remove_rownames %>% column_to_rownames(var="species")
abundance <- df2vec(abundance, colID=1)

#make tree
contMap <- contMap(pruned.tree.forfig, abundance, res=100, plot=FALSE)
contMap <- setMap(contMap, viridisLite::viridis(n=8))
plot(contMap)

##test log value of abundance###
abundance_log <- log(abundance, base = 10)

contMap_log <- contMap(pruned.tree.forfig, abundance_log, res=100, plot=FALSE)
contMap_log <- setMap(contMap_log, viridisLite::viridis(n=8))
plot(contMap_log)

#range size phylogeny-------------
##make community data matrix#### 
allsitesmatrix_rs <- read.table("comm_phylo_analyses/Phylogenetic_signal/rangesize_comm_matrix.txt", 
                             sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree.forfig.rs <- treedata(SBtree, unlist(allsitesmatrix_rs[1,allsitesmatrix_rs[1,]>0]),warnings = F)$phy
plot(pruned.tree.forfig.rs)
is.rooted(pruned.tree.rs)

#trait data
rangesize_df <- read.csv("comm_phylo_analyses/Phylogenetic_signal/range_size_results.csv")
colnames(rangesize_df)[colnames(rangesize_df) == 'AOO..km2.'] <- 'AOO'
rangesize_df <- rangesize_df %>% 
  rename(trait_value = AOO)
rangesize_df <- rangesize_df %>% remove_rownames %>% column_to_rownames(var="Species")
rangesize_df <- df2vec(rangesize_df, colID=1)


##test log value of abundance###
rangesize_log <- log(rangesize_df, base = 10)

contMap_log_rs <- contMap(pruned.tree.forfig.rs, rangesize_log, res=100, plot=FALSE)
contMap_log_rs <- setMap(contMap_log_rs, viridisLite::viridis(n=8))
plot(contMap_log_rs)

#combine phylogenies
patchwork <- contMap_log_rs + contMap_log
