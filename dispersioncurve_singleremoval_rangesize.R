#code to generate curves of phylogenetic dispersion while removing 1 species at a time 
#will result in an individual impact on PD/MPD/MNTD for each species 

library(tidyverse)
library(dplyr)
library(picante)
library(geiger)
library(ape)
library(vegan)
library(forcats)
library(broom)
library(janitor)

#import S&B phylogeny#
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM----------
##PD####
##make community data matrix#### 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses/Removing1species_atatime")
PBM_range_matrix <- read.table("PBMrangesize_comm_matrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_range_matrix[33,PBM_range_matrix[33,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_PBM_range_removal <- ses.pd(PBM_range_matrix, pruned.tree, null.model = c("sample.pool"),
                    runs = 5000, include.root=TRUE)

PD_PBM_range_removal$Range_size_rank <- c(1:33)
#all PBM PD is SES is 0.36
PD_PBM_range_removal <- PD_PBM_range_removal[-c(32,33),]

PD_PBM_rangesize_removal <- ggplot(data= PD_PBM_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0.36, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size PD") +
  ylim(-1,1) +
  theme_classic(14) +
  geom_hline(yintercept = 0.36, col = "lightgrey") +
  xlim(0,32) 



##MPD###
##MNTD####
#Pfeiler--------
##PD####
##MPD###
##MNTD####
#Road-------
##PD####
##MPD###
##MNTD####