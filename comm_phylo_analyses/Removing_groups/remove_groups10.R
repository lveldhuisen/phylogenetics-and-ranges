#removing groups of 10 species by abundance 

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

#import S&B phylogeny--------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM------------
##make community data matrix#### 
PBM_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/PBM_removing10_abundance.txt", 
                                sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(PBM_groups_matrix[24,PBM_groups_matrix[24,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_PBM_group_removal <- ses.pd(PBM_groups_matrix, pruned.tree, 
                               null.model = c("sample.pool"),
                                   runs = 5000, include.root=TRUE) #all PBM PD SES is 0.89


PD_PBM_group_removal <- PD_PBM_group_removal[-c(23,24),]
PD_PBM_group_removal$Group_removed <- c(1:22)

PD_PBM_group_removal_fig <- ggplot(data= PD_PBM_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=0.89, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=pd.obs.z), size = 2) +
  xlab("Group of species removed") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(0, 1, 0.5, 1, 1.5,2),limits=c(0, 2))+
  theme_classic(14) +
  geom_hline(yintercept = 0.89, col = "lightgrey") +
  xlim(0,22) 
plot(PD_PBM_group_removal_fig)


#Pfeiler-----------
###PD#####
##make community data matrix#### 
PBM_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/PBM_removing10_abundance.txt", 
                                sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(PBM_groups_matrix[24,PBM_groups_matrix[24,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_PBM_group_removal <- ses.pd(PBM_groups_matrix, pruned.tree, 
                               null.model = c("sample.pool"),
                               runs = 5000, include.root=TRUE) #all PBM PD SES is 0.89


PD_PBM_group_removal <- PD_PBM_group_removal[-c(23,24),]
PD_PBM_group_removal$Group_removed <- c(1:22)

PD_PBM_group_removal_fig <- ggplot(data= PD_PBM_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=0.89, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=pd.obs.z), size = 2) +
  xlab("Group of species removed") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(0, 1, 0.5, 1, 1.5,2),limits=c(0, 2))+
  theme_classic(14) +
  geom_hline(yintercept = 0.89, col = "lightgrey") +
  xlim(0,22) 
plot(PD_PBM_group_removal_fig)
