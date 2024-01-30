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

#code to generate curves of phylogenetic dispersion while removing 1 species at a time, ranked by decreasing abundance 

#import S&B phylogeny--------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM------------
##make community data matrix#### 
PBM_abundance_matrix <- read.table("comm_phylo_analyses/Removing1species_atatime/PBM_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_abundance_matrix[32,PBM_abundance_matrix[32,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_PBM_abundance_removal <- ses.pd(PBM_abundance_matrix, pruned.tree, null.model = c("sample.pool"),
                               runs = 5000, include.root=TRUE) #all PBM PD is SES is 0.57

PD_PBM_abundance_removal$Abundance_rank <- c(1:31)
PD_PBM_abundance_removal <- PD_PBM_abundance_removal[-c(32,33),]

PD_PBM_abundance_removal_fig <- ggplot(data= PD_PBM_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.57, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species") +
  ylab("Standard effect size PD") +
  ylim(-0.5,1) +
  theme_classic(14) +
  geom_hline(yintercept = 0.57, col = "lightgrey") +
  xlim(0,32) 
plot(PD_PBM_abundance_removal_fig)

###MPD########
MPD_PBM_abundance_removal <- ses.mpd(PBM_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), 
                                 abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all PBM MPD is SES is 0.9

MPD_PBM_abundance_removal <- MPD_PBM_abundance_removal[-c(32,33),]
MPD_PBM_abundance_removal$Abundance_rank <- c(1:31)

MPD_PBM_abundance_removal_fig <- ggplot(data= MPD_PBM_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.9, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mpd.obs.z), size = 2) +
  xlab("Abundance rank of removed species") +
  ylab("Standard effect size MPD") +
  ylim(0,1.5) +
  theme_classic(14) +
  geom_hline(yintercept = 0.9, col = "lightgrey") +
  xlim(0,32) 
plot(MPD_PBM_abundance_removal_fig)

###MNTD######
MNTD_PBM_abundance_removal <- ses.mntd(PBM_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                                   abundance.weighted=FALSE, runs = 5000, iterations = 5000) #all PBM MNTD is SES is 0.12

MNTD_PBM_abundance_removal <- MNTD_PBM_abundance_removal[-c(32,33),]
MNTD_PBM_abundance_removal$Abundance_rank <- c(1:31)

MNTD_PBM_abundance_removal_fig <- ggplot(data= MNTD_PBM_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.12, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mntd.obs.z), size = 2) +
  xlab("Abundance rank of removed species") +
  ylab("Standard effect size MNTD") +
  ylim(-1,1) +
  theme_classic(14) +
  geom_hline(yintercept = 0.12, col = "lightgrey") +
  xlim(0,32) 
plot(MNTD_PBM_abundance_removal_fig)

#Pfeiler--------
###PD#####
###MPD########
###MNTD######
#Road--------------
###PD#####
###MPD########
###MNTD######