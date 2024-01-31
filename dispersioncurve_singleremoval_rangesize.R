#code to generate curves of phylogenetic dispersion while removing 1 species at a time, ranked by decreasing range size 
#will result in an individual impact on PD/MPD/MNTD for each species 

install.packages("patchwork")

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
library(ggplot2)

#import S&B phylogeny#
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM----------
##make community data matrix#### 
PBM_range_matrix <- read.table("comm_phylo_analyses/Removing1species_atatime/PBMrangesize_comm_matrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_range_matrix[33,PBM_range_matrix[33,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

co_matrix <- cophenetic(pruned.tree)

###PD#####
PD_PBM_range_removal <- ses.pd(PBM_range_matrix, pruned.tree, null.model = c("sample.pool"),
                    runs = 5000, include.root=TRUE)

PD_PBM_range_removal$Range_size_rank <- c(1:33)
#all PBM PD is SES is 0.91
PD_PBM_range_removal <- PD_PBM_range_removal[-c(32,33),]

PD_PBM_rangesize_removal_fig <- ggplot(data= PD_PBM_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0.90, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size PD") +
  ylim(-0.5,2) +
  theme_classic(14) +
  geom_hline(yintercept = 0.9, col = "lightgrey") +
  xlim(0,32) 
plot(PD_PBM_rangesize_removal_fig)


###MPD#####
MPD_PBM_range_removal <- ses.mpd(PBM_range_matrix, co_matrix, null.model = c("sample.pool"), 
        abundance.weighted = FALSE, runs = 5000, iterations = 5000)

#all PBM MPD is SES is 1.25
MPD_PBM_range_removal <- MPD_PBM_range_removal[-c(32,33),]
MPD_PBM_range_removal$Range_size_rank <- c(1:31)

MPD_PBM_range_removal_fig <- ggplot(data= MPD_PBM_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=1.25, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mpd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size MPD") +
  ylim(0,2) +
  theme_classic(14) +
  geom_hline(yintercept = 1.25, col = "lightgrey") +
  xlim(0,32) 

###MNTD####
MNTD_PBM_range_removal <- ses.mntd(PBM_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                        abundance.weighted=FALSE, runs = 5000, iterations = 5000)
#all PBM MNTD is SES is 0.3
MNTD_PBM_range_removal <- MNTD_PBM_range_removal[-c(32,33),]
MNTD_PBM_range_removal$Range_size_rank <- c(1:31)

MNTD_PBM_range_removal_fig <- ggplot(data= MNTD_PBM_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0.3, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mntd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size MNTD") +
  ylim(-1,1) +
  theme_classic(14) +
  geom_hline(yintercept = 0.3, col = "lightgrey") +
  xlim(0,32) 
plot(MNTD_PBM_range_removal_fig)

#Pfeiler--------
Pfeiler_range_matrix <- read.table("comm_phylo_analyses/Removing1species_atatime/Pfeiler_ranges_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Pfeiler_range_matrix[28,Pfeiler_range_matrix[28,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Pfeiler_range_removal <- ses.pd(Pfeiler_range_matrix, pruned.tree, null.model = c("sample.pool"),
                               runs = 5000, include.root=TRUE)
PD_Pfeiler_range_removal <- PD_Pfeiler_range_removal[-c(27),]
PD_Pfeiler_range_removal$Range_size_rank <- c(1:27)
#all Pfeiler  PD is SES is -0.6


PD_Pfeiler_rangesize_removal_fig <- ggplot(data= PD_Pfeiler_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-0.6, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size PD") +
  ylim(-1,0.5) +
  theme_classic(14) +
  geom_hline(yintercept = -0.6, col = "lightgrey") +
  xlim(0,28) 

plot(PD_Pfeiler_rangesize_removal_fig)

###MPD#####
MPD_Pfeiler_range_removal <- ses.mpd(Pfeiler_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), 
                                 abundance.weighted = FALSE, runs = 5000, iterations = 5000)

#all Pfeiler MPD is SES is -0.1
MPD_Pfeiler_range_removal <- MPD_Pfeiler_range_removal[-c(27),]
MPD_Pfeiler_range_removal$Range_size_rank <- c(1:26)

MPD_Pfeiler_range_removal_fig <- ggplot(data= MPD_Pfeiler_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-0.1, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mpd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size MPD") +
  ylim(-1,1) +
  theme_classic(14) +
  geom_hline(yintercept = -0.1, col = "lightgrey") +
  xlim(0,28) 
plot(MPD_Pfeiler_range_removal_fig)

###MNTD####
MNTD_Pfeiler_range_removal <- ses.mntd(Pfeiler_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                                   abundance.weighted=FALSE, runs = 5000, iterations = 5000)
#all Pfeiler MNTD is SES is -1.06
MNTD_Pfeiler_range_removal <- MNTD_Pfeiler_range_removal[-c(27:29),]
MNTD_Pfeiler_range_removal$Range_size_rank <- c(1:26)

MNTD_Pfeiler_range_removal_fig <- ggplot(data= MNTD_Pfeiler_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-1.06, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mntd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size MNTD") +
  ylim(-2,1) +
  theme_classic(14) +
  geom_hline(yintercept = -1.06, col = "lightgrey") +
  xlim(0,28) 
plot(MNTD_Pfeiler_range_removal_fig)

#Road-------
Road_range_matrix <- read.table("comm_phylo_analyses/Removing1species_atatime/Road_range_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Road_range_matrix[34,Road_range_matrix[34,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Road_range_removal <- ses.pd(Road_range_matrix, pruned.tree, null.model = c("sample.pool"),
                                   runs = 5000, include.root=TRUE)

PD_Road_range_removal <- PD_Road_range_removal[-c(33,34),]
PD_Road_range_removal$Range_size_rank <- c(1:32)
#all Road  PD is SES is 0


PD_Road_rangesize_removal_fig <- ggplot(data= PD_Road_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size PD") +
  ylim(-1,0.5) +
  theme_classic(14) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  xlim(0,33) 

plot(PD_Road_rangesize_removal_fig)

###MPD####
MPD_Road_range_removal <- ses.mpd(Road_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), 
                                     abundance.weighted = FALSE, runs = 5000, iterations = 5000)

#all Road MPD is SES is -0.89
MPD_Road_range_removal <- MPD_Road_range_removal[-c(33,34),]
MPD_Road_range_removal$Range_size_rank <- c(1:32)

MPD_Road_range_removal_fig <- ggplot(data= MPD_Road_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-0.89, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mpd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size MPD") +
  ylim(-1.5,0.5) +
  theme_classic(14) +
  geom_hline(yintercept = -0.89, col = "lightgrey") +
  xlim(0,32) 
plot(MPD_Road_range_removal_fig)

###MNTD#####
MNTD_Road_range_removal <- ses.mntd(Road_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                                       abundance.weighted=FALSE, runs = 5000, iterations = 5000)
#all Road MNTD is SES is 0.15
MNTD_Road_range_removal <- MNTD_Road_range_removal[-c(33,34),]
MNTD_Road_range_removal$Range_size_rank <- c(1:32)

MNTD_Road_range_removal_fig <- ggplot(data= MNTD_Road_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0.15, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mntd.obs.z), size = 2) +
  xlab("Range size rank of removed species") +
  ylab("Standard effect size MNTD") +
  ylim(-1,1) +
  theme_classic(14) +
  geom_hline(yintercept = 0.15, col = "lightgrey") +
  xlim(0,32) 
plot(MNTD_Road_range_removal_fig)

#make combined fig with all sites and metrics-----------------
PBM_fig <- (PD_PBM_rangesize_removal_fig + MPD_PBM_range_removal_fig + MNTD_PBM_range_removal_fig)+ plot_layout(axis_titles = "collect")
plot(PBM_fig)

Pfeiler_fig <- (PD_Pfeiler_rangesize_removal_fig+MPD_Pfeiler_range_removal_fig+MNTD_Pfeiler_range_removal_fig)+ plot_layout(axis_titles = "collect")
plot(Pfeiler_fig)

Road_fig <- (PD_Road_rangesize_removal_fig+MPD_Road_range_removal_fig+MNTD_Road_range_removal_fig) + plot_layout(axis_titles = "collect")

PBM_fig/Pfeiler_fig/Road_fig

PD_fig <- (PD_PBM_rangesize_removal_fig | PD_Pfeiler_rangesize_removal_fig|PD_Road_rangesize_removal_fig) + plot_layout(axis_titles = "collect")
plot(PD_fig)

MPD_fig <- (MPD_PBM_range_removal_fig|MPD_Pfeiler_range_removal_fig|MPD_Road_range_removal_fig) + plot_layout(axis_titles = "collect")
plot(MPD_fig)

MNTD_fig <- (MNTD_PBM_range_removal_fig|MNTD_Pfeiler_range_removal_fig|MNTD_Road_range_removal_fig) + plot_layout(axis_titles = "collect")
plot(MNTD_fig)

(PD_fig/MPD_fig/MNTD_fig) + plot_layout(axis_titles = "collect")
