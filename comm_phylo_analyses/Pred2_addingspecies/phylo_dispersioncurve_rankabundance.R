library(tidyverse)
library(dplyr)
library(picante)
library(geiger)
library(ape)
library(vegan)
library(forcats)
library(broom)
library(janitor)

#calculate phylogenetic dispersion curves based on rank abundance-------

##import S&B phylogeny###
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM---------
##make community data matrix#### 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
PBM_abun_matrix <- read.table("PBM_RA_comm_matrix.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_abun_matrix[30,PBM_abun_matrix[30,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_PBM_RA <- ses.pd(PBM_abun_matrix, pruned.tree, null.model = c("sample.pool"),
                       runs = 5000, include.root=TRUE)

df_PD = subset(PD_PBM_RA, select = -c(pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,pd.obs.p,runs) )

ggplot() + geom_point(data = df_PD, mapping = aes(x=ntaxa, y=pd.obs.z)) + xlab("Number of species") + ylab("Standard effect size PD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

###MPD####
MPD_PBM_RA <- ses.mpd(PBM_abun_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                         abundance.weighted = FALSE, runs = 5000, iterations = 5000)

df_MPD = subset(MPD_PBM_RA, select = -c(mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,mpd.obs.p,runs) )

ggplot() + geom_point(data = df_MPD, mapping = aes(x=ntaxa, y=mpd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MPD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

###MNTD###
MNTD_PBM_RA <- ses.mntd(PBM_abun_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                           abundance.weighted=FALSE, runs = 5000, iterations = 5000)

df_MNTD = subset(MNTD_PBM_RA, select = -c(mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,mntd.obs.p,runs) )

ggplot() + geom_point(data = df_MNTD, mapping = aes(x=ntaxa, y=mntd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MNTD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.3, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

##Pfeiler#####
##make community data matrix#### 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
Pfeiler_abun_matrix <- read.table("Pfeiler_RA_comm_matrix.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Pfeiler_abun_matrix[26,Pfeiler_abun_matrix[26,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Pfeiler_RA <- ses.pd(Pfeiler_abun_matrix, pruned.tree, null.model = c("sample.pool"),
                    runs = 5000, include.root=TRUE)

df_PD = subset(PD_Pfeiler_RA, select = -c(pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,pd.obs.p,runs) )

ggplot() + geom_point(data = df_PD, mapping = aes(x=ntaxa, y=pd.obs.z)) + xlab("Number of species") + ylab("Standard effect size PD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)


###MPD####
MPD_Pfeiler_RA <- ses.mpd(Pfeiler_abun_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                      abundance.weighted = FALSE, runs = 5000, iterations = 5000)

df_MPD = subset(MPD_Pfeiler_RA, select = -c(mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,mpd.obs.p,runs) )

ggplot() + geom_point(data = df_MPD, mapping = aes(x=ntaxa, y=mpd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MPD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

###MNTD###
MNTD_Pfeiler_RA <- ses.mntd(Pfeiler_abun_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                        abundance.weighted=FALSE, runs = 5000, iterations = 5000)

df_MNTD = subset(MNTD_Pfeiler_RA, select = -c(mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,mntd.obs.p,runs) )

ggplot() + geom_point(data = df_MNTD, mapping = aes(x=ntaxa, y=mntd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MNTD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.3, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

##Road#####
##make community data matrix#### 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
Road_abun_matrix <- read.table("Road_RA_comm_matrix.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Road_abun_matrix[32,Road_abun_matrix[32,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Road_RA <- ses.pd(Road_abun_matrix, pruned.tree, null.model = c("sample.pool"),
                        runs = 5000, include.root=TRUE)

df_PD = subset(PD_Road_RA, select = -c(pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,pd.obs.p,runs) )

ggplot() + geom_point(data = df_PD, mapping = aes(x=ntaxa, y=pd.obs.z)) + xlab("Number of species") + ylab("Standard effect size PD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

###MPD####
MPD_Road_RA <- ses.mpd(Road_abun_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                          abundance.weighted = FALSE, runs = 5000, iterations = 5000)

df_MPD = subset(MPD_Road_RA, select = -c(mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,mpd.obs.p,runs) )

ggplot() + geom_point(data = df_MPD, mapping = aes(x=ntaxa, y=mpd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MPD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1.5,2) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.25, col = "lightgrey") + xlim(0,32)

###MNTD###

MNTD_Road_RA <- ses.mntd(Road_abun_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                            abundance.weighted=FALSE, runs = 5000, iterations = 5000)

df_MNTD = subset(MNTD_Road_RA, select = -c(mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,mntd.obs.p,runs) )

ggplot() + geom_point(data = df_MNTD, mapping = aes(x=ntaxa, y=mntd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MNTD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.3, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)
