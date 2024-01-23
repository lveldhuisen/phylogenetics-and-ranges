install.packages("broom")
install.packages("janitor")

library(tidyverse)
library(dplyr)
library(picante)
library(geiger)
library(ape)
library(vegan)
library(forcats)
library(broom)
library(janitor)

##import S&B phylogeny###
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM---------------------
##make community data matrix#### 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
PBM_range_matrix <- read.table("PBMrangesize_comm_matrix.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_range_matrix[31,PBM_range_matrix[31,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

##make distance matrix for MPD and MNTD calculations####
dist.mat <- cophenetic(pruned.tree)

#phylogenetic analyses-----------

##PD calculation###
PD_PBM_range <- ses.pd(PBM_range_matrix, pruned.tree, null.model = c("sample.pool"),
       runs = 5000, include.root=TRUE)

df_PD = subset(PD_PBM_range, select = -c(pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,pd.obs.p,runs) )

plot(df$ntaxa,df$pd.obs.z)
ggplot() + geom_point(data = df_PD, mapping = aes(x=ntaxa, y=pd.obs.z)) + xlab("Number of species") + ylab("Standard effect size PD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.15, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

##MPD calculation####
MPD_PBM_range <- ses.mpd(PBM_range_matrix, dist.mat, null.model = c("sample.pool"),
        abundance.weighted = FALSE, runs = 5000, iterations = 5000)

df_MPD = subset(MPD_PBM_range, select = -c(mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,mpd.obs.p,runs) )

ggplot() + geom_point(data = df_MPD, mapping = aes(x=ntaxa, y=mpd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MPD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.15, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

##MNTD calculation####
MNTD_PBM_range <- ses.mntd(PBM_range_matrix, dist.mat, null.model = c("sample.pool"),
         abundance.weighted=FALSE, runs = 5000, iterations = 5000)

df_MNTD = subset(MNTD_PBM_range, select = -c(mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,mntd.obs.p,runs) )

ggplot() + geom_point(data = df_MNTD, mapping = aes(x=ntaxa, y=mntd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MNTD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.3, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

#Pfeiler--------------------
##make community data matrix#### 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
Pfeiler_range_matrix <- read.table("Pfeilerrangesize_comm_matrix.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Pfeiler_range_matrix[26,Pfeiler_range_matrix[26,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

##PD calculation###
PD_Pfeiler_range <- ses.pd(Pfeiler_range_matrix, pruned.tree, null.model = c("sample.pool"),
                       runs = 5000, include.root=TRUE)

df_PD = subset(PD_Pfeiler_range, select = -c(pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,pd.obs.p,runs) )

plot(df$ntaxa,df$pd.obs.z)
ggplot() + geom_point(data = df_PD, mapping = aes(x=ntaxa, y=pd.obs.z)) + xlab("Number of species") + ylab("Standard effect size PD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

##MPD calculations#####
MPD_Pfeiler_range <- ses.mpd(Pfeiler_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                         abundance.weighted = FALSE, runs = 5000, iterations = 5000)

df_MPD = subset(MPD_Pfeiler_range, select = -c(mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,mpd.obs.p,runs) )

ggplot() + geom_point(data = df_MPD, mapping = aes(x=ntaxa, y=mpd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MPD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1,2) +theme_classic() + geom_hline(yintercept = 1.15, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,30)

##MNTD calculation####
MNTD_Pfeiler_range <- ses.mntd(Pfeiler_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                           abundance.weighted=FALSE, runs = 5000, iterations = 5000)

df_MNTD = subset(MNTD_Pfeiler_range, select = -c(mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,mntd.obs.p,runs) )

ggplot() + geom_point(data = df_MNTD, mapping = aes(x=ntaxa, y=mntd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MNTD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-1.5,2) +theme_classic() + geom_hline(yintercept = 1.3, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,30)

#Road---------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
Road_range_matrix <- read.table("Roadrangesize_comm_matrix.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Road_range_matrix[32,Road_range_matrix[32,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

##PD calculation###
PD_Road_range <- ses.pd(Road_range_matrix, pruned.tree, null.model = c("sample.pool"),
                           runs = 5000, include.root=TRUE)

df_PD = subset(PD_Road_range, select = -c(pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,pd.obs.p,runs) )

plot(df$ntaxa,df$pd.obs.z)
ggplot() + geom_point(data = df_PD, mapping = aes(x=ntaxa, y=pd.obs.z)) + xlab("Number of species") + ylab("Standard effect size PD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-3,1) +theme_classic() + geom_hline(yintercept = 1.2, col = "lightgrey") + geom_hline(yintercept = -1.45, col = "lightgrey") + xlim(0,32)

##MPD calculations#####
MPD_Road_range <- ses.mpd(Road_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                             abundance.weighted = FALSE, runs = 5000, iterations = 5000)

df_MPD = subset(MPD_Road_range, select = -c(mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,mpd.obs.p,runs) )

ggplot() + geom_point(data = df_MPD, mapping = aes(x=ntaxa, y=mpd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MPD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-3,1) +theme_classic() + geom_hline(yintercept = 1.15, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)

##MNTD calculation####
MNTD_Road_range <- ses.mntd(Road_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                               abundance.weighted=FALSE, runs = 5000, iterations = 5000)

df_MNTD = subset(MNTD_Road_range, select = -c(mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,mntd.obs.p,runs) )

ggplot() + geom_point(data = df_MNTD, mapping = aes(x=ntaxa, y=mntd.obs.z)) + xlab("Number of species") + ylab("Standard effect size MNTD") +
  geom_hline(yintercept = 0, col = "darkgrey") + ylim(-2.5,1) +theme_classic() + geom_hline(yintercept = 1.3, col = "lightgrey") + geom_hline(yintercept = -1.32, col = "lightgrey") + xlim(0,32)
