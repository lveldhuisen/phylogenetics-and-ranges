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

#MPD---------------
#PBM------------
##make community data matrix#### 
PBM_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/PBM_removing10_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(PBM_groups_matrix[29,PBM_groups_matrix[29,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###MPD#####
MPD_PBM_group_removal <- ses.mpd(PBM_groups_matrix, cophenetic(pruned.tree), 
                               null.model = c("sample.pool"),abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all PBM MPD SES is 0.84


MPD_PBM_group_removal <- MPD_PBM_group_removal[-c(28,29),]#remove all and all pbm rows
MPD_PBM_group_removal$Group_removed <- c(1:27) #add column for groups
MPD_PBM_group_removal$Site <- c("High elevation (3380 m)") #add column for site 
MPD_PBM_group_removal$Baseline_MPD <- c(0.84) #add column for baseline PD for grey line in fig
MPD_PBM_group_removal = subset(MPD_PBM_group_removal, 
                              select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs)) #remove extra columns
MPD_PBM_group_removal <- MPD_PBM_group_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 


#make individual figure 
MPD_PBM_group_removal_fig <- ggplot(data= MPD_PBM_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=0.84, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(-0.5,0, 1, 0.5, 1, 1.5,2),limits=c(-0.7, 2))+
  theme_classic(14) +
  geom_hline(yintercept = 0.84, col = "lightgrey") +
  xlim(0,27) 
plot(MPD_PBM_group_removal_fig)


#Pfeiler-----------
###MPD#####
##make community data matrix#### 
Pfeiler_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/Pfeiler_group10_matrix_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(Pfeiler_groups_matrix[22,Pfeiler_groups_matrix[22,]>0]), 
                        warnings = F)$phy

write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

mpd_pfeiler_tree <- cophenetic(pruned.tree)

###MPD#####
MPD_Pfeiler_group_removal <- ses.mpd(Pfeiler_groups_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted = FALSE, runs = 5000, iterations = 5000)
#all Pfeiler MPD SES is -0.51

MPD_Pfeiler_group_removal <- MPD_Pfeiler_group_removal[-c(21,22),]
MPD_Pfeiler_group_removal$Group_removed <- c(1:20)
MPD_Pfeiler_group_removal$Site <- c("Middle elevation (3165 m)") #add column for site 
MPD_Pfeiler_group_removal$Baseline_MPD <- c(-0.51) #add column for baseline PD for grey line in fig
MPD_Pfeiler_group_removal = subset(MPD_Pfeiler_group_removal, 
                                  select = -c(ntaxa,mpd.obs,mpd.rand.mean,
                                              mpd.rand.sd,mpd.obs.rank,runs)) #remove extra columns
MPD_Pfeiler_group_removal <- MPD_Pfeiler_group_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

#make individual figure
MPD_Pfeiler_group_removal_fig <- ggplot(data= MPD_Pfeiler_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=-0.51, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(-1,-0.5,0, 1, 0.5, 1),limits=c(-1.5, 1))+
  theme_classic(14) +
  geom_hline(yintercept = -0.51, col = "lightgrey") +
  xlim(0,20) 
plot(MPD_Pfeiler_group_removal_fig)

#Road-------------
##make community data matrix#### 
Road_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/Road_groups10_matrix_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(Road_groups_matrix[32,Road_groups_matrix[32,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###MPD#######
MPD_Road_group_removal <- ses.mpd(Road_groups_matrix, cophenetic(pruned.tree), 
                                null.model = c("sample.pool"),
                                abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all Road PD SES is -0.85


MPD_Road_group_removal <- MPD_Road_group_removal[-c(31,32),]
MPD_Road_group_removal$Group_removed <- c(1:30)
MPD_Road_group_removal$Baseline_MPD <- c(-0.85) #add column for baseline PD for grey line in fig
MPD_Road_group_removal$Site <- c("Low elevation (2815 m)") #add column for site 
MPD_Road_group_removal = subset(MPD_Road_group_removal, 
                               select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs)) #remove extra columns
MPD_Road_group_removal <- MPD_Road_group_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

#make individual figure
MPD_Road_group_removal_fig <- ggplot(data= MPD_Road_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=-0.85, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(-2,-1.5,-1,-0.5,0, 1, 0.5, 1,1.5),limits=c(-2.2, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.85, col = "lightgrey") +
  xlim(0,30) 
plot(MPD_Road_group_removal_fig)

#combine into one dataset---------
MPD_groups_allsites <- rbind(MPD_PBM_group_removal,MPD_Pfeiler_group_removal,MPD_Road_group_removal)
MPD_groups_allsites$Site = factor(MPD_groups_allsites$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

#make faceted figure
allsites_MPD_groups <- ggplot(data= MPD_groups_allsites) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, yend=SES, y=Baseline_MPD), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Abundance rank of group of removed species (most to least)") +
  ylab("Standard effect size") +
  scale_y_continuous(name="Standard effect size", breaks = c(-2,-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5,2),limits=c(-2.2, 2))+
  theme_bw(base_size = 20) +
  xlim(0,30) +
  geom_abline(data = MPD_groups_allsites, aes(intercept = Baseline_MPD, slope = 0)) +
  facet_grid(.~Site)
plot(allsites_MPD_groups)


#MNTD--------------
#PBM------------
##make community data matrix#### 
PBM_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/PBM_removing10_abundance.txt",  sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(PBM_groups_matrix[24,PBM_groups_matrix[24,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###MNTD#####
MNTD_PBM_group_removal <- ses.mntd(PBM_groups_matrix, cophenetic(pruned.tree), 
                                 null.model = c("sample.pool"),abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all PBM mntd SES is 0

MNTD_PBM_group_removal <- MNTD_PBM_group_removal[-c(28,29),]#remove all and all pbm rows
MNTD_PBM_group_removal$Group_removed <- c(1:27) #add column for groups
MNTD_PBM_group_removal$Site <- c("High elevation (3380 m)") #add column for site 
MNTD_PBM_group_removal$Baseline_MNTD <- c(0) #add column for baseline PD for grey line in fig
MNTD_PBM_group_removal = subset(MNTD_PBM_group_removal, 
                               select = -c(ntaxa,mntd.obs,mntd.rand.mean,
                                           mntd.rand.sd,mntd.obs.rank,runs)) #remove extra columns
MNTD_PBM_group_removal <- MNTD_PBM_group_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 


#make individual figure 
MNTD_PBM_group_removal_fig <- ggplot(data= MNTD_PBM_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=0, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-0.5,0, 1, 0.5, 1, 1.5,2),limits=c(-0.5, 2))+
  theme_classic(14) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  xlim(0,27) 
plot(MNTD_PBM_group_removal_fig)


#Pfeiler-----------
###MNTD#####
##make community data matrix#### 
Pfeiler_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/Pfeiler_group10_matrix_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(Pfeiler_groups_matrix[19,Pfeiler_groups_matrix[19,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###MNTD#####
MNTD_Pfeiler_group_removal <- ses.mntd(Pfeiler_groups_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted = FALSE, runs = 5000, iterations = 5000)
#all Pfeiler MPD SES is -1.35

MNTD_Pfeiler_group_removal <- MNTD_Pfeiler_group_removal[-c(21, 22),]
MNTD_Pfeiler_group_removal$Group_removed <- c(1:20)
MNTD_Pfeiler_group_removal$Site <- c("Middle elevation (3165 m)") #add column for site 
MNTD_Pfeiler_group_removal$Baseline_MNTD <- c(-1.35) #add column for baseline PD for grey line in fig
MNTD_Pfeiler_group_removal = subset(MNTD_Pfeiler_group_removal, 
                                  select = -c(ntaxa,mntd.obs,mntd.rand.mean,
                                              mntd.rand.sd,mntd.obs.rank,runs)) #remove extra columns
MNTD_Pfeiler_group_removal <- MNTD_Pfeiler_group_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

#make individual figure
MNTD_Pfeiler_group_removal_fig <- ggplot(data= MNTD_Pfeiler_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=-1.35, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-1.5,-1,-0.5,0, 1, 0.5, 1),limits=c(-1.7, 1))+
  theme_classic(14) +
  geom_hline(yintercept = -1.35, col = "lightgrey") +
  xlim(0,20) 
plot(MNTD_Pfeiler_group_removal_fig)

#Road-------------
##make community data matrix#### 
Road_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/Road_groups10_matrix_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(Road_groups_matrix[25,Road_groups_matrix[25,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###MPD#######
MNTD_Road_group_removal <- ses.mntd(Road_groups_matrix, cophenetic(pruned.tree), 
                                  null.model = c("sample.pool"),
                                  abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all Road MNTD SES is 0.32


MNTD_Road_group_removal <- MNTD_Road_group_removal[-c(31,32),]
MNTD_Road_group_removal$Group_removed <- c(1:30)
MNTD_Road_group_removal$Baseline_MNTD <- c(0.32) #add column for baseline PD for grey line in fig
MNTD_Road_group_removal$Site <- c("Low elevation (2815 m)") #add column for site 
MNTD_Road_group_removal = subset(MNTD_Road_group_removal, 
                                select = -c(ntaxa,mntd.obs,mntd.rand.mean,
                                            mntd.rand.sd,mntd.obs.rank,runs)) #remove extra columns
MNTD_Road_group_removal <- MNTD_Road_group_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

#make individual figure
MNTD_Road_group_removal_fig <- ggplot(data= MNTD_Road_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=0.32, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-2,-1.5,-1,-0.5,0, 1, 0.5, 1,1.5),limits=c(-2, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.32, col = "lightgrey") +
  xlim(0,30) 
plot(MNTD_Road_group_removal_fig)

#combine into one dataset---------
MNTD_groups_allsites <- rbind(MNTD_PBM_group_removal,MNTD_Pfeiler_group_removal,MNTD_Road_group_removal)
MNTD_groups_allsites$Site = factor(MNTD_groups_allsites$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

#make faceted figure
allsites_MNTD_groups <- ggplot(data= MNTD_groups_allsites) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, yend=SES, y=Baseline_MNTD), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Abundance rank of group of removed species (most to least)") +
  ylab("Standard effect size") +
  scale_y_continuous(name="Standard effect size", breaks = c(-2,-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5,2),limits=c(-2, 2))+
  theme_bw(base_size = 20) +
  xlim(0,30) +
  geom_abline(data = MNTD_groups_allsites, aes(intercept = Baseline_MNTD, slope = 0)) +
  facet_grid(.~Site)
plot(allsites_MNTD_groups)

#combine MPD and MNTD figures with patchwork------------
MPD_MNTD_fig <- allsites_MPD_groups / allsites_MNTD_groups + 
  plot_annotation(tag_levels = c('A'), tag_suffix = ')')+
  plot_layout(guides = 'collect')
plot(MPD_MNTD_fig)
