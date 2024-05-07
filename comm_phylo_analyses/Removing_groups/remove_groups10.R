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


PD_PBM_group_removal <- PD_PBM_group_removal[-c(23,24),]#remove all and all pbm rows
PD_PBM_group_removal$Group_removed <- c(1:22) #add column for groups
PD_PBM_group_removal$Site <- c("High elevation (3380 m)") #add column for site 
PD_PBM_group_removal$Baseline_PD <- c(0.89) #add column for baseline PD for grey line in fig
PD_PBM_group_removal = subset(PD_PBM_group_removal, 
                              select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs)) #remove extra columns
PD_PBM_group_removal <- PD_PBM_group_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 


#make individual figure 
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
Pfeiler_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/Pfeiler_group10_matrix_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(Pfeiler_groups_matrix[19,Pfeiler_groups_matrix[19,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Pfeiler_group_removal <- ses.pd(Pfeiler_groups_matrix, pruned.tree, 
                               null.model = c("sample.pool"),
                               runs = 5000, include.root=TRUE) #all Pfeiler PD SES is -0.36


PD_Pfeiler_group_removal <- PD_Pfeiler_group_removal[-c(18,19),]
PD_Pfeiler_group_removal$Group_removed <- c(1:17)
PD_Pfeiler_group_removal$Site <- c("Middle elevation (3165 m)") #add column for site 
PD_Pfeiler_group_removal$Baseline_PD <- c(-0.36) #add column for baseline PD for grey line in fig
PD_Pfeiler_group_removal = subset(PD_Pfeiler_group_removal, 
                              select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs)) #remove extra columns
PD_Pfeiler_group_removal <- PD_Pfeiler_group_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

#make individual figure
PD_Pfeiler_group_removal_fig <- ggplot(data= PD_Pfeiler_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=-0.36, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=pd.obs.z), size = 2) +
  xlab("Group of species removed") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1,-0.5,0, 1, 0.5, 1),limits=c(-1, 1))+
  theme_classic(14) +
  geom_hline(yintercept = -0.36, col = "lightgrey") +
  xlim(0,18) 
plot(PD_Pfeiler_group_removal_fig)

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

###PD#######
PD_Road_group_removal <- ses.pd(Road_groups_matrix, pruned.tree, 
                                   null.model = c("sample.pool"),
                                   runs = 5000, include.root=TRUE) #all Road PD SES is -0.17


PD_Road_group_removal <- PD_Road_group_removal[-c(24,25),]
PD_Road_group_removal$Group_removed <- c(1:23)
PD_Road_group_removal$Baseline_PD <- c(-0.17) #add column for baseline PD for grey line in fig
PD_Road_group_removal$Site <- c("Low elevation (2815 m)") #add column for site 
PD_Road_group_removal = subset(PD_Road_group_removal, 
                                  select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs)) #remove extra columns
PD_Road_group_removal <- PD_Road_group_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

#make individual figure
PD_Road_group_removal_fig <- ggplot(data= PD_Road_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=-0.17, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-2,-1.5,-1,-0.5,0, 1, 0.5, 1,1.5),limits=c(-2, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.17, col = "lightgrey") +
  xlim(0,23) 
plot(PD_Road_group_removal_fig)

#combine into one dataset---------
PD_groups_allsites <- rbind(PD_PBM_group_removal,PD_Pfeiler_group_removal,PD_Road_group_removal)
PD_groups_allsites$Site = factor(PD_groups_allsites$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

#make faceted figure
allsites_PD_groups <- ggplot(data= PD_groups_allsites) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, yend=SES, y=Baseline_PD), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Abundance rank of group of removed species (most to least)") +
  ylab("Standard effect size") +
  scale_y_continuous(name="Standard effect size", breaks = c(-2,-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5,2),limits=c(-2, 2))+
  theme_bw(14) +
  xlim(0,23) +
  geom_abline(data = PD_groups_allsites, aes(intercept = Baseline_PD, slope = 0)) +
  facet_grid(.~Site)
plot(allsites_PD_groups)
