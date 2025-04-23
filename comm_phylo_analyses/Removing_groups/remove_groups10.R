#removing groups of 10 species by abundance to test Faith's PD

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
PBM_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/PBM_removing10_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(PBM_groups_matrix[29,PBM_groups_matrix[29,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_PBM_group_removal <- ses.pd(PBM_groups_matrix, pruned.tree,null.model = c("sample.pool"),runs = 5000, include.root=TRUE) #all PBM PD SES is 0.62


PD_PBM_group_removal <- PD_PBM_group_removal[-c(28,29),]#remove all and all pbm rows
PD_PBM_group_removal$Group_removed <- c(1:27) #add column for groups
PD_PBM_group_removal$Site <- c("High elevation (3380 m)") #add column for site 
PD_PBM_group_removal$Baseline_PD <- c(0.62) #add column for baseline PD for grey line in fig
PD_PBM_group_removal = subset(PD_PBM_group_removal, 
                              select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs)) #remove extra columns
PD_PBM_group_removal <- PD_PBM_group_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 


#make individual figure 
PD_PBM_group_removal_fig <- ggplot(data= PD_PBM_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=0.62, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(0, 1, 0.5, 1, 1.5,2),limits=c(-0.2, 2))+
  theme_classic(14) +
  geom_hline(yintercept = 0.62, col = "lightgrey") +
  xlim(0,27) 
plot(PD_PBM_group_removal_fig)


#Pfeiler-----------
##make community data matrix#### 
Pfeiler_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/Pfeiler_group10_matrix_abundance.txt", sep = "\t", header = TRUE, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree,unlist(Pfeiler_groups_matrix[22,Pfeiler_groups_matrix[22,]>0]), 
                        warnings = F)$phy
#check tree
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Pfeiler_group_removal <- ses.pd(Pfeiler_groups_matrix, pruned.tree, 
                               null.model = c("sample.pool"),
                               runs = 5000, include.root=TRUE) #all Pfeiler PD SES is -1.3


PD_Pfeiler_group_removal <- PD_Pfeiler_group_removal[-c(21,22),] #remove 'all' lines
PD_Pfeiler_group_removal$Group_removed <- c(1:20) #number in order of abundance
PD_Pfeiler_group_removal$Site <- c("Middle elevation (3165 m)") #add column for site 
PD_Pfeiler_group_removal$Baseline_PD <- c(-1.3) #add column for baseline PD for grey line in fig
PD_Pfeiler_group_removal = subset(PD_Pfeiler_group_removal, 
                              select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs)) #remove extra columns
PD_Pfeiler_group_removal <- PD_Pfeiler_group_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

#make individual figure
PD_Pfeiler_group_removal_fig <- ggplot(data= PD_Pfeiler_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=-1.3, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5,-1,-0.5,0, 1, 0.5, 1),limits=c(-1.7, 1))+
  theme_classic(14) +
  geom_hline(yintercept = -1.3, col = "lightgrey") +
  xlim(0,20) 

plot(PD_Pfeiler_group_removal_fig)

#Road-------------
##make community data matrix#### 
Road_groups_matrix <- read.table("comm_phylo_analyses/Removing_groups/comm_matrices/Road_groups10_matrix_abundance.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, 
                        unlist(Road_groups_matrix[32,Road_groups_matrix[32,]>0]), 
                        warnings = F)$phy
#check tree
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#######
PD_Road_group_removal <- ses.pd(Road_groups_matrix, pruned.tree, 
                                   null.model = c("sample.pool"),
                                   runs = 5000, include.root=TRUE) #all Road PD SES is -0.03


PD_Road_group_removal <- PD_Road_group_removal[-c(31,32),] #remove 'all' lines
PD_Road_group_removal$Group_removed <- c(1:30) #number in order of abundance
PD_Road_group_removal$Baseline_PD <- c(-0.03) #add column for baseline PD for grey line in fig
PD_Road_group_removal$Site <- c("Low elevation (2815 m)") #add column for site 
PD_Road_group_removal = subset(PD_Road_group_removal, 
                                  select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs)) #remove extra columns
PD_Road_group_removal <- PD_Road_group_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

#make individual figure
PD_Road_group_removal_fig <- ggplot(data= PD_Road_group_removal) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, y=-0.03, yend=SES), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Group of species removed") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-2,-1.5,-1,-0.5,0, 1, 0.5, 1,1.5),limits=c(-2, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.03, col = "lightgrey") +
  xlim(0,30) 

plot(PD_Road_group_removal_fig)

#combine into one dataset---------

#combine
PD_groups_allsites <- rbind(PD_PBM_group_removal,PD_Pfeiler_group_removal,PD_Road_group_removal)

#order sites by elevation
PD_groups_allsites$Site = factor(PD_groups_allsites$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

#make faceted figure
allsites_PD_groups <- ggplot(data= PD_groups_allsites) + 
  geom_segment( aes(x=Group_removed, xend=Group_removed, yend=SES, y=Baseline_PD), color="grey")+
  geom_point(mapping = aes(x=Group_removed, y=SES), size = 2) +
  xlab("Abundance rank of group of removed species (most to least)") +
  ylab("Standard effect size") +
  scale_y_continuous(name="Standard effect size", breaks = c(-2,-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5,2),limits=c(-2, 2))+
  theme_bw(base_size = 20) +
  xlim(0,30) +
  geom_abline(data = PD_groups_allsites, aes(intercept = Baseline_PD, slope = 0)) +
  facet_grid(.~Site)

plot(allsites_PD_groups)
ggsave("figures/AJB_final/Fig5.png", dpi = 600, width = 14.5, height = 6)
