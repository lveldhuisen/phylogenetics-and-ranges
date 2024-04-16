install.packages("car")

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

#code to generate curves of phylogenetic dispersion while removing 1 species at a time, ranked by decreasing abundance 

#import S&B phylogeny--------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM------------
##make community data matrix#### 
PBM_abundance_matrix <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/PBM_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_abundance_matrix[32,PBM_abundance_matrix[32,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_PBM_abundance_removal <- ses.pd(PBM_abundance_matrix, pruned.tree, null.model = c("sample.pool"),
                               runs = 5000, include.root=TRUE) #all PBM PD is SES is 0.57

PD_PBM_abundance_removal <- PD_PBM_abundance_removal[-c(32,33),]
PD_PBM_abundance_removal$Abundance_rank <- c(1:31)

PD_PBM_abundance_removal_fig <- ggplot(data= PD_PBM_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.57, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
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
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5)) +
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
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.12, col = "lightgrey") +
  xlim(0,32) 
plot(MNTD_PBM_abundance_removal_fig)

#Pfeiler--------
##make community data matrix#### 
Pfeiler_abundance_matrix <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/Pfeiler_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Pfeiler_abundance_matrix[28,Pfeiler_abundance_matrix[28,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Pfeiler_abundance_removal <- ses.pd(Pfeiler_abundance_matrix, pruned.tree, null.model = c("sample.pool"),
                                   runs = 5000, include.root=TRUE) #all Pfeiler PD is SES is -0.6

PD_Pfeiler_abundance_removal$Abundance_rank <- c(1:28)
PD_Pfeiler_abundance_removal <- PD_Pfeiler_abundance_removal[-c(27,28),]

PD_Pfeiler_abundance_removal_fig <- ggplot(data= PD_Pfeiler_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.6, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.6, col = "lightgrey") +
  xlim(0,28) 
plot(PD_Pfeiler_abundance_removal_fig)

###MPD########
MPD_Pfeiler_abundance_removal <- ses.mpd(Pfeiler_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), 
                                     abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all Pfeiler MPD is SES is 0.02

MPD_Pfeiler_abundance_removal <- MPD_Pfeiler_abundance_removal[-c(27,28),]
MPD_Pfeiler_abundance_removal$Abundance_rank <- c(1:26)

MPD_Pfeiler_abundance_removal_fig <- ggplot(data= MPD_Pfeiler_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.02, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mpd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  scale_y_continuous(name="SES MPD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.02, col = "lightgrey") +
  xlim(0,28) 
plot(MPD_Pfeiler_abundance_removal_fig)

###MNTD######
MNTD_Pfeiler_abundance_removal <- ses.mntd(Pfeiler_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                                       abundance.weighted=FALSE, runs = 5000, iterations = 5000) #all Pfeiler MNTD is SES is -1.04

MNTD_Pfeiler_abundance_removal <- MNTD_Pfeiler_abundance_removal[-c(27,28),]
MNTD_Pfeiler_abundance_removal$Abundance_rank <- c(1:26)

MNTD_Pfeiler_abundance_removal_fig <- ggplot(data= MNTD_Pfeiler_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-1.05, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mntd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -1.05, col = "lightgrey") +
  xlim(0,28) 
plot(MNTD_Pfeiler_abundance_removal_fig)

#Road--------------
Road_abundance_matrix <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/Road_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Road_abundance_matrix[34,Road_abundance_matrix[34,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Road_abundance_removal <- ses.pd(Road_abundance_matrix, pruned.tree, null.model = c("sample.pool"),
                                       runs = 5000, include.root=TRUE) #all Road PD is SES is -0.01

PD_Road_abundance_removal$Abundance_rank <- c(1:34)
PD_Road_abundance_removal <- PD_Road_abundance_removal[-c(33,34),]

PD_Road_abundance_removal_fig <- ggplot(data= PD_Road_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  xlim(0,32) 
plot(PD_Road_abundance_removal_fig)

###MPD########
MPD_Road_abundance_removal <- ses.mpd(Road_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), 
                                         abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all Road MPD is SES is -0.99

MPD_Road_abundance_removal <- MPD_Road_abundance_removal[-c(33,34),]
MPD_Road_abundance_removal$Abundance_rank <- c(1:32)

MPD_Road_abundance_removal_fig <- ggplot(data= MPD_Road_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-1, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mpd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -1, col = "lightgrey") +
  xlim(0,32) 
plot(MPD_Road_abundance_removal_fig)

###MNTD######
MNTD_Road_abundance_removal <- ses.mntd(Road_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                                           abundance.weighted=FALSE, runs = 5000, iterations = 5000) #all Road MNTD is SES is 0.015

MNTD_Road_abundance_removal <- MNTD_Road_abundance_removal[-c(33,34),]
MNTD_Road_abundance_removal$Abundance_rank <- c(1:32)

MNTD_Road_abundance_removal_fig <- ggplot(data= MNTD_Road_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.015, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mntd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.015, col = "lightgrey") +
  xlim(0,32) 
plot(MNTD_Road_abundance_removal_fig)

#Durbin-Watson test for autocorrelation----------
model <- lm(pd.obs.z ~ Abundance_rank, data = PD_Road_abundance_removal)
durbinWatsonTest(model, max.lag = 4)

#combine figures using patchwork
fig3_pbm <- (PD_PBM_abundance_removal_fig | MPD_PBM_abundance_removal_fig | MNTD_PBM_abundance_removal_fig) +
  plot_layout(axis_titles = "collect")+
  plot_annotation(title = "High elevation (3380 m)")
fig3_pbm

fig3_pfeiler <- (PD_Pfeiler_abundance_removal_fig | MPD_Pfeiler_abundance_removal_fig | MNTD_Pfeiler_abundance_removal_fig)+
  plot_layout(axis_titles = "collect")+
  plot_annotation(title = "Middle elevation (3165 m)")
fig3_pfeiler

fig3_road <- (PD_Road_abundance_removal_fig | MPD_Road_abundance_removal_fig | MNTD_Road_abundance_removal_fig) +
  plot_layout(axis_titles = "collect")+
  plot_annotation(title = "Low elevation (2815 m)")

fig3_all <- (fig3_pbm / fig3_pfeiler / fig3_road) +
  plot_layout(axis_titles = "collect")+
  plot_annotation(tag_levels = 'A')
fig3_all


#all figures together in ggplot-------------
##clean up data sets to have matching column names and columns for site and metric type#####

###PBM#####
###PD####
PD_PBM_abundance_removal = subset(PD_PBM_abundance_removal, select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs) ) #remove unnecessary columns

PD_PBM_abundance_removal <- PD_PBM_abundance_removal %>% 
  rename(SES = pd.obs.z,
    P_value = pd.obs.p) #rename columns to match other datasets 

PD_PBM_abundance_removal$Type <- c("PD") #add column for metric type 
PD_PBM_abundance_removal$Site <- c("High elevation (3380 m)") #add column for site name 

###MPD####
MPD_PBM_abundance_removal = subset(MPD_PBM_abundance_removal, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns

MPD_PBM_abundance_removal <- MPD_PBM_abundance_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

MPD_PBM_abundance_removal$Type <- c("MPD") #add column for metric type 
MPD_PBM_abundance_removal$Site <- c("High elevation (3380 m)") #add column for site name 

###MNTD####
MNTD_PBM_abundance_removal = subset(MNTD_PBM_abundance_removal, select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns

MNTD_PBM_abundance_removal <- MNTD_PBM_abundance_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

MNTD_PBM_abundance_removal$Type <- c("MNTD") #add column for metric type 
MNTD_PBM_abundance_removal$Site <- c("High elevation (3380 m)") #add column for site name 

###Pfeiler###
###PD####
PD_Pfeiler_abundance_removal = subset(PD_Pfeiler_abundance_removal, select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs) ) #remove unnecessary columns
PD_Pfeiler_abundance_removal <- PD_Pfeiler_abundance_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

PD_Pfeiler_abundance_removal$Type <- c("PD") #add column for metric type 
PD_Pfeiler_abundance_removal$Site <- c("Middle elevation (3165 m)") #add column for site name

###MPD####
MPD_Pfeiler_abundance_removal = subset(MPD_Pfeiler_abundance_removal, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns

MPD_Pfeiler_abundance_removal <- MPD_Pfeiler_abundance_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

MPD_Pfeiler_abundance_removal$Type <- c("MPD") #add column for metric type 
MPD_Pfeiler_abundance_removal$Site <- c("Middle elevation (3165 m)") #add column for site name

###MNTD#####
MNTD_Pfeiler_abundance_removal = subset(MNTD_Pfeiler_abundance_removal, select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns

MNTD_Pfeiler_abundance_removal <- MNTD_Pfeiler_abundance_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

MNTD_Pfeiler_abundance_removal$Type <- c("MNTD") #add column for metric type 
MNTD_Pfeiler_abundance_removal$Site <- c("Middle elevation (3165 m)") #add column for site name

###Road####
###PD####
PD_Road_abundance_removal = subset(PD_Road_abundance_removal, select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs) ) #remove unnecessary columns
PD_Road_abundance_removal <- PD_Road_abundance_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

PD_Road_abundance_removal$Type <- c("PD") #add column for metric type 
PD_Road_abundance_removal$Site <- c("Low elevation (2815 m)") #add column for site name 

###MPD###
MPD_Road_abundance_removal = subset(MPD_Road_abundance_removal, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns
MPD_Road_abundance_removal <- MPD_Road_abundance_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

MPD_Road_abundance_removal$Type <- c("MPD") #add column for metric type 
MPD_Road_abundance_removal$Site <- c("Low elevation (2815 m)") #add column for site name 

###MNTD#####
MNTD_Road_abundance_removal = subset(MNTD_Road_abundance_removal, select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns
MNTD_Road_abundance_removal <- MNTD_Road_abundance_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

MNTD_Road_abundance_removal$Type <- c("MNTD") #add column for metric type 
MNTD_Road_abundance_removal$Site <- c("Low elevation (2815 m)") #add column for site name 

##combine into one big dataset####
pbm <- rbind(PD_PBM_abundance_removal, MPD_PBM_abundance_removal, MNTD_PBM_abundance_removal)

pfeiler <- rbind(PD_Pfeiler_abundance_removal, MPD_Pfeiler_abundance_removal, MNTD_Pfeiler_abundance_removal)

road <- rbind(PD_Road_abundance_removal, MPD_Road_abundance_removal, MNTD_Road_abundance_removal)

all_sites_abundance_df <- rbind(pbm, pfeiler, road)

##make one figure with facet wrapping####
dummy <- data.frame(Site = c("High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)", "High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)", "High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)"),Type = c("PD","PD","PD","MPD","MPD","MPD","MNTD","MNTD","MNTD"), Z = c(0.57, -0.6, -0.01, 0.9, 0.02, -0.99, 0.12, -1.04, 0.015, 0.57, -0.6, -0.01, 0.9, 0.02, -0.99, 0.12, -1.04, 0.015, 0.57, -0.6, -0.01, 0.9, 0.02, -0.99, 0.12, -1.04, 0.015))

dummy <- dummy[-c(10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27), ] 

test <- left_join(all_sites_abundance_df, dummy, by = c("Site","Type"))
test$Site <- factor(test$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

ggplot(data= test) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, yend=SES, y=Z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=SES), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES") +
  scale_y_continuous(name="SES", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_abline(data = dummy, aes(intercept = Z, slope = 0)) +
  xlim(0,32) +
  facet_wrap(Type~Site)

