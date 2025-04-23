#code for removing single species ranked by abundance, Fig 4

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
PBM_abundance_matrix <- read.table("comm_phylo_analyses/Removing_species/comm_matrices/PBM_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_abundance_matrix[37,PBM_abundance_matrix[37,]>0]), warnings = F)$phy

#check tree
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)
specieslist <- pruned.tree$tip.label
specieslist <- as.data.frame(specieslist)

###PD#####
PD_PBM_abundance_removal <- ses.pd(PBM_abundance_matrix, pruned.tree, null.model = c("sample.pool"),runs = 5000, include.root=TRUE) #all PBM PD is SES is 0.612

PD_PBM_abundance_removal <- PD_PBM_abundance_removal[-c(37,38),] #remove 'all' rows
PD_PBM_abundance_removal$Abundance_rank <- c(1:36) #rank species by abundance in new column

#figure
PD_PBM_abundance_removal_fig <- ggplot(data= PD_PBM_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.61, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.61, col = "lightgrey") +
  xlim(0,36) 

plot(PD_PBM_abundance_removal_fig)

###MPD########
MPD_PBM_abundance_removal <- ses.mpd(PBM_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all PBM MPD is SES is 0.84

MPD_PBM_abundance_removal <- MPD_PBM_abundance_removal[-c(37,38),] #remove 'all' rows
MPD_PBM_abundance_removal$Abundance_rank <- c(1:36) #rank species by abundance in new column

#figure
MPD_PBM_abundance_removal_fig <- ggplot(data= MPD_PBM_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.84, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mpd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5)) +
  theme_classic(14) +
  geom_hline(yintercept = 0.84, col = "lightgrey") +
  xlim(0,36) 
plot(MPD_PBM_abundance_removal_fig)

###MNTD######
MNTD_PBM_abundance_removal <- ses.mntd(PBM_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted=FALSE, runs = 5000, iterations = 5000) #all PBM MNTD is SES is 0.04

MNTD_PBM_abundance_removal <- MNTD_PBM_abundance_removal[-c(37,38),] #remove 'all' rows
MNTD_PBM_abundance_removal$Abundance_rank <- c(1:36) #rank species by abundance in new column

#figure
MNTD_PBM_abundance_removal_fig <- ggplot(data= MNTD_PBM_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.04, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mntd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.04, col = "lightgrey") +
  xlim(0,36) 
plot(MNTD_PBM_abundance_removal_fig)

#Pfeiler--------
##make community data matrix#### 
Pfeiler_abundance_matrix <- read.table("comm_phylo_analyses/Removing_species/comm_matrices/Pfeiler_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Pfeiler_abundance_matrix[31,Pfeiler_abundance_matrix[31,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Pfeiler_abundance_removal <- ses.pd(Pfeiler_abundance_matrix, pruned.tree, null.model = c("sample.pool"),runs = 5000, include.root=TRUE) #all Pfeiler PD is SES is -1.3

PD_Pfeiler_abundance_removal <- PD_Pfeiler_abundance_removal[-c(30,31),] #remove 'all' rows
PD_Pfeiler_abundance_removal$Abundance_rank <- c(1:29) #rank species by abundance in new column

#figure
PD_Pfeiler_abundance_removal_fig <- ggplot(data= PD_Pfeiler_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-1.3, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.8, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -1.3, col = "lightgrey") +
  xlim(0,29) 
plot(PD_Pfeiler_abundance_removal_fig)

###MPD########
MPD_Pfeiler_abundance_removal <- ses.mpd(Pfeiler_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all Pfeiler MPD is SES is -0.5

MPD_Pfeiler_abundance_removal <- MPD_Pfeiler_abundance_removal[-c(30,31),] #remove 'all' rows
MPD_Pfeiler_abundance_removal$Abundance_rank <- c(1:29) #rank species by abundance in new column

#figure
MPD_Pfeiler_abundance_removal_fig <- ggplot(data= MPD_Pfeiler_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.5, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mpd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  scale_y_continuous(name="SES MPD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.5, col = "lightgrey") +
  xlim(0,29) 
plot(MPD_Pfeiler_abundance_removal_fig)

###MNTD######
MNTD_Pfeiler_abundance_removal <- ses.mntd(Pfeiler_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted=FALSE, runs = 5000, iterations = 5000) #all Pfeiler MNTD is SES is -1.3

MNTD_Pfeiler_abundance_removal <- MNTD_Pfeiler_abundance_removal[-c(30,31),] #remove 'all' rows
MNTD_Pfeiler_abundance_removal$Abundance_rank <- c(1:29) #rank species by abundance in new column

#figure
MNTD_Pfeiler_abundance_removal_fig <- ggplot(data= MNTD_Pfeiler_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-1.3, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mntd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.9, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -1.3, col = "lightgrey") +
  xlim(0,29) 
plot(MNTD_Pfeiler_abundance_removal_fig)

#Road--------------
Road_abundance_matrix <- read.table("comm_phylo_analyses/Removing_species/comm_matrices/Road_abundance_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Road_abundance_matrix[41,Road_abundance_matrix[41,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Road_abundance_removal <- ses.pd(Road_abundance_matrix, pruned.tree, null.model = c("sample.pool"),runs = 5000, include.root=TRUE) #all Road PD is SES is -0.05

PD_Road_abundance_removal <- PD_Road_abundance_removal[-c(40,41),] #remove 'all' rows
PD_Road_abundance_removal$Abundance_rank <- c(1:39) #rank species by abundance in new column

#figure
PD_Road_abundance_removal_fig <- ggplot(data= PD_Road_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.05, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.05, col = "lightgrey") +
  xlim(0,39) 
plot(PD_Road_abundance_removal_fig)

###MPD########
MPD_Road_abundance_removal <- ses.mpd(Road_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted = FALSE, runs = 5000, iterations = 5000) #all Road MPD is SES is -0.86

MPD_Road_abundance_removal <- MPD_Road_abundance_removal[-c(40,41),] #remove 'all' rows
MPD_Road_abundance_removal$Abundance_rank <- c(1:39) #rank species by abundance in new column

#figure
MPD_Road_abundance_removal_fig <- ggplot(data= MPD_Road_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.86, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mpd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.86, col = "lightgrey") +
  xlim(0,39) 
plot(MPD_Road_abundance_removal_fig)

###MNTD######
MNTD_Road_abundance_removal <- ses.mntd(Road_abundance_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted=FALSE, runs = 5000, iterations = 5000) #all Road MNTD is SES is 0.33

MNTD_Road_abundance_removal <- MNTD_Road_abundance_removal[-c(40,41),] #remove 'all' rows
MNTD_Road_abundance_removal$Abundance_rank <- c(1:39) #rank species by abundance in new column

#figure
MNTD_Road_abundance_removal_fig <- ggplot(data= MNTD_Road_abundance_removal) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.33, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mntd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MNTD") +
  scale_y_continuous(name="SES MNTD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.33, col = "lightgrey") +
  xlim(0,39) 

plot(MNTD_Road_abundance_removal_fig)

#Durbin-Watson test for autocorrelation----------
model <- lm(pd.obs.z ~ Abundance_rank, data = PD_Pfeiler_abundance_removal)
durbinWatsonTest(model, max.lag = 3)

#combine figures using patchwork, not used in final figures
fig3_pbm <- (PD_PBM_abundance_removal_fig | PD_PBM_abundance_removal_fig | MNTD_PBM_abundance_removal_fig) +
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
pbm <- rbind(PD_PBM_abundance_removal, MPD_PBM_abundance_removal, 
             MNTD_PBM_abundance_removal)

pfeiler <- rbind(PD_Pfeiler_abundance_removal, MPD_Pfeiler_abundance_removal, MNTD_Pfeiler_abundance_removal)

road <- rbind(PD_Road_abundance_removal, MPD_Road_abundance_removal, 
              MNTD_Road_abundance_removal)

all_sites_abundance_df <- rbind(pbm, pfeiler, road)

##make one figure with facet wrapping####

#make dataframe with sites, phylogenetic metrics and the whole community value
#of that metric to use as baseline in figure

dummy <- data.frame(Site = c("High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)", "High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)", "High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)"),Type = c("PD","PD","PD","MPD","MPD","MPD","MNTD","MNTD","MNTD"), Z = c(0.61, -1.3, -0.05, 0.84, -0.5, -0.86, 0.04, -1.3, 0.33, 0.61, -1.3, -0.05, 0.84, -0.5, -0.86, 0.04, -1.3, 0.33, 0.61, -1.3, -0.05, 0.84, -0.5, -0.86, 0.04, -1.3, 0.33))

#remove extra rows
dummy <- dummy[-c(10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27), ] 

#join baseline 'dummy' data with species removal data
test <- left_join(all_sites_abundance_df, dummy, by = c("Site","Type"))
test$Site = factor(test$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

#reorder sites by elevation
test$Site_f = factor(test$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

#figure
ggplot(data= test) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, yend=SES, y=Z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=SES), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("Standard effect size") +
  scale_y_continuous(name="Standard effect size", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),
                     limits=c(-1.9, 1.5))+
  theme_bw(base_size = 22) +
  xlim(0,39) +
  geom_abline(data = test, aes(intercept = Z, slope = 0)) +
  facet_grid(Type~Site_f)

ggsave("figures/AJB_final/Fig4.pdf", height = 14.5, width = 14.5)
