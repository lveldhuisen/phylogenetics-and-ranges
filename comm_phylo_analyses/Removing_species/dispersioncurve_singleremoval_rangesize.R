#code to generate curves of phylogenetic dispersion while removing 1 species at a time, ranked by decreasing range size 
#will result in an individual impact on PD/MPD/MNTD for each species 

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

#import S&B phylogeny#
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#PBM----------
##make community data matrix#### 
PBM_range_matrix <- read.table("comm_phylo_analyses/Removing_species/comm_matrices/PBMrangesize_comm_matrix_removal.txt", 
                               sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(PBM_range_matrix[38,PBM_range_matrix[38,]>0]), 
                        warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

co_matrix <- cophenetic(pruned.tree)

###PD#####
PD_PBM_range_removal <- ses.pd(PBM_range_matrix, pruned.tree, 
                               null.model = c("sample.pool"),
                    runs = 5000, include.root=TRUE) #all PBM PD is SES is 0.61

PD_PBM_range_removal <- PD_PBM_range_removal[-c(37,38),]
PD_PBM_range_removal$Range_size_rank <- c(1:36)

PD_PBM_rangesize_removal_fig <- ggplot(data= PD_PBM_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0.61, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES PD") +
  ylim(-0.5,2) +
  theme_classic(14) +
  geom_hline(yintercept = 0.61, col = "lightgrey") +
  xlim(0,36) 
plot(PD_PBM_rangesize_removal_fig)


###MPD#####
MPD_PBM_range_removal <- ses.mpd(PBM_range_matrix, co_matrix, null.model = c("sample.pool"), 
        abundance.weighted = FALSE, runs = 5000, iterations = 5000)
#all PBM MPD is SES is 0.93

MPD_PBM_range_removal <- MPD_PBM_range_removal[-c(37,38),]
MPD_PBM_range_removal$Range_size_rank <- c(1:36)

MPD_PBM_range_removal_fig <- ggplot(data= MPD_PBM_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0.93,yend= mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mpd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES MPD") +
  ylim(0,2) +
  theme_classic(14) +
  geom_hline(yintercept = 0.93, col = "lightgrey") +
  xlim(0,36)
plot(MPD_PBM_range_removal_fig)

###MNTD####
MNTD_PBM_range_removal <- ses.mntd(PBM_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
                        abundance.weighted=FALSE, runs = 5000, iterations = 5000)
#all PBM MNTD is SES is -0.03
MNTD_PBM_range_removal <- MNTD_PBM_range_removal[-c(37,38),]
MNTD_PBM_range_removal$Range_size_rank <- c(1:36)

MNTD_PBM_range_removal_fig <- ggplot(data= MNTD_PBM_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-0.03, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mntd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES MNTD") +
  ylim(-1,1) +
  theme_classic(14) +
  geom_hline(yintercept = -0.03, col = "lightgrey") +
  xlim(0,36) 
plot(MNTD_PBM_range_removal_fig)

#Pfeiler--------
Pfeiler_range_matrix <- read.table("comm_phylo_analyses/Removing_species/comm_matrices/Pfeiler_ranges_comm_matrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Pfeiler_range_matrix[31,Pfeiler_range_matrix[31,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Pfeiler_range_removal <- ses.pd(Pfeiler_range_matrix, pruned.tree, null.model = c("sample.pool"), runs = 5000, include.root=TRUE) #all Pfeiler  PD is SES is -1.3

PD_Pfeiler_range_removal <- PD_Pfeiler_range_removal[-c(30,31),]
PD_Pfeiler_range_removal$Range_size_rank <- c(1:29)


PD_Pfeiler_rangesize_removal_fig <- ggplot(data= PD_Pfeiler_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-1.3, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES PD") +
  ylim(-2,0.5) +
  theme_classic(14) +
  geom_hline(yintercept = -1.3, col = "lightgrey") +
  xlim(0,29) 

plot(PD_Pfeiler_rangesize_removal_fig)

###MPD#####
MPD_Pfeiler_range_removal <- ses.mpd(Pfeiler_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted = FALSE, runs = 5000, iterations = 5000)

#all Pfeiler MPD is SES is -0.6

MPD_Pfeiler_range_removal <- MPD_Pfeiler_range_removal[-c(30,31),]
MPD_Pfeiler_range_removal$Range_size_rank <- c(1:29)

MPD_Pfeiler_range_removal_fig <- ggplot(data= MPD_Pfeiler_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-0.6, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mpd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES MPD") +
  ylim(-1.2,1) +
  theme_classic(14) +
  geom_hline(yintercept = -0.6, col = "lightgrey") +
  xlim(0,29) 
plot(MPD_Pfeiler_range_removal_fig)

###MNTD####
MNTD_Pfeiler_range_removal <- ses.mntd(Pfeiler_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted=FALSE, runs = 5000, iterations = 5000)

#all Pfeiler MNTD is SES is -1.41

MNTD_Pfeiler_range_removal <- MNTD_Pfeiler_range_removal[-c(30,31),]
MNTD_Pfeiler_range_removal$Range_size_rank <- c(1:29)

MNTD_Pfeiler_range_removal_fig <- ggplot(data= MNTD_Pfeiler_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-1.41, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mntd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES MNTD") +
  ylim(-2,1) +
  theme_classic(14) +
  geom_hline(yintercept = -1.41, col = "lightgrey") +
  xlim(0,29) 
plot(MNTD_Pfeiler_range_removal_fig)

#Road-------
Road_range_matrix <- read.table("comm_phylo_analyses/Removing_species/comm_matrices/Road_range_commmatrix_removal.txt", sep = "\t", header = T, row.names = 1)

##prune tree#####
pruned.tree <- treedata(SBtree, unlist(Road_range_matrix[41,Road_range_matrix[41,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

###PD#####
PD_Road_range_removal <- ses.pd(Road_range_matrix, pruned.tree, null.model = c("sample.pool"),
                                   runs = 5000, include.root=TRUE)
#all Road  PD is SES is 0

PD_Road_range_removal <- PD_Road_range_removal[-c(40,41),]
PD_Road_range_removal$Range_size_rank <- c(1:39)


PD_Road_rangesize_removal_fig <- ggplot(data= PD_Road_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=pd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES PD") +
  ylim(-1,0.5) +
  theme_classic(14) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  xlim(0,39) 

plot(PD_Road_rangesize_removal_fig)

###MPD####
MPD_Road_range_removal <- ses.mpd(Road_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted = FALSE, runs = 5000, iterations = 5000)
#all Road MPD is SES is -0.76

MPD_Road_range_removal <- MPD_Road_range_removal[-c(40,41),]
MPD_Road_range_removal$Range_size_rank <- c(1:39)

MPD_Road_range_removal_fig <- ggplot(data= MPD_Road_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=-0.76, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mpd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES MPD") +
  ylim(-1.5,0.5) +
  theme_classic(14) +
  geom_hline(yintercept = -0.76, col = "lightgrey") +
  xlim(0,39) 
plot(MPD_Road_range_removal_fig)

###MNTD#####
MNTD_Road_range_removal <- ses.mntd(Road_range_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),)
 #all Road MNTD is SES is 0.42

MNTD_Road_range_removal <- MNTD_Road_range_removal[-c(40,41),]
MNTD_Road_range_removal$Range_size_rank <- c(1:39)

MNTD_Road_range_removal_fig <- ggplot(data= MNTD_Road_range_removal) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, y=0.42, yend=mntd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=mntd.obs.z), size = 2) +
  xlab("Range size rank of removed species (largest to smallest)") +
  ylab("SES MNTD") +
  ylim(-1,1) +
  theme_classic(14) +
  geom_hline(yintercept = 0.42, col = "lightgrey") +
  xlim(0,39) 
plot(MNTD_Road_range_removal_fig)


#Durbin-Watson test for autocorrelation----------
model <- lm(pd.obs.z ~ Range_size_rank, data = PD_Road_range_removal)
durbinWatsonTest(model, max.lag = 3)

#make combined figure using patchwork

fig2_PBM <- (PD_PBM_rangesize_removal_fig | MPD_PBM_range_removal_fig | MNTD_PBM_range_removal_fig) +
  plot_layout(axis_titles = 'collect')

fig2_PBM

fig2_pfeiler <- (PD_Pfeiler_rangesize_removal_fig | MPD_Pfeiler_range_removal_fig | MNTD_Pfeiler_range_removal_fig) +
  plot_layout(axis_titles = 'collect')+
  plot_annotation(title = "Middle elevation (3165 m)")
fig2_pfeiler

fig2_road <- (PD_Road_rangesize_removal_fig | MPD_Road_range_removal_fig | MNTD_Road_range_removal_fig) +
  plot_layout(axis_titles = 'collect')+
  plot_annotation(title = "Low elevation (2815 m)")
fig2_road

fig2_all <- fig2_PBM / fig2_pfeiler / fig2_road +
  plot_layout(axis_titles = "collect")+
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = 'A')
fig2_all

#all figures together in ggplot-------------
##clean up data sets to have matching column names and columns for site and metric type#####

##PBM#####
###PD####
PD_PBM_range_removal = subset(PD_PBM_range_removal, select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs) ) #remove unnecessary columns

PD_PBM_range_removal <- PD_PBM_range_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

PD_PBM_range_removal$Type <- c("PD") #add column for metric type 
PD_PBM_range_removal$Site <- c("High elevation (3380 m)") #add column for site name 

###MPD####
MPD_PBM_range_removal = subset(MPD_PBM_range_removal, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns

MPD_PBM_range_removal <- MPD_PBM_range_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

MPD_PBM_range_removal$Type <- c("MPD") #add column for metric type 
MPD_PBM_range_removal$Site <- c("High elevation (3380 m)") #add column for site name 

###MNTD####
MNTD_PBM_range_removal = subset(MNTD_PBM_range_removal, select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns

MNTD_PBM_range_removal <- MNTD_PBM_range_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

MNTD_PBM_range_removal$Type <- c("MNTD") #add column for metric type 
MNTD_PBM_range_removal$Site <- c("High elevation (3380 m)") #add column for site name 

##Pfeiler#####
###PD####
PD_Pfeiler_range_removal = subset(PD_Pfeiler_range_removal, select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs) ) #remove unnecessary columns
PD_Pfeiler_range_removal <- PD_Pfeiler_range_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

PD_Pfeiler_range_removal$Type <- c("PD") #add column for metric type 
PD_Pfeiler_range_removal$Site <- c("Middle elevation (3165 m)") #add column for site name

###MPD####
MPD_Pfeiler_range_removal = subset(MPD_Pfeiler_range_removal, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns

MPD_Pfeiler_range_removal <- MPD_Pfeiler_range_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

MPD_Pfeiler_range_removal$Type <- c("MPD") #add column for metric type 
MPD_Pfeiler_range_removal$Site <- c("Middle elevation (3165 m)") #add column for site name

###MNTD#####
MNTD_Pfeiler_range_removal = subset(MNTD_Pfeiler_range_removal, select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns

MNTD_Pfeiler_range_removal <- MNTD_Pfeiler_range_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

MNTD_Pfeiler_range_removal$Type <- c("MNTD") #add column for metric type 
MNTD_Pfeiler_range_removal$Site <- c("Middle elevation (3165 m)") #add column for site name

##Road####
###PD####
PD_Road_range_removal = subset(PD_Road_range_removal, select = -c(ntaxa,pd.obs,pd.rand.mean,pd.rand.sd,pd.obs.rank,runs) ) #remove unnecessary columns
PD_Road_range_removal <- PD_Road_range_removal %>% 
  rename(SES = pd.obs.z,
         P_value = pd.obs.p) #rename columns to match other datasets 

PD_Road_range_removal$Type <- c("PD") #add column for metric type 
PD_Road_range_removal$Site <- c("Low elevation (2815 m)") #add column for site name 

###MPD###
MPD_Road_range_removal = subset(MPD_Road_range_removal, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns
MPD_Road_range_removal <- MPD_Road_range_removal %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

MPD_Road_range_removal$Type <- c("MPD") #add column for metric type 
MPD_Road_range_removal$Site <- c("Low elevation (2815 m)") #add column for site name 

###MNTD#####
MNTD_Road_range_removal = subset(MNTD_Road_range_removal, select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns
MNTD_Road_range_removal <- MNTD_Road_range_removal %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

MNTD_Road_range_removal$Type <- c("MNTD") #add column for metric type 
MNTD_Road_range_removal$Site <- c("Low elevation (2815 m)") #add column for site name 

##combine into one big dataset####
pbm_range <- rbind(PD_PBM_range_removal, MPD_PBM_range_removal, MNTD_PBM_range_removal)

pfeiler_range <- rbind(PD_Pfeiler_range_removal, MPD_Pfeiler_range_removal, MNTD_Pfeiler_range_removal)

road_range <- rbind(PD_Road_range_removal, MPD_Road_range_removal, MNTD_Road_range_removal)

all_sites_range_df <- rbind(pbm_range, pfeiler_range, road_range)

##make one figure with facet wrapping####
dummy_range <- data.frame(Site = c("High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)", "High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)", "High elevation (3380 m)", "Middle elevation (3165 m)","Low elevation (2815 m)"),Type = c("PD","PD","PD","MPD","MPD","MPD","MNTD","MNTD","MNTD"), Z = c(0.61, -1.3, 0,
        0.93, -0.6, -0.76, -0.03, -1.41, 0.42,0.61, -1.3, 0,
        0.93, -0.6, -0.76, -0.03, -1.41, 0.42,0.61, -1.3, 0,
        0.93, -0.6, -0.76, -0.03, -1.41, 0.42))

dummy_range <- dummy_range[-c(10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27), ] 

test_range <- left_join(all_sites_range_df, dummy_range, by = c("Site","Type"))
test_range$Site = factor(test_range$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

test_range$Site_f = factor(test_range$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"
))

ggplot(data= test_range) + 
  geom_segment( aes(x=Range_size_rank, xend=Range_size_rank, yend=SES, y=Z), color="grey")+
  geom_point(mapping = aes(x=Range_size_rank, y=SES), size = 2) +
  xlab("Range size rank of removed species (biggest to smallest)") +
  ylab("Standard effect size") +
  scale_y_continuous(name="Standard effect size", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-2, 1.4))+
  theme_bw(14) +
  xlim(0,39) +
  geom_abline(data = test_range, aes(intercept = Z, slope = 0)) +
  facet_grid(Type~Site_f)

