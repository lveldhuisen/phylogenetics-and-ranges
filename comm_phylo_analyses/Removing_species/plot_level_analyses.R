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
library(forstringr)
library(reshape2)

#baseline pd for all plots
all_matrix <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/all_plots_baseline.csv")

#remove extra rows at the bottom
all_matrix <- na.omit(all_matrix)

#switch first column to row names
all_matrix <- all_matrix %>% remove_rownames %>% column_to_rownames(var='X')

#calculate baseline pd for all plots-----
all_plots_pd <- ses.pd(all_matrix, pruned.tree, null.model = c("sample.pool"),
                       runs = 5000, include.root=TRUE) 
#abundance------
#PBM-------

##plot 1########
PBM1 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pbm1_a.csv")
#switch first column to row names
PBM1 <- PBM1 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pbm1_pd <- ses.pd(PBM1, pruned.tree, null.model = c("sample.pool"),
                       runs = 5000, include.root=TRUE) 

pbm1_pd$Abundance_rank <- c(1:17) #rank species by abundance in new column

#figure
pbm1_fig <- ggplot(data= pbm1_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.06, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = 0.06, col = "lightgrey") +
  xlim(0,17) 

plot(pbm1_fig)

##plot 2#####
PBM2 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pbm2_a.csv")
#switch first column to row names
PBM2 <- PBM2 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pbm2_pd <- ses.pd(PBM2, pruned.tree, null.model = c("sample.pool"),
                  runs = 5000, include.root=TRUE) 

pbm2_pd$Abundance_rank <- c(1:7) #rank species by abundance in new column

#figure
pbm2_fig <- ggplot(data= pbm2_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.54, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-2, 1.7))+
  theme_classic(14) +
  geom_hline(yintercept = -0.54, col = "lightgrey") +
  xlim(0,7) 

plot(pbm2_fig)

##plot 3####
PBM3 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pbm3_a.csv")
#switch first column to row names
PBM3 <- PBM3 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pbm3_pd <- ses.pd(PBM3, pruned.tree, null.model = c("sample.pool"),
                  runs = 5000, include.root=TRUE) 

pbm3_pd$Abundance_rank <- c(1:16) #rank species by abundance in new column

#figure
pbm3_fig <- ggplot(data= pbm3_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=1.58, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1, 2))+
  theme_classic(14) +
  geom_hline(yintercept = 1.58, col = "lightgrey") +
  xlim(0,16) 

plot(pbm3_fig)

##plot 4 #####
PBM4 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pbm4_a.csv")
#switch first column to row names
PBM4 <- PBM4 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pbm4_pd <- ses.pd(PBM4, pruned.tree, null.model = c("sample.pool"),
                  runs = 5000, include.root=TRUE) 

pbm4_pd$Abundance_rank <- c(1:15) #rank species by abundance in new column

#figure
pbm4_fig <- ggplot(data= pbm4_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.23, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.5))+
  theme_classic(14) +
  geom_hline(yintercept = -0.23, col = "lightgrey") +
  xlim(0,15) 

plot(pbm4_fig)

##plot 5####

PBM5 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pbm5_a.csv")
#switch first column to row names
PBM5 <- PBM5 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pbm5_pd <- ses.pd(PBM5, pruned.tree, null.model = c("sample.pool"),
                  runs = 5000, include.root=TRUE) 

pbm5_pd$Abundance_rank <- c(1:10) #rank species by abundance in new column

#figure
pbm5_fig <- ggplot(data= pbm5_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.46, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = 0.46, col = "lightgrey") +
  xlim(0,11) 

plot(pbm5_fig)

##combine all PBM figures into one wtih patchwork###
PBM_PD <- pbm1_fig + pbm2_fig + pbm3_fig + pbm4_fig + pbm5_fig
plot(PBM_PD)

#Pfeiler--------

##plot 1###
pfeiler1 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pfeiler1_a.csv")
#switch first column to row names
pfeiler1 <- pfeiler1 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pfeiler1_pd <- ses.pd(pfeiler1, pruned.tree, null.model = c("sample.pool"),
                  runs = 5000, include.root=TRUE) 

pfeiler1_pd$Abundance_rank <- c(1:13) #rank species by abundance in new column

#figure
pfeiler1_fig <- ggplot(data= pfeiler1_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.26, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = -0.26, col = "lightgrey") +
  xlim(0,13) 

plot(pfeiler1_fig)

##plot 2####
pfeiler2 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pfeiler2_a.csv")
#switch first column to row names
pfeiler2 <- pfeiler2 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pfeiler2_pd <- ses.pd(pfeiler2, pruned.tree, null.model = c("sample.pool"),
                      runs = 5000, include.root=TRUE) 

pfeiler2_pd$Abundance_rank <- c(1:19) #rank species by abundance in new column

#figure
pfeiler2_fig <- ggplot(data= pfeiler2_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.2, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-2, 1.3))+
  theme_classic(14) +
  geom_hline(yintercept = -0.2, col = "lightgrey") +
  xlim(0,19) 

plot(pfeiler2_fig)

##plot 3###
pfeiler3 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pfeiler3_a.csv")
#switch first column to row names
pfeiler3 <- pfeiler3 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pfeiler3_pd <- ses.pd(pfeiler3, pruned.tree, null.model = c("sample.pool"),
                      runs = 5000, include.root=TRUE) 

pfeiler3_pd$Abundance_rank <- c(1:19) #rank species by abundance in new column

#figure
pfeiler3_fig <- ggplot(data= pfeiler3_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.57, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = -0.57, col = "lightgrey") +
  xlim(0,19) 

plot(pfeiler3_fig)

##plot 4###
pfeiler4 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pfeiler4_a.csv")
#switch first column to row names
pfeiler4 <- pfeiler4 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pfeiler4_pd <- ses.pd(pfeiler4, pruned.tree, null.model = c("sample.pool"),
                      runs = 5000, include.root=TRUE) 

pfeiler4_pd$Abundance_rank <- c(1:19) #rank species by abundance in new column

#figure
pfeiler4_fig <- ggplot(data= pfeiler4_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=-0.11, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = -0.11, col = "lightgrey") +
  xlim(0,19) 

plot(pfeiler4_fig)

##plot 5####
pfeiler5 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/pfeiler5_a.csv")
#switch first column to row names
pfeiler5 <- pfeiler5 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
pfeiler5_pd <- ses.pd(pfeiler5, pruned.tree, null.model = c("sample.pool"),
                      runs = 5000, include.root=TRUE) 

pfeiler5_pd$Abundance_rank <- c(1:21) #rank species by abundance in new column

#figure
pfeiler5_fig <- ggplot(data= pfeiler5_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.28, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = 0.28, col = "lightgrey") +
  xlim(0,21) 

plot(pfeiler5_fig)

#Road--------

##plot 1####
road1 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/road1_a.csv")
#switch first column to row names
road1 <- road1 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
road1_pd <- ses.pd(road1, pruned.tree, null.model = c("sample.pool"),
                      runs = 5000, include.root=TRUE) 

road1_pd$Abundance_rank <- c(1:20) #rank species by abundance in new column

#figure
road1_fig <- ggplot(data= road1_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=1.49, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.5, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = 1.49, col = "lightgrey") +
  xlim(0,20) 

plot(road1_fig)

##plot 2####
road2 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/road2_a.csv")
#switch first column to row names
road2 <- road2 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
road2_pd <- ses.pd(road2, pruned.tree, null.model = c("sample.pool"),
                   runs = 5000, include.root=TRUE) 

road2_pd$Abundance_rank <- c(1:19) #rank species by abundance in new column

#figure
road2_fig <- ggplot(data= road2_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.73, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.7, 1.3))+
  theme_classic(14) +
  geom_hline(yintercept = 0.73, col = "lightgrey") +
  xlim(0,19) 

plot(road2_fig)

##plot 3####

road3 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/road3_a.csv")
#switch first column to row names
road3 <- road3 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
road3_pd <- ses.pd(road3, pruned.tree, null.model = c("sample.pool"),
                   runs = 5000, include.root=TRUE) 

road3_pd$Abundance_rank <- c(1:17) #rank species by abundance in new column

#figure
road3_fig <- ggplot(data= road3_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.84, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.5, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = 0.84, col = "lightgrey") +
  xlim(0,17) 

plot(road3_fig)

##plot 4####
road4 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/road4_a.csv")
#switch first column to row names
road4 <- road4 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
road4_pd <- ses.pd(road4, pruned.tree, null.model = c("sample.pool"),
                   runs = 5000, include.root=TRUE) 

road4_pd$Abundance_rank <- c(1:16) #rank species by abundance in new column

#figure
road4_fig <- ggplot(data= road4_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=0.05, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.5, 1.8))+
  theme_classic(14) +
  geom_hline(yintercept = 0.05, col = "lightgrey") +
  xlim(0,16) 

plot(road4_fig)

##plot 5####

road5 <- read.csv("comm_phylo_analyses/Removing_species/plot_level_matrices/road5_a.csv")
#switch first column to row names
road5 <- road5 %>% remove_rownames %>% column_to_rownames(var='X')

#calculate pd
road5_pd <- ses.pd(road5, pruned.tree, null.model = c("sample.pool"),
                   runs = 5000, include.root=TRUE) 

road5_pd$Abundance_rank <- c(1:12) #rank species by abundance in new column

#figure
road5_fig <- ggplot(data= road5_pd) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=1.44, yend=pd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=pd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES PD") +
  scale_y_continuous(name="SES PD", breaks = c(-1.5, -1, -0.5, 0, 1, 0.5, 1, 1.5),limits=c(-1.5, 2.1))+
  theme_classic(14) +
  geom_hline(yintercept = 1.44, col = "lightgrey") +
  xlim(0,12) 

plot(road5_fig)
