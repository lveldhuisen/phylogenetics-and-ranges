install.packages("fabricatr")
library(fabricatr)
library(ggplot2)
library(gridExtra)
library(rgdal)
library(sf)
library(terra)
library(sp)
library(rgeos)
library(dplyr)
library(picante)
library(geiger)
library(ape)
library(vegan)
library(forcats)




#group species by range size-------------------------------------------------
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")
all_df <- read.csv("results_all.csv")

#import S&B phylogeny---------------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")

##import S&B18 tree and check data## 
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

##make community data matrix 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
matrix <- read.table("community_matrix_rangesize.txt", sep = "\t", header = T, row.names = 1)

#Faith's PD-------------------------------------------------------------------
pd(matrix, pruned.tree, include.root = T)


###SES###############################
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")

pruned.tree <- treedata(SBtree, unlist(matrix[16,matrix[16,]>0]), warnings = F)$phy
plot(pruned.tree)

#standard effect size using picante 
library(picante)
ses.pd(matrix, pruned.tree, null.model = c("sample.pool"),
       runs = 5000, include.root=TRUE)

#make distance matrix for MPD and MNTD
dist.mat <- cophenetic(pruned.tree)

ses.mpd(matrix, dist.mat, null.model = c("sample.pool"),
        abundance.weighted = FALSE, runs = 5000, iterations = 5000)

ses.mntd(matrix, dist.mat, null.model = c("sample.pool"),
         abundance.weighted=FALSE, runs = 5000, iterations = 5000)

#figures###
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/comm_phylo_analyses")
phylo_df <- read.csv("phylo_metrics_rangesize.csv")

##subset data by site-----
subset_road <- subset(phylo_df, 
                      Site %in% c("Road"))
subset_road <- subset_road[-c(7,8,9), ]

subset_pfeiler <- subset(phylo_df, 
                         Site %in% c("Pfeiler"))
subset_pfeiler <- subset_pfeiler[-c(10,11,12), ]

subset_PBM <- subset(phylo_df, 
                     Site %in% c("PBM"))
subset_PBM <- subset_PBM[-c(10,11,12), ]

##figures#####
PBM_fig <- ggplot(subset_PBM, aes(fill=Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("PBM - high elevation")
plot(PBM_fig)

pfeiler_fig <- ggplot(subset_pfeiler, aes(fill=Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("Pfeiler - middle elevation")
plot(pfeiler_fig)

road_fig <- ggplot(subset_road, aes(fill=Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("Road - low elevation")
plot(road_fig)

all_fig <- ggplot(phylo_df, aes(fill=Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("All sites together")
plot(all_fig)

#fig without range size categories
general_fig <- ggplot(phylo_df, aes(fill=Type, y=SES, x=fct_relevel(Site, c("all")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-3,3) 
plot(general_fig)  
ß