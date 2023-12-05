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

#group species by suitability-------------------------------------------------
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")
suitability_df <- read.csv("specieslist_gbifcites_suitability.csv")
suitability_df <- suitability_df[-c(90, 91, 92, 93, 94, 95, 96, 97, 98, 99), ]

hist_SD<- ggplot(suitability_df, aes(x=suitability_df$Suitability_SD)) + 
  geom_histogram(binwidth = 25)
ggplot(suitability_df, aes(x=Suitability_median)) + 
  geom_histogram(binwidth = 100)
ggplot(suitability_df, aes(x=Suitability_mean)) + 
  geom_histogram(binwidth = 100)

model <- lm(suitability_SD ~ occurance.data..n., data = suitability_df)
summary(model)

ggplot(suitability_df, aes(x=occurance.data..n., y=suitability.SD)) + 
  geom_point() + stat_smooth(method = "lm")
summary(suitability_df$suitability.SD)

#import S&B phylogeny---------------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")

##import S&B18 tree and check data## 
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

##make community data matrix 
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")
matrix <- read.table("community_matrix.txt", sep = "\t", header = T, row.names = 1)
range_matrix <- read.table("community_matrix_rangesize.txt", sep = "\t", header = T, row.names = 1)

#Faith's PD-------------------------------------------------------------------
pd(matrix, pruned.tree, include.root = T)


###SES###############################
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 1")

pruned.tree <- treedata(SBtree, unlist(matrix[16,matrix[16,]>0]), warnings = F)$phy
plot(pruned.tree)

#standard effect size using picante 
install.packages("picante")
library(picante)
ses.pd(matrix, pruned.tree, null.model = c("sample.pool"),
       runs = 5000, include.root=TRUE)

dist.mat <- cophenetic(pruned.tree)

ses.mpd(matrix, dist.mat, null.model = c("sample.pool"),
        abundance.weighted = FALSE, runs = 5000, iterations = 5000)



ses.mpd(matrix, pruned.tree, null.model = c("sample.pool"),
        runs = 5000)

ses.mntd(matrix, dist.mat, null.model = c("sample.pool"),
         abundance.weighted=FALSE, runs = 5000, iterations = 5000)

#figures###
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")
phylo_df <- read.csv("phylo_results.csv")
phylo_df_test <- phylo_df[-c(1,2,3), ]
phylo_df <- phylo_df_test

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
PBM_fig <- ggplot(subset_PBM, aes(fill=Type, y=SES, x=fct_relevel(Suitability, c("Low","Medium","High")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Suitability category") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("PBM - high elevation")
plot(PBM_fig)

pfeiler_fig <- ggplot(subset_pfeiler, aes(fill=Type, y=SES, x=fct_relevel(Suitability, c("Low","Medium","High")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Suitability category") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("Pfeiler - middle elevation")
plot(pfeiler_fig)

road_fig <- ggplot(subset_road, aes(fill=Type, y=SES, x=fct_relevel(Suitability, c("Medium","High")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Suitability category") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("Road - low elevation")
plot(road_fig)

all_fig <- ggplot(phylo_df, aes(fill=Type, y=SES, x=fct_relevel(Suitability, c("Low","Medium","High")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Suitability category") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + ggtitle("All sites together")
plot(all_fig)

#phylogenetics for range size categories------------------------
#standard effect size using picante 
library(picante)

pruned.tree <- treedata(SBtree, unlist(range_matrix[16,range_matrix[16,]>0]), warnings = F)$phy

ses.pd(range_matrix, pruned.tree, null.model = c("sample.pool"),
       runs = 5000, include.root=TRUE)

dist.mat <- cophenetic(pruned.tree)

ses.mpd(range_matrix, dist.mat, null.model = c("sample.pool"),
        abundance.weighted = FALSE, runs = 5000, iterations = 5000)

ses.mpd(range_matrix, pruned.tree, null.model = c("sample.pool"),
        runs = 5000)

ses.mntd(range_matrix, dist.mat, null.model = c("sample.pool"),
         abundance.weighted=FALSE, runs = 5000, iterations = 5000)
