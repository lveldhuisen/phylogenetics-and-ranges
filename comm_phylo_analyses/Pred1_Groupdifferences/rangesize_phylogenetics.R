install.packages("fabricatr")
install.packages("ggpattern")
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
library(ggpubr)
library(ggpattern)

#group species by range size-------------------------------------------------
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")
all_df <- read.csv("results/results_all.csv")
big_results <- read.csv("big_results_all.csv")

#range size histograms
ggplot(data = big_results, aes(x = AOO..km2.)) +
  geom_histogram() + facet_wrap(~Site, scales = "free_x")

#import S&B phylogeny---------------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")

##import S&B18 tree and check data## 
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

##make community data matrix 
matrix <- read.table("comm_phylo_analyses/Pred2_addingspecies/comm_matrices/community_matrix_rangesize.txt", sep = "\t", header = T, row.names = 1)

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
phylo_df_rs <- read.csv("comm_phylo_analyses/Pred1_Groupdifferences/phylo_metrics_rangesize.csv")

##subset data by site-----
subset_road <- subset(phylo_df_rs, 
                      Site %in% c("Road"))

subset_pfeiler <- subset(phylo_df_rs, 
                         Site %in% c("Pfeiler"))
subset_pfeiler <- subset_pfeiler[-c(10,11,12), ]

subset_PBM <- subset(phylo_df_rs, 
                     Site %in% c("PBM"))
subset_PBM <- subset_PBM[-c(10,11,12), ]

##figures#####
PBM_fig <- ggplot(subset_PBM, aes(y=SES,alpha = Type, x=fct_relevel(Range_Size, c("small","medium","large")))) +
  geom_bar(position = "dodge",stat = "identity", color = "black")+
  xlab("Range size group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Range size group"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + 
  ylim(-2,2) + 
  ggtitle("High elevation (3380 m)")
  
plot(PBM_fig)

#colorful pattern fig, probbably wont use 
ggplot(subset_PBM, aes(x=fct_relevel(Range_Size, c("small","medium","large")),y=SES)) +
  geom_col_pattern(
    aes(pattern_type = Type, 
      pattern_fill = Type),
    pattern       = 'magick',
    pattern_key_scale_factor = 0.7,
    fill          = "white",
    colour        = 'black', 
    position = "dodge",
  ) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'cm')) +
  scale_pattern_type_discrete(choices = c('horizontal3', 'gray35', 'right45')) +
  coord_fixed(ratio = 1/2)+
  xlab("Range size cateogory")+
  ylim(-2,2)



pfeiler_fig <- ggplot(subset_pfeiler, aes(alpha = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity", color = "black") +
  xlab("Range size group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-2,2) + ggtitle("Middle elevation (3165 m)")
plot(pfeiler_fig)

road_fig <- ggplot(subset_road, aes(alpha = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity", color = "black") +
  xlab("Range size group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-2.5,1) + ggtitle("Low elevation (2815 m)")
plot(road_fig)

all_fig <- ggplot(phylo_df, aes(fill=Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-5,2) + 
  ggtitle("All sites together")+
  facet_wrap(~Site)
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

#code all sites together but without all sites combined
subset_allsites <- subset(phylometrics_df, 
                         Site %in% c("Pfeiler","PBM","Road"))

#test difference in phylo diversity between range size groups---------
phylometrics_df <- read.csv("results/phylo_metrics_rangesize.csv")

general_fig <- ggplot(subset_allsites, aes(fill = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  facet_wrap(~Site)
plot(general_fig) 

##all sites together###
kruskal.test(SES ~ Range_Size, data = phylo_df_rs)
anova(phylometrics_df)
pairwise.wilcox.test(phylo_df_rs$SES, phylo_df_rs$Range_Size,
                     p.adjust.method = "BH")

##road###
kruskal.test(SES ~ Range_Size, data = subset_road)
pairwise.wilcox.test(subset_road$SES, subset_road$Range_Size, p.adjust.method = "none")

##PBM###
kruskal.test(SES ~ Range_Size, data = subset_PBM)
pairwise.wilcox.test(subset_PBM$SES, subset_PBM$Range_Size)

##Pfeiler###
kruskal.test(SES ~ Range_Size, data = subset_pfeiler)


