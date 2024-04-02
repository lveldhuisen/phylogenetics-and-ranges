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
library(patchwork)
library(egg)

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

all_fig_rs <- ggplot(phylo_df_rs, aes(fill = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  ggtitle("All sites combined")
  

plot(all_fig_rs)

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
PBM_fig_rs <- ggplot(subset_PBM, aes(fill = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  ggtitle("High elevation (3380 m)")+
  stat_compare_means(method = "kruskal.test")
  
plot(PBM_fig_rs)


pfeiler_fig_rs <- ggplot(subset_pfeiler, aes(fill = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  ggtitle("Middle elevation (3165 m)")
plot(pfeiler_fig_rs)

road_fig_rs <- ggplot(subset_road, aes(fill = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  ggtitle("Low elevation (2815 m)")
plot(road_fig_rs)

#use patchwork to combine all figures together
combined_fig_rs <- (all_fig | road_fig | pfeiler_fig | PBM_fig) + 
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = 'collect')+
  plot_layout(axes = "collect")

plot(combined_fig_rs)

#code all sites together but without all sites combined
subset_allsites <- subset(phylometrics_df, 
                         Site %in% c("Pfeiler","PBM","Road"))


#test difference in phylo diversity between range size groups---------
phylometrics_df <- read.csv("results/phylo_metrics_rangesize.csv")

phylometrics_df$Site <- factor(phylometrics_df$Site, levels = c("All", "Road", "Pfeiler","PBM"))

general_fig_withcombined <- ggplot(phylometrics_df, aes(fill = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  facet_wrap(factor(Site, levels = c("All", "Road", "Pfeiler","PBM"), 
                    labels = c("All", "Low (2815 m)", "Middle (3165 m)","High (3380 m)")) ~., ncol = 4) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "grey")+
  geom_hline(yintercept=-1.3, linetype="dashed", color = "grey")
plot(general_fig_withcombined) 

#plot all three sites but not with combined 
fig_individualsites <- ggplot(subset_allsites, aes(fill = Type, y=SES, x=fct_relevel(Range_Size, c("small","medium","large")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Range size") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  facet_wrap(~Site)+
  geom_hline(yintercept=1.3, linetype="dashed", color = "grey")+
  geom_hline(yintercept=-1.3, linetype="dashed", color = "grey")+
  tag_facet()
  
plot(fig_individualsites) 

##stats--------
##all sites together###
kruskal.test(SES ~ Range_Size, data = phylo_df_rs)
anova(phylometrics_df)
pairwise.wilcox.test(phylo_df_rs$SES, phylo_df_rs$Range_Size,
                     p.adjust.method = "BH")
###make dataframe to add p values to figures###
pairwise_tests_allsites <- phylo_df_rs %>% 
  dunn_test(SES ~ Range_Size, p.adjust.method = "bonferroni") 

pairwise_tests_allsites$y.position <- c(2, 2, 2)

##road###
kruskal.test(SES ~ Range_Size, data = subset_road)
subset_road %>% 
  dunn_test(SES ~ Range_Size, p.adjust.method = "bonferroni") 

##PBM###
kruskal.test(SES ~ Range_Size, data = subset_PBM)
subset_PBM %>% 
  dunn_test(SES ~ Range_Size, p.adjust.method = "bonferroni") 


##Pfeiler###
kruskal.test(SES ~ Range_Size, data = subset_pfeiler)
subset_pfeiler %>% 
  dunn_test(SES ~ Range_Size, p.adjust.method = "bonferroni") 



