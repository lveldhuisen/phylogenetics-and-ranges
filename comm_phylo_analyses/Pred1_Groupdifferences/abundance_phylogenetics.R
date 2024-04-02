install.packages("rstatix")

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
library(rstatix)
library(tidyverse)
library(patchwork)

#make histograms of abundance
ggplot(data = big_results, aes(x = Mean_abundance)) +
  geom_histogram()+ facet_wrap(~Site, scales = "free_x")



#import S&B phylogeny---------------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")

##import S&B18 tree and check data## 
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

##make community data matrix 
abundance_matrix <- read.table("comm_phylo_analyses/Pred1_Groupdifferences/abundance_groups_matrix.txt", sep = "\t", header = T, row.names = 1)

#import S&B phylogeny---------------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")

##import S&B18 tree and check data## 
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

pruned.tree <- treedata(SBtree, unlist(abundance_matrix[16,abundance_matrix[16,]>0]), warnings = F)$phy
plot(pruned.tree)

###SES###############################
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")
library(picante)

#standard effect size using picante 
pd <- ses.pd(abundance_matrix, pruned.tree, null.model = c("sample.pool"),
       runs = 5000, include.root=TRUE)

as.table(pd)

#make distance matrix for MPD and MNTD
dist.mat <- cophenetic(pruned.tree)

mpd <- ses.mpd(abundance_matrix, dist.mat, null.model = c("sample.pool"),
        abundance.weighted = FALSE, runs = 5000, iterations = 5000)

mntd <- ses.mntd(abundance_matrix, dist.mat, null.model = c("sample.pool"),
         abundance.weighted=FALSE, runs = 5000, iterations = 5000)


#figures###
phylo_df_a <- read.csv("comm_phylo_analyses/Pred1_Groupdifferences/abundance_phylo_metrics.csv")

all_sites_a <- ggplot(phylo_df_a, aes(fill = Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) 

plot(all_sites_a)

##subset data by site-----
subset_road <- subset(phylo_df_a, 
                      Site %in% c("Road"))
subset_road <- subset_road[-c(10,11,12), ] 

subset_pfeiler <- subset(phylo_df_a, 
                         Site %in% c("Pfeiler"))
subset_pfeiler <- subset_pfeiler[-c(10,11,12), ]

subset_PBM <- subset(phylo_df_a, 
                     Site %in% c("PBM"))
subset_PBM <- subset_PBM[-c(10,11,12), ]

## individual site figures#####
PBM_fig_a <- ggplot(subset_PBM, aes(fill = Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2)
plot(PBM_fig_a)

pfeiler_fig_a <- ggplot(subset_pfeiler, aes(fill = Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) 
plot(pfeiler_fig_a)

road_fig_a <- ggplot(subset_road, aes(fill = Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) 
  
plot(road_fig_a)

#use patchwork to make combined fig
combined_fig_a <- all_fig_a | road_fig_a | pfeiler_fig_a | PBM_fig_a +
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = 'collect')+
  plot_layout(axes = "collect")

plot(combined_fig_a)

#abundance and range size together figure
combined_fig <- all_fig_rs + road_fig_rs + pfeiler_fig_rs + PBM_fig_rs | all_fig_a + road_fig_a + pfeiler_fig_a + PBM_fig_a +
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = 'collect')+
  plot_layout(axes = "collect", ncol = 2)

plot(combined_fig)

##combine into one big fig with patchwork-------
fig1 <- combined_fig_a / combined_fig_rs +  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = 'collect')+
  plot_layout(axes = "collect")

plot(fig1)

#code for all sites together in one fig, use these######
#subset to get rid of metrics for sites as wholes
phylo_df_a <- subset(phylo_df_a, 
                            Abundance_group %in% c("low","medium","high"))


fig_abundance_withcombined <- ggplot(phylo_df_a, aes(fill = Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-2.5,2) +
  facet_wrap(factor(Site, levels = c("All", "Road", "Pfeiler","PBM"), 
                    labels = c("All", "Low (2815 m)", "Middle (3165 m)","High (3380 m)")) ~., ncol = 4) +
  geom_hline(yintercept=1.5, linetype="dashed", color = "grey")+
  geom_hline(yintercept=-1.5, linetype="dashed", color = "grey")

plot(fig_abundance_withcombined) 


#plot all sites in one set of code without all combined
subset_allsites_a <- subset(phylo_df_a, 
                          Site %in% c("Pfeiler","PBM","Road"))


fig_abundance_indsites <- ggplot(subset_allsites_a, aes(fill = Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-1.7,2) +
  facet_wrap(~Site) +
  geom_hline(yintercept=1.2, linetype="dashed", color = "grey")+
  geom_hline(yintercept=-1.4, linetype="dashed", color = "grey")

plot(fig_abundance_indsites)

#test difference in phylo diversity between abundance groups---------
phylo_df_a <- read.csv("comm_phylo_analyses/Pred1_Groupdifferences/abundance_phylo_metrics.csv")
phylo_df_a = subset(phylo_df_a, select = -c(Value) )

##all sites together###
kruskal.test(SES ~ Abundance_group, data = phylo_df_a)
pwc <- phylo_df_a %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni") 
pwc

##road###
kruskal.test(SES ~ Abundance_group, data = subset_road)

subset_road %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni")

##PBM###
kruskal.test(SES ~ Abundance_group, data = subset_PBM)
subset_PBM %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni") 

##Pfeiler###
kruskal.test(SES ~ Abundance_group, data = subset_pfeiler)
subset_pfeiler %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni") 



