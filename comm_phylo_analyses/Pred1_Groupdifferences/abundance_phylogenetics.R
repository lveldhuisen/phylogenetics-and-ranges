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

##figures#####
PBM_fig <- ggplot(subset_PBM, aes(y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")), alpha = Type)) + 
  geom_bar(position = "dodge",stat = "identity", color = "black") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Abundance group"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-1,2) +
  ggtitle("PBM - high elevation")
plot(PBM_fig)

pfeiler_fig <- ggplot(subset_pfeiler, aes(alpha=Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-1,1.5) + ggtitle("Pfeiler - middle elevation")
plot(pfeiler_fig)

road_fig <- ggplot(subset_road, aes(alpha=Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance gorup") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-2,1) + ggtitle("Road - low elevation")
plot(road_fig)

all_fig <- ggplot(phylo_df_a, aes(alpha=Type, y=SES, x=fct_relevel(Abundance_group, c("low","medium","high")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + 
  ylim(-5,2) + 
  ggtitle("All sites together")+
  facet_wrap(~Site)
plot(all_fig)

#fig without range size categories
general_fig <- ggplot(phylo_df, aes(alpha=Type, y=SES, x=fct_relevel(Site, c("all")))) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Abundance group") + 
  theme_light() + 
  guides(fill=guide_legend(title="Phylogenetic metric"))+
  scale_fill_manual(values=c("#c385b3",
                             "#cdd870",
                             "#4ea6c4"))  + ylim(-3,3) 
plot(general_fig)  

#test difference in phylo diversity between range size groups---------
phylo_df_a <- read.csv("comm_phylo_analyses/Pred1_Groupdifferences/abundance_phylo_metrics.csv")
phylo_df_a = subset(phylo_df_a, select = -c(Value) )

##all sites together###
kruskal.test(SES ~ Abundance_group, data = phylo_df_a)
pwc <- phylo_df_a %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni") 
pwc

all_facet_site <- ggboxplot(phylo_df_a, x = "Abundance_group", y = "SES",
          color = "Abundance_group", palette = c("#00AFBB", "#E7B800", "#FC4E07", "grey"),
          order = c("low", "medium", "high", "ALL"),
          ylab = "SES", xlab = "Abundance group") +
  facet_wrap(~Site)+
  geom_hline(yintercept = 0, linetype = "dotted")

plot(all_facet_site)

##all sites grouped by abundance, no facet wrap###
All_together <- ggboxplot(phylo_df_a, x = "Abundance_group", y = "SES",
                            color = "Abundance_group", palette = c("#00AFBB", "#E7B800", "#FC4E07", "grey"),
                            order = c("low", "medium", "high", "ALL"),
                            ylab = "SES", xlab = "Abundance group") +
  geom_hline(yintercept = 0, linetype = "dotted")+
  stat_pvalue_manual(pwc, hide.ns = TRUE, y.position = 2.1, label = "p.adj")

plot(All_together)

##road###
kruskal.test(SES ~ Abundance_group, data = subset_road)
pairwise.wilcox.test(subset_road$SES, subset_road$Abundance_group, p.adjust.method = "BH")

pwc <- subset_road %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni") 
pwc

Road_boxplot <- ggboxplot(subset_road, x = "Abundance_group", y = "SES",
          color = "Abundance_group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("low", "medium", "high"),
          ylab = "SES", xlab = "Abundance_group") +
  stat_pvalue_manual(pwc, hide.ns = TRUE, y.position = 0.55, label = "p.adj")+
  geom_hline(yintercept = 0, linetype = "dotted")

plot(Road_boxplot)

##PBM###
kruskal.test(SES ~ Abundance_group, data = subset_PBM)
pwc <- subset_PBM %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni") 
pwc

PBM_boxplot <- ggboxplot(subset_PBM, x = "Abundance_group", y = "SES",
          color = "Abundance_group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("low", "medium", "high"),
          ylab = "SES", xlab = "Abundance_group")+
  stat_pvalue_manual(pwc, hide.ns = TRUE, y.position = 0.55, label = "p.adj")+
  geom_hline(yintercept = 0, linetype = "dotted")

plot(PBM_boxplot)

##Pfeiler###
kruskal.test(SES ~ Abundance_group, data = subset_pfeiler)
pwc <- subset_pfeiler %>% 
  dunn_test(SES ~ Abundance_group, p.adjust.method = "bonferroni") 
pwc

Pfeiler_boxplot <- ggboxplot(subset_pfeiler, x = "Abundance_group", y = "SES",
          color = "Abundance_group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("low", "medium", "high"),
          ylab = "SES", xlab = "Abundance_group")+
  stat_pvalue_manual(pwc, hide.ns = TRUE, y.position = 0.55, label = "p.adj")+
  geom_hline(yintercept = 0, linetype = "dotted")

plot(Pfeiler_boxplot)

combined <- Road_boxplot + Pfeiler_boxplot + PBM_boxplot + theme(legend.position = "none")
plot(combined)
