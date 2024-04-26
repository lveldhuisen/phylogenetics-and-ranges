#update community data matrices with abundance data
library(picante)
library(ape)
library(geiger)
library(tidyverse)

#abundance data----------
##PBM####
##make community data matrix#### 
PBM_abundance_matrix_weighted <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/PBM_abundance_commmatrix_removal_weighted.txt", sep = "\t", header = T, row.names = 1)

MPD_PBM_abundance_removal_weighted <- ses.mpd(PBM_abundance_matrix_weighted, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted = TRUE, runs = 5000, iterations = 5000)




MPD_PBM_abundance_removal_weighted <- MPD_PBM_abundance_removal_weighted[-c(32,33),]
MPD_PBM_abundance_removal_weighted$Abundance_rank <- c(1:31)

test <- ggplot(data= MPD_PBM_abundance_removal_weighted) + 
  geom_segment( aes(x=Abundance_rank, xend=Abundance_rank, y=1.58, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Abundance_rank, y=mpd.obs.z), size = 2) +
  xlab("Abundance rank of removed species (most to least)") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(1, 1.5,2),limits=c(1, 2)) +
  theme_classic(14) +
  geom_hline(yintercept = 1.58, col = "lightgrey") +
  xlim(0,32) 
plot(test)


data(phylocom)
mpd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)

##Pfeiler####
##Road####

#range size data----------
##PBM####
pruned.tree <- treedata(SBtree, unlist(PBM_rangesize_matrix_weighted[33,PBM_rangesize_matrix_weighted[33,]>0]), warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

PBM_rangesize_matrix_weighted <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/PBM_rangesize_removal_weighted.txt", sep = "\t", header = T, row.names = 1)

MPD_PBM_range_removal_weighted <- ses.mpd(PBM_rangesize_matrix_weighted, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted = TRUE, runs = 5000, iterations = 5000)


##Pfeiler####
##Road####