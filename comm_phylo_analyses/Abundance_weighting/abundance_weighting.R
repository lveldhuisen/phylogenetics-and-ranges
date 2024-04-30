#update community data matrices with abundance data
library(picante)
library(ape)
library(geiger)
library(tidyverse)

#calculate MPD and MNTD weighted and unweighted--------
#import S&B phylogeny
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")

##import S&B18 tree and check data## 
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#prune tree
pruned.tree <- treedata(SBtree, unlist(community_matrix[4,community_matrix[4,]>0]), warnings = F)$phy
plot(pruned.tree)

##make community data matrices 
#unweighted comparison
community_matrix <- read.table("comm_phylo_analyses/Abundance_weighting/community_matrix_unweighted.txt", 
                                        sep = "\t", header = T, row.names = 1)

#weighted by abundance
abundance_matrix_weighted <- read.table("comm_phylo_analyses/Abundance_weighting/community_matrix_weighted_abundance.txt", 
                                        sep = "\t", header = T, row.names = 1)

#weighted by range size


##unweighted####
###mpd####
mpd_unweighted <- ses.mpd(community_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),
               abundance.weighted = FALSE, runs = 5000, iterations = 5000)
###format dataframe#####
mpd_unweighted = subset(mpd_unweighted, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns
mpd_unweighted <- mpd_unweighted[-c(4),]
mpd_unweighted <- mpd_unweighted %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

mpd_unweighted$Type <- c("MPD") #add column for metric type 
mpd_unweighted$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                         "High elevation (3380 m)") #add column for site name 
mpd_unweighted$Weighting <- c("Unweighted") #add column for weighting

###mntd###
mntd_unweighted <- ses.mntd(community_matrix, cophenetic(pruned.tree), 
                            null.model = c("sample.pool"),
                 abundance.weighted=FALSE, runs = 5000, iterations = 5000)

###format dataframe#####
mntd_unweighted = subset(mntd_unweighted, 
                        select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns
mntd_unweighted <- mntd_unweighted[-c(4),]
mntd_unweighted <- mntd_unweighted %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

mntd_unweighted$Type <- c("MNTD") #add column for metric type 
mntd_unweighted$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                          "High elevation (3380 m)") #add column for site name 
mntd_unweighted$Weighting <- c("Unweighted") #add column for weighting

###combine#####
all_unweighted <- rbind(mpd_unweighted,mntd_unweighted)

##abundance####
###mpd####
mpd_weighted_a <- ses.mpd(abundance_matrix_weighted, cophenetic(pruned.tree), null.model = c("sample.pool"),
                          abundance.weighted = TRUE, runs = 5000, iterations = 5000)
###format dataframe#####
mpd_weighted_a = subset(mpd_weighted_a, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns
mpd_weighted_a <- mpd_weighted_a[-c(4),]
mpd_weighted_a <- mpd_weighted_a %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

mpd_weighted_a$Type <- c("MPD") #add column for metric type 
mpd_weighted_a$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                         "High elevation (3380 m)") #add column for site name 
mpd_weighted_a$Weighting <- c("Weighted - abundance") #add column for weighting

mntd_unweighted <- ses.mntd(community_matrix, cophenetic(pruned.tree), 
                            null.model = c("sample.pool"),
                            abundance.weighted=FALSE, runs = 5000, iterations = 5000)
###mntd####
mntd_weighted_a <- ses.mntd(abundance_matrix_weighted, cophenetic(pruned.tree), 
                            null.model = c("sample.pool"),
                            abundance.weighted=TRUE, runs = 5000, iterations = 5000)

###format dataframe#####
mntd_weighted_a = subset(mntd_weighted_a, 
                         select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,mntd.obs.rank,runs) ) #remove unnecessary columns
mntd_weighted_a <- mntd_weighted_a[-c(4),]
mntd_weighted_a <- mntd_weighted_a %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

mntd_weighted_a$Type <- c("MNTD") #add column for metric type 
mntd_weighted_a$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                          "High elevation (3380 m)") #add column for site name 
mntd_weighted_a$Weighting <- c("Weighted - abundance") #add column for weighting

###combine abundance weighted dataframes####
all_weighted_a <- rbind(mntd_weighted_a,mpd_weighted_a)
all_df <- rbind(all_unweighted,all_weighted_a)

##range size####


#make figure to compare diversity with weighting----------
all_df$Site <- factor(all_df$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"))

ggplot(all_df, aes(fill = Weighting, y=SES, x=Type)) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Phylogenetic metric") + 
  ylab("Standard effect size")+
  theme_light() + 
  scale_fill_viridis_d(begin = 0.1) + 
  ylim(-1,2) +
  facet_wrap(Site ~ .)


#abundance weighting for pred 2 removing species, not used--------------------
#abundance data----------
##PBM####
##make community data matrix#### 
PBM_abundance_matrix_weighted <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/PBM_abundance_commmatrix_removal_weighted.txt", sep = "\t", header = T, row.names = 1)

MPD_PBM_abundance_removal_weighted <- ses.mpd(PBM_abundance_matrix_weighted, 
                                              cophenetic(pruned.tree), null.model = c("sample.pool"), 
                                              abundance.weighted = TRUE, runs = 5000, iterations = 5000)


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
  xlim(0,32)+
  ggtitle("High elevation site")
plot(test)

##Pfeiler####
##Road####

#range size data----------
##PBM####
pruned.tree.rs <- treedata(SBtree, unlist(PBM_rangesize_matrix_weighted[32,PBM_rangesize_matrix_weighted[32,]>0]), 
                           warnings = F)$phy
write.tree(pruned.tree)
plot(pruned.tree)
is.rooted(pruned.tree)

PBM_rangesize_matrix_weighted <- read.table("comm_phylo_analyses/Pred3_removingspecies/comm_matrices/PBM_rangesize_removal_weighted.txt", sep = "\t", header = T, row.names = 1)

MPD_PBM_range_removal_weighted <- ses.mpd(PBM_rangesize_matrix_weighted, 
                                          cophenetic(pruned.tree.rs), 
                                          null.model = c("sample.pool"),
                                          abundance.weighted = TRUE, 
                                          runs = 5000, iterations = 5000)


MPD_PBM_range_removal_weighted <- MPD_PBM_range_removal_weighted[-c(32,33),]
MPD_PBM_range_removal_weighted$Range_rank <- c(1:31)

test <- ggplot(data= MPD_PBM_range_removal_weighted) + 
  geom_segment( aes(x=Range_rank, xend=Range_rank, y=0.98, yend=mpd.obs.z), color="grey")+
  geom_point(mapping = aes(x=Range_rank, y=mpd.obs.z), size = 2) +
  xlab("Range size rank of removed species (most to least)") +
  ylab("SES MPD") +
  scale_y_continuous(name="SES MPD", breaks = c(0.5,1, 1.5),limits=c(0.5, 1.5)) +
  theme_classic(14) +
  geom_hline(yintercept = 0.98, col = "lightgrey") +
  xlim(0,32)+
  ggtitle("High elevation site")
plot(test)
##Pfeiler####
##Road####