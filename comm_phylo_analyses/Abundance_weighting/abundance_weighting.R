#test how weighting by abundance and range size changes community phylogenetic
#diversity metrics, code for Figure 3

library(picante)
library(ape)
library(geiger)
library(tidyverse)

#calculate MPD and MNTD weighted and unweighted--------

#import S&B phylogeny
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")

##import S&B18 tree and prune## 
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)

#prune tree
pruned.tree <- treedata(SBtree, unlist(community_matrix[4,community_matrix[4,]>0]), warnings = F)$phy
plot(pruned.tree)

#check species in tree
specieslist <- SBtree$tip.label
specieslist <- as.data.frame(specieslist)

##bring in community data matrices####

#unweighted
community_matrix <- read.table("comm_phylo_analyses/Abundance_weighting/community_matrix_unweighted.txt", sep = "\t", header = T, row.names = 1)

#weighted by abundance
abundance_matrix_weighted <- read.table("comm_phylo_analyses/Abundance_weighting/community_matrix_weighted_abundance.txt", sep = "\t", header = T, row.names = 1)

#weighted by range size
rangesize_matrix <- read.table("comm_phylo_analyses/Abundance_weighting/community_matrix_weighted_rangesize.txt", sep = "\t", header = T, row.names = 1)

##unweighted####

###mpd####
mpd_unweighted <- ses.mpd(community_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"), abundance.weighted = FALSE, runs = 5000, iterations = 5000)

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
                        select = -c(ntaxa,mntd.obs,mntd.rand.mean,mntd.rand.sd,
                                    mntd.obs.rank,runs) ) #remove unnecessary columns
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
mpd_weighted_a <- ses.mpd(abundance_matrix_weighted, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted = TRUE, runs = 5000, iterations = 5000)

###format dataframe#####
mpd_weighted_a = subset(mpd_weighted_a, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns
mpd_weighted_a <- mpd_weighted_a[-c(4),]
mpd_weighted_a <- mpd_weighted_a %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

mpd_weighted_a$Type <- c("MPD") #add column for metric type 
mpd_weighted_a$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                         "High elevation (3380 m)") #add column for site name 
mpd_weighted_a$Weighting <- c("Weighted") #add column for weighting

mntd_unweighted <- ses.mntd(community_matrix, cophenetic(pruned.tree), 
                            null.model = c("sample.pool"),
                            abundance.weighted=FALSE, runs = 5000, iterations = 5000)

###mntd####
mntd_weighted_a <- ses.mntd(abundance_matrix_weighted, cophenetic(pruned.tree), 
                            null.model = c("sample.pool"),
                            abundance.weighted=TRUE, runs = 5000, iterations = 5000)

###format dataframe#####
mntd_weighted_a = subset(mntd_weighted_a, 
                         select = -c(ntaxa,mntd.obs,mntd.rand.mean,
                                     mntd.rand.sd,mntd.obs.rank,
                                     runs) ) #remove unnecessary columns
mntd_weighted_a <- mntd_weighted_a[-c(4),]
mntd_weighted_a <- mntd_weighted_a %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

mntd_weighted_a$Type <- c("MNTD") #add column for metric type 
mntd_weighted_a$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                          "High elevation (3380 m)") #add column for site name 
mntd_weighted_a$Weighting <- c("Weighted") #add column for weighting

###combine abundance weighted dataframes####
all_weighted_a <- rbind(mntd_weighted_a,mpd_weighted_a)
all_df <- rbind(all_unweighted,all_weighted_a)

##range size####

###mpd####
mpd_weighted_rs <- ses.mpd(rangesize_matrix, cophenetic(pruned.tree), null.model = c("sample.pool"),abundance.weighted = TRUE, runs = 5000, iterations = 5000)

###format dataframe#####
mpd_weighted_rs = subset(mpd_weighted_rs, select = -c(ntaxa,mpd.obs,mpd.rand.mean,mpd.rand.sd,mpd.obs.rank,runs) ) #remove unnecessary columns
mpd_weighted_rs <- mpd_weighted_rs[-c(4),]
mpd_weighted_rs <- mpd_weighted_rs %>% 
  rename(SES = mpd.obs.z,
         P_value = mpd.obs.p) #rename columns to match other datasets 

mpd_weighted_rs$Type <- c("MPD") #add column for metric type 
mpd_weighted_rs$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                         "High elevation (3380 m)") #add column for site name 
mpd_weighted_rs$Weighting <- c("Weighted") #add column for weighting

###mntd####
mntd_weighted_rs <- ses.mntd(rangesize_matrix, cophenetic(pruned.tree), 
                            null.model = c("sample.pool"),
                            abundance.weighted=TRUE, runs = 5000, iterations = 5000)

###format dataframe#####
mntd_weighted_rs = subset(mntd_weighted_rs, 
                         select = -c(ntaxa,mntd.obs,mntd.rand.mean,
                                     mntd.rand.sd,mntd.obs.rank,
                                     runs) ) #remove unnecessary columns
mntd_weighted_rs <- mntd_weighted_rs[-c(4),]
mntd_weighted_rs <- mntd_weighted_rs %>% 
  rename(SES = mntd.obs.z,
         P_value = mntd.obs.p) #rename columns to match other datasets 

mntd_weighted_rs$Type <- c("MNTD") #add column for metric type 
mntd_weighted_rs$Site <- c("Low elevation (2815 m)","Middle elevation (3165 m)",
                          "High elevation (3380 m)") #add column for site name 
mntd_weighted_rs$Weighting <- c("Weighted") #add column for weighting

#combine range size and unweighted datasets
all_weighted_rs <- rbind(mpd_weighted_rs,mntd_weighted_rs)
combo_rs <- rbind(all_unweighted,all_weighted_rs)

#make figure to compare diversity with weighting----------
###abundance fig#####
all_df$Site <- factor(all_df$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)"))  #order sites low to high elevation

fig_weighting_a <- ggplot(all_df, aes(fill = Weighting, y=SES, x=Type)) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Phylogenetic metric") + 
  ylab("Standard effect size")+
  theme_light(base_size = 20) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8) + 
  ylim(-1,2) +
  facet_wrap(Site ~ .)+
  theme(strip.text = element_text(color = "black"))

plot(fig_weighting_a)

###range size fig####

combo_rs$Site <- factor(combo_rs$Site, levels = c("Low elevation (2815 m)","Middle elevation (3165 m)","High elevation (3380 m)")) #order sites low to high elevation

fig_weighting_rs <- ggplot(combo_rs, aes(fill = Weighting, y=SES, x=Type)) + 
  geom_bar(position = "dodge",stat = "identity") +
  xlab("Phylogenetic metric") + 
  ylab("Standard effect size")+
  theme_light(base_size = 20) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8) + 
  ylim(-1,2) +
  facet_wrap(Site ~ .)+
  theme(strip.text = element_text(color = "black"))

plot(fig_weighting_rs)

#combine to make one figure with range size and abundance 
weighting_fig <- fig_weighting_a / fig_weighting_rs + 
  plot_annotation(tag_levels = c('A'), tag_suffix = ')')+
  plot_layout(guides = 'collect')
plot(weighting_fig)
