#calculate phylogenetic diversity of abundance and range size
library(tidyverse)
library(dplyr)
library(picante)
library(geiger)
library(ape)

#import S&B phylogeny--------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/RMBL phylogeny/Smith&Brown18")
SBtree <- read.tree(file = "ALLMB.tre")
write.tree(SBtree)
is.rooted(SBtree)