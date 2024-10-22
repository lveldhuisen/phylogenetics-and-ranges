# testing how rare species contribute to phylogenetic diversity in a subalpine plant community 

ACCESS INFORMATION
1. MIT License
2. Data derived from other sources: phylogeny from Smith & Brown 2018 "ALLMB.tre" file


DATA & CODE FILE OVERVIEW

Code for all analyses is in the folder titled "comm_phylo_analyses."

This data repository consists of 5 folders containing 7 code scripts (not including the 6 unused analysis scripts) within the folder titled "comm_phylo_analyses," and this README document.


Data files and variables (listed in the order they appear in the repository, not in the order they appear in the manuscript)

    1. "Removing_species" folder contains the files "dispersioncurve_singleremoval_abundance.R" and "dispersioncurve_singleremoval_rangesize.R" The abundace file has the code for Fig. 4, and the range size file has code for the same analysis with range size, which is in the supporting information. 
    
    2. "abundance_rangesize_rank_curves.R" has the code for Figure 1. 

    3. "Removing_groups" folder contains the file "remove_groups10.R," which has the code for removing sets of 10 species and Fig. 5

    4. "Abundance_weighting" folder contains the file "abundance_weighting.R," which has the code for testing abundance and range size weighted phylogenetic diversity metrics. This file also has the code for Fig. 3. 

    5. "Phylogenetic_signal" folder contains the file "phylosignal.R," which has the code for testing for phylogenetic signal in abundance and range size. Also has code for Fig. 2. 

    6. "rangesize_calculations.R" has the code to calculate range size for a species. I did not functionalize this, so you have to input each individual species name to get the GBIF data and corresponding range size. 



Code scripts and workflow - code should be run in order listed here. 

    1. First, calculate the range size of each species in file #7, "rangesize_calculations.R." 

    2. Then, visualize the distributions of abundance and range size (Fig. 1) in file #2, "abundance_rangesize_rank_curves.R." 

    3. Calculate phylogenetic signal for abundance and range size using file #6, "phylosignal.R." 

    4. Calculate weighted (by abundance and range) and unweighted phylogenetic diversity metrics using file #5, "abundance_weighting.R." 

    5. Test the impact of removing individual species by abundance with file #1, "dispersioncurve_singleremoval_abundance.R." 

    6. Test the impact of removing sets of 10 species grouped by abundance with file #4, "remove_groups10.R." 


SOFTWARE VERSIONS

R: 4.2.3,

biomod2: 4.2-4

RBIF: 3.7.7

red: 1.6.1

picante: 1.8.2

geiger: 2.0.11

ape: 5.8 

phytools: 2.1-1

tidverse: 2.0.0

viridis: 0.6.5


REFERENCES

Smith, S. A., and J. W. Brown. 2018. Constructing a broadly inclusive seed plant phylogeny. American journal of Botany 105:302–314.


 
