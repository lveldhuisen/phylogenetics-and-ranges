# testing how rare species contribute to phylogenetic diversity in a subalpine plant community 
# code for all analyses in Veldhuisen, L. N., V. Zepeda, B. J. Enquist, and K. M. Dlugosch. 2025. Rare species do not disproportionately contribute to phylogenetic diversity in a subalpine plant community. American Journal of Botany e70061. https://doi.org/10.1002/ajb2.70061

ACCESS INFORMATION
1. MIT License
2. Data derived from other sources: phylogeny from Smith & Brown 2018 "ALLMB.tre" file


DATA & CODE FILE OVERVIEW

The "comm_phylo_analyses" folder has 4 subfolders with 8 toal code scripts for all analyses and figures for the manuscript "Rare species do not disproportionately contribute to phylogenetic diversity in a subalpine plant community." Outside of the "comm_phylo_analsyes" folder, this repository contains this README document, a folder for old, unsed analyses and empty folders for figures and results spreadsheets.


Data files and variables (listed in the order they appear in the repository, not in the order they appear in the manuscript)

    1. "Abundance_weighting" folder contains the file "abundance_weighting.R," which has the code for testing abundance and range size weighted phylogenetic diversity metrics. This file also has the code for Fig. 3. 
    
    2. "Phylogenetic_signal" folder contains the file "phylosignal.R," which has the code for testing for phylogenetic signal in abundance and range size. Also has code for Fig. 2. 

    3. "Removing_groups" folder contains the files "remove_groups10.R" and "removing_groups10_MPD_MNTD.R." "remove_groups10.R" has the code for removing sets of 10 species and testing Faith's PD and Fig. 5, while "removing_groups10_MPD_MNTD.R" has the same analyses for MPD and MNTD (Appendix 5, Fig. S5). 

    4.  "Removing_species" folder contains the files "dispersioncurve_singleremoval_abundance.R" and "dispersioncurve_singleremoval_rangesize.R" The abundace file has the code for Fig. 4, and the range size file has code for the same analysis with range size, which is in Appendix S3. 

    5. "abundance_rangesize_rank_curves.R" has the code for Figure 1. 

    6. "rangesize_calculations.R" has the code to calculate range size for a species. I did not functionalize this, so you have to input each individual species name to get the GBIF data and corresponding range size. 



Code scripts and workflow for manuscript analyses - code should be run in order listed here. 

    1. First, calculate the range size of each species in file #6, "rangesize_calculations.R." 

    2. Then, visualize the distributions of abundance and range size (Fig. 1) in file #5, "abundance_rangesize_rank_curves.R." 

    3. Calculate phylogenetic signal for abundance and range size using the script in folder #5, "phylosignal.R." Generate the phylogenies in Fig. 2 with this same script. 

    4. Calculate weighted (by abundance and range) and unweighted phylogenetic diversity metrics using the script in folder #1, "abundance_weighting.R." Generate Fig. 3 in this script also. 

    5. Test the impact of removing individual species by abundance with the script in folder #4, "dispersioncurve_singleremoval_abundance.R." Generate Fig. 4 and test for autocorrelation between points (Appendix 4, Fig. S4) in this same script. 

    6. Test the impact of removing sets of 10 species grouped by abundance with the script "remove_groups10.R" in folder #3. Generate Fig. 5 with the same script. This script only assesses Faith's PD, not MPD or MNTD (see below). 


Code scripts and workflow for supplemental analyses and figures - code should be run in order listed here. 

    1. Test the impact of removing individual species by range size with the script in folder #4, "dispersioncurve_singleremoval_rangesize.R." Generate Fig. S3 in Appendix 3 in this same script. 

    2. Test the impact of removing sets of 10 species grouped by abundance with the script "remove_groups10_MPD_MNTD.R" in folder #3. Generate Fig. S5 in Appendix 5 with the same script.


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


 
