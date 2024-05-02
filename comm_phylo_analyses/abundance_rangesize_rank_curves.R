#generate curves for rank abundance and rank range size 

library(tidyverse)

#bring in data
all_df <- read.csv("big_results_all.csv")

#abundance---------
fig_ab <- ggplot(all_df, aes(x=Mean_abundance)) + 
  geom_density()+
  facet_grid(~Site)+
  theme_bw()
plot(fig_ab)

#range size--------
fig_rs <- ggplot(all_df, aes(x=AOO..km2.))+
  geom_density()+
  facet_grid(~Site)+
  xlab('Area of occupancy (km2)')+
  ylab("")
plot(fig_rs)
