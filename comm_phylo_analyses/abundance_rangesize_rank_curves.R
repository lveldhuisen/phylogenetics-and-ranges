#generate curves for rank abundance and rank range size 

library(tidyverse)
library(patchwork)

#bring in data
all_df <- read.csv("big_results_all.csv")
all_df$Site[all_df$Site=="Road"]<-"Low elevation (2815 m)"
all_df$Site[all_df$Site=="Pfeiler"]<-"Middle elevation (3165 m)"
all_df$Site[all_df$Site=="PBM"]<-"High elevation (3380 m)"
all_df$Site <- factor(all_df$Site,levels=c("Low elevation (2815 m)",
                                "Middle elevation (3165 m)",
                                "High elevation (3380 m)")) 


#abundance---------
fig_ab <- ggplot(all_df, aes(x=Mean_abundance)) + 
  geom_density()+
  facet_grid(~Site)+
  theme_bw()+
  xlab("Mean abundance 2021-2022")
plot(fig_ab)

#range size--------
fig_rs <- ggplot(all_df, aes(x=AOO..km2.))+
  geom_density()+
  facet_grid(~Site)+
  scale_x_continuous(name = 'Area of occupancy (km2)', 
                     limits = c(0,500000),
                   breaks = c(10,500,1000,100000,500000))+
  theme_bw()
  
plot(fig_rs)

#combine
curve <- fig_ab / fig_rs +  
  plot_annotation(tag_levels = c('A'), tag_suffix = ')')
plot(curve)
