#generate curves for rank abundance and rank range size 
#code for Figure 1

library(tidyverse)
library(patchwork)
library(ggExtra)
library(RColorBrewer)
library(cowplot)
library(ggpubr)


#bring in data
all_df <- read.csv("comm_phylo_analyses/results_all.csv")

#rename sites with elevations
all_df$Site[all_df$Site=="Road"]<-"Low elevation (2815 m)"
all_df$Site[all_df$Site=="Pfeiler"]<-"Middle elevation (3165 m)"
all_df$Site[all_df$Site=="PBM"]<-"High elevation (3380 m)"
all_df$Site <- factor(all_df$Site,levels=c("Low elevation (2815 m)",
                                "Middle elevation (3165 m)",
                                "High elevation (3380 m)")) 

#abundance figure---------
fig_ab <- ggplot(all_df, aes(x=Mean_abundance, fill = Site, alpha = 0.2)) + 
  geom_density()+
  scale_fill_manual(values = c("#AADC32FF","#21908CFF","#440154"))+
  facet_wrap(~Site, scales = "free")+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        strip.background = element_blank(), 
        strip.text.x.top = element_blank())+
  xlab("Mean abundance 2021-2022")

plot(fig_ab)

#range size figure--------
fig_rs <- ggplot(all_df, aes(x=AOO..km2., fill = Site, alpha = 0.2))+
  geom_density()+
  facet_wrap(~Site, scales = "free")+
  scale_fill_manual(values = c("#AADC32FF","#21908CFF","#440154"))+
  scale_x_continuous(name = 'Area of occupancy (km2)')+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        strip.background = element_blank(), 
        strip.text.x.top = element_blank())
  
plot(fig_rs)

#combine into one figure with both range size and abundance
curve <- fig_ab / fig_rs +  
  plot_annotation(tag_levels = c('A'), tag_suffix = ')')
plot(curve)

#Plot range size and abundance------

#figures for each site separately
df_low <- all_df %>%
  filter(Site %in% "Low elevation (2815 m)")

rs_ab_fig_low <- ggplot(df_low, aes(x=AOO..km2., y = Mean_abundance))+
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))+
  ylab("Mean abundance 2021-2022")+
  xlab(" ")

low_fig <- ggMarginal(rs_ab_fig_low, type = "density", 
                      fill = "#AADC32FF", alpha = 0.5, color = "#AADC32FF")
plot(low_fig)

#mid elevation
df_mid <- all_df %>%
  filter(Site %in% "Middle elevation (3165 m)")

rs_ab_fig_mid <- ggplot(df_mid, aes(x=AOO..km2., y = Mean_abundance))+
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))+
  ylab(" ")+
  xlab("Area of occupancy (km2)")

mid_fig <- ggMarginal(rs_ab_fig_mid, type = "density", fill = "#21908CFF", 
                      alpha = 0.5, color = "#21908CFF")


#high
df_high <- all_df %>%
  filter(Site %in% "High elevation (3380 m)")

rs_ab_fig_high <- ggplot(df_high, aes(x=AOO..km2., y = Mean_abundance))+
  geom_point(size = 3)+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))+
  ylab(" ")+
  xlab(" ")

high_fig <- ggMarginal(rs_ab_fig_high, type = "density", fill = "#440154", 
                       alpha = 0.5, color = "#440154")

#combine
all <-patchwork::wrap_elements(low_fig) + 
  patchwork::wrap_elements(mid_fig) + 
  patchwork::wrap_elements(high_fig)+
  plot_layout(guides = 'collect', axes = "collect")+ 
  plot_annotation(tag_levels = c('A'), tag_suffix = ')')&
  theme(plot.tag = element_text(face = 'bold'))+
  theme(plot.tag = element_text(size = 18))

plot(all) 

#all points on one fig-----
fig_all <- ggplot(all_df, aes(x=log(AOO..km2.), y = log(Mean_abundance),
                              color = Site))+
  geom_point()+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))+
  ylab("Log mean local abundance")+
  xlab("Log area of occupancy (km2)")

ggMarginal(fig_all, type = "density", 
                      groupFill = TRUE, alpha = 0.2, groupColor = TRUE)

#faceted figure------
rs_ab_fig <- ggplot(all_df, aes(x=log(AOO..km2.), y = log(Mean_abundance)))+
  geom_point()+
  facet_grid(~Site)+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))+
  xlab("Log area of occupancy (km2)")+
  ylab("Log mean local abundance")+
  geom_density(aes(x=log(AOO..km2.), y=(-2+(after_stat(scaled)))), 
             position="identity")+
  geom_density(aes(y=log(Mean_abundance), x=(3+(after_stat(scaled)))), 
               position="identity")
  
plot(rs_ab_fig)


#extra bits that worked
stat_density(aes(x=log(AOO..km2.), y=(-2+(after_stat(scaled)))), 
             position="identity", geom="line")+
  stat_density(aes(y=log(Mean_abundance), x=(2+(after_stat(scaled)))), 
               position="identity", geom="line")

#add to overall fi---------
threepanel_fig <- all / fig_ab / fig_rs  +  
  plot_annotation(tag_levels = c('A'), tag_suffix = ')')+plot_layout(guides = 'collect')
  

plot(threepanel_fig)


