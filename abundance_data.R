library(dplyr)
library(tidyverse)

#looking at RMBL abundance data to sort 

#bring in data
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research")

abundance <- read.csv("RMBL_abundance_2021-2022.csv")
abundance <- as.data.frame(abundance)

#split by year
abundance2021 <- abundance[ which(abundance$Year==2021), ]
abundance2022 <- abundance[ which(abundance$Year==2022), ]

#check distributions of abundance values 
hist2021 <- ggplot(abundance2021, aes(x=Total)) + 
  geom_histogram(binwidth = 50) + ggtitle(2021)
hist2021

hist2022 <- ggplot(abundance2022, aes(x=Total)) + 
  geom_histogram(binwidth = 50) + ggtitle(2022)
hist2022

#average abundance values across both years
abundance_avg <- abundance %>%
  group_by(Site, Species) %>%
  summarize(mean_abundance = mean(Total))

write.csv(abundance_avg, file="phylogenetics-and-ranges/abundance_avg.csv")

