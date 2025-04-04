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

write.csv(abundance_avg, file="/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/abundance_avg.csv")

#test for differences in abundance based on range size 
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2")

results <- read.csv("results_all.csv")

#remove last 2 empty rows 
results <- results[-c(90,91),]

results_certain <- results[-c(69:90),]

#Kruskal wallis test 
##only using abundance values from Veronica ---- p value is 0.8365
kruskal.test(Mean_abundance ~ Range_size, data = results_certain)
ggplot(results_certain, aes(x=Range_size, y=Mean_abundance)) + geom_boxplot()

##including added abundance values from from my flowering data---- p value is 0.5635
kruskal.test(Mean_abundance ~ Range_size, data = results)
ggplot(results, aes(x=Range_size, y=Mean_abundance)) + geom_boxplot()
