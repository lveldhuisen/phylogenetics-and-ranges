library(rgbif)
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rgdal)
library(sf)
library(terra)
library(sp)
library(rgeos)
library(dplyr)


install.packages("usethis")
usethis::edit_r_environ()

#get taxon key 
backbone("Agoseris_glauca")
key1<-name_backbone(name='Vicia_americana',rank='species')$usageKey

#download occurrence data from GBIF directly--------------------------------

occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"),
  pred_gte("year",1990),
  pred("taxonKey",2975204),
  format = "SIMPLE_CSV",
  user="leah.veldhuisen", 
  pwd="Columbia2305", 
  email="leah.veldhuisen@gmail.com"
)

d <- occ_download_get('0020092-231120084113126') %>%
  occ_download_import()
library(dplyr)

##clean up data list to only have species name, lat and long##############
df1 <- d %>%
  select(species, decimalLatitude, decimalLongitude,coordinateUncertaintyInMeters)
df1 <- df1[df1$coordinateUncertaintyInMeters < 1000, ]
df1 <- na.omit(df1)
df1 <- df1 %>%
  select(species, decimalLatitude, decimalLongitude)
df1$species <- sub(" ",".", df1$species)

#download WorldClim environmental data###################
dir.create("WorldClim_data", showWarnings = F)
download.file(url = 
                "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_esri.zip",
              destfile = "WorldClim_data/current_bioclim_10min.zip",method = "auto")

##unzip climatic files
###bioclim variables-------------
unzip(zipfile = "WorldClim_data/current_bioclim_10min.zip", exdir = "WorldClim_data/current", overwrite = T)
list.files("WorldClim_data/current/bio/")

#soil variables---------------
setwd("/Users/leahvedlhuisen/Downloads")
unzip(zipfile = "NACP_MsTMIP_Unified_NA_SoilMap_1242.zip", exdir = "/Users/leahvedlhuisen/Downloads", overwrite = T)
list.files("NACP_MsTMIP_Unified_NA_SoilMap_1242/data")

#extract species occurrences
species_occ <- df1

###stack bioclim variables into one file#################### 
install.packages("raster")
library(raster)
bioclim_world <- stack(list.files("WorldClim_data/current/bio/", 
                                  pattern = "bio_", full.names = T), RAT = FALSE)


#code to stack soil variables with bioclim variables#####################
soil_vars <- stack(list.files("NACP_MsTMIP_Unified_NA_SoilMap_1242/data", 
                              pattern = "Unified_NA_Soil_Map_", full.names = T), RAT = FALSE)

##resample soil files to match bioclim##############
soil1 <- resample(soil_vars, bioclim_world, method="bilinear")

###check resolution and extent of bioclim vars and soil vars###################
bioclim_NA
soil1

#combine bioclim vars with resampled soil vars#############
combo_env_vars <- stack(bioclim_world,soil1)
brick(combo_env_vars)

#north america shape file---------------------
setwd("/Users/leahvedlhuisen/Downloads")

load("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/distribution_modeling/mask_of_ NA.rda")
bioclim_NA <- mask(bioclim_world, mask_of_NA[ mask_of_NA,]) 
bioclim_NA <- mask(bioclim_world, mask_of_NA) #alternate option if top line doesnt work 
bioclim_NA <- crop(bioclim_NA, mask_of_NA)

#get occurrence points to align with env data
points_species <- data.frame(species_occ[,c("decimalLongitude", "decimalLatitude")])
species_cell_id <- cellFromXY(subset(bioclim_NA,1), points_species)

#pca to test correlation between env data vars--------------------- 
install.packages("ade4")
library(ade4)

##bioclim vars#############
bioclim_NA_df <- na.omit(as.data.frame(bioclim_NA))
head(bioclim_NA_df)

pca_NA <- dudi.pca(bioclim_NA_df, scannf = F, nf = 2)

plot(pca_NA$li[,1:2])

###check PCA outliers#########
###tail of distributions$###########
sort(pca_NA$li[,1])[1:10]

###ID points to remove############
to_remove <- which((pca_NA$li[, 1] < -10))

###remove points and redo PCA############
if(length(to_remove)){
  bioclim_NA_df <- bioclim_NA_df[- to_remove,]
  pca_NA<- dudi.pca(bioclim_NA_df,scannf = F, nf=2)
}

###pca with removed outliers to test collinearity###############
par(mfrow=c(1,2))
s.class(pca_NA$li[,1:2],fac = factor(rownames(bioclim_NA_df)
                                     %in% species_cell_id, levels = c("FALSE","TRUE"),
                                     labels=c("background","species")), col = c("red","blue"),
        csta = 0, cellipse = 2, cpoint = .3, pch = 16)
mtext("(a)",side = 3, line = 3, adj = 0)
s.corcircle(pca_NA$co, clabel = .5)
mtext("(b)",side = 3, line = 3, adj = 0)


#soil variables PCA---------------------
soil_NA_df <- na.omit(as.data.frame(soil1))
head(soil_NA_df)

pca_NA_soil <- dudi.pca(soil_NA_df, scannf = F, nf = 2)

plot(pca_NA_soil$li[,1:2])

###check PCA outliers###########
###tail of distributions##############
sort(pca_NA_soil$li[,1])[1:10]

###ID points to remove#################
to_remove_soil <- which((pca_NA_soil$li[, 1] < -10))

###remove points and redo PCA#####################
if(length(to_remove_soil)){
  soil_NA_df <- soil_NA_df[- to_remove_soil,]
  pca_NA_soil<- dudi.pca(soil_NA_df,scannf = F, nf=2)
}

###pca with removed outliers to test collinearity################
par(mfrow=c(1,2))
s.class(pca_NA_soil$li[,1:2],fac = factor(rownames(soil_NA_df)
                                          %in% species_cell_id, levels = c("FALSE","TRUE"),
                                          labels=c("background","species")), col = c("red","blue"),
        csta = 0, cellipse = 2, cpoint = .3, pch = 16)
mtext("(a)",side = 3, line = 3, adj = 0)
s.corcircle(pca_NA_soil$co, clabel = .5)
mtext("(b)",side = 3, line = 3, adj = 0)

#use chosen vars to run model -- cutoff VIF at 4-----------------
env_vars_NA_sub <- stack(subset(combo_env_vars, c("bio_4","bio_10","bio_12","Unified_NA_Soil_Map_Topsoil_Gravel_Content","Unified_NA_Soil_Map_Topsoil_Silt_Fraction","Unified_NA_Soil_Map_Topsoil_Cation_Exchange_Capacity", "Unified_NA_Soil_Map_Topsoil_Organic_Carbon")))

#run models---------------------
library(biomod2)
install.packages("devtools")
library(devtools)


#switch format of env vars to see if it runs 
env_vars_NA_sub_test <- as(env_vars_NA_sub, "SpatRaster")

expl.var <- droplevels(env_vars_NA_sub_test)

species_occ_df <- as.data.frame(species_occ)

resp_var <- rep(1,nrow(species_occ))

#format and run models--------------
library(data.table)
library(terra)
library(biomod2)

##format data for biomod2---------
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/distribution modeling")


env_vars_NA_sub_test = unwrap(env_vars_NA_sub_test)
env_vars_NA_sub_test$bio_4 = as.numeric(env_vars_NA_sub_test$bio_4)
env_vars_NA_sub_test$bio_10 = as.numeric(env_vars_NA_sub_test$bio_10)

env_vars_NA_sub_test$bio_12 = as.numeric(env_vars_NA_sub_test$bio_12)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Gravel_Content = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Gravel_Content)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Silt_Fraction = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Silt_Fraction)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Cation_Exchange_Capacity = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Cation_Exchange_Capacity)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Organic_Carbon = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Organic_Carbon)

###formatting data############
species_data <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(species_occ)), 
  expl.var = env_vars_NA_sub_test,
  resp.xy = species_occ[,c('decimalLongitude','decimalLatitude')],
  resp.name = "Vicia_americana",
  PA.nb.rep = 3,
  PA.nb.absences = 500, 
  PA.strategy = 'random',
  filter.raster = TRUE
)

species_occ[,c('decimalLongitude','decimalLatitude')]

species_data

###optimizing data##############
install.packages("biomod2", dependencies = TRUE)
library(devtools)


species_opt <- BIOMOD_ModelingOptions(
  GLM = list(type ='quadratic', interaction.level=1),
  GBM = list(n.trees=1000), 
  GAM = list(algo="GAM_mgcv"))


AgoAur_opt2 <- bm_ModelingOptions(data.type = 'binary',  #######this is the code for the 4.2-5 newer version of biomod2 
                                  bm.format = AgoAur_data, 
                                  models = c('GLM','GBM','GAM','RF'),
                                  strategy = 'default')

###running models##############
species_models <- BIOMOD_Modeling(species_data, 
                                 models = c("GBM","RF","GAM","GLM"),
                                 bm.options= species_opt, 
                                 CV.nb.rep = 4, 
                                 CV.perc =0.8,
                                 var.import = 3, 
                                 CV.do.full.models = F, 
                                 modeling.id = "Vicia_americana", 
                                 do.progress = T, 
                                 metric.eval=c('TSS','ROC'))

###get model scores#######
species_models_scores <- get_evaluations(species_models)
dim(species_models_scores)
dimnames(species_models_scores)

bm_PlotEvalMean(bm.out = species_models)


#ensemble modeling-----------------
species_ensemble_model <- BIOMOD_EnsembleModeling(species_models, 
                                                 em.by = 'all',
                                                 metric.select = 'TSS',
                                                 metric.eval=c("KAPPA","TSS","ROC"), 
                                                 metric.select.thresh = 0.8,
                                                 em.algo = c('EMcv', 'EMca', 'EMwmean'),
                                                 var.import = 0)
species_ensemble_model_scores <- get_evaluations(species_ensemble_model)
species_ensemble_model_scores
bm_PlotEvalMean(bm.out = species_ensemble_model)

#make maps using ensemble modeling------------------------------
species_models_proj_current <- BIOMOD_Projection(bm.mod = species_models, 
                                                new.env = env_vars_NA_sub_test,
                                                proj.name = 'Vicia_americana',
                                                metric.binary = 'TSS',
                                                output.format=".img",
                                                do.stack = FALSE,
                                                build.clamping.mask = FALSE,
                                                output.format = '.img',
                                                type = "RasterLayer")

species_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(bm.em = species_ensemble_model,
                                                                  bm.proj =species_models_proj_current,
                                                                  proj.name = 'Vicia_americana',
                                                                  metric.binary = "TSS",
                                                                  output.format = ".img",
                                                                  do.stack = FALSE,
                                                                  models.chosen = "all"
)


saveRDS(species_ensemble_models_proj_current, file = "phylogenetics&ranges/H.hoopesii.rds")
plot(species_models_proj_current)
plot(species_ensemble_models_proj_current)

#test ways of calculating range size-----------------
install.packages("BiodiversityR")
install.packages(pkgs=c("BiodiversityR", "vegan","Rcmdr", "MASS", "mgcv","cluster", "RODBC", "rpart", "effects", "multcomp","ellipse", "maptree", "sp", "splancs", "spatial","akima", "nnet", "dismo", "raster", "rgdal","bootstrap", "PresenceAbsence","maxlike", "gbm", "randomForest", "gam", "earth", "mda","kernlab", "e1071", "glmnet", "sem", "rgl", "relimp","lmtest", "leaps", "Hmisc", "colorspace", "aplpack","abind", "XLConnect", "car", "markdown", "knitr","geosphere", "maptools", "rgeos", "ENMeval", "red"),dependencies=c("Depends", "Imports"))
install.packages("tcltk")
install.packages("red")
library(BiodiversityR)
library(tcltk2)
library(red)
library(raster)


#use BiodiversityR to calculate AOO and EOO 
ensemble.red(test)
aoo(points_species)
eoo(points_species)

#get habitat suitability value at specific point---------------------------- 
install.packages("dismo")
library(dismo)
library(raster)
library(rgdal)
setwd("/Users/leahvedlhuisen/Downloads/WorldClim_data/current/bio")

##make individual layers#######
env_layer1 <- raster(env_vars_NA_sub_test$bio_4)
env_layer2 <- raster(env_vars_NA_sub_test$bio_10)
env_layer3 <- raster(env_vars_NA_sub_test$bio_12)
env_layer4 <- raster(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Gravel_Content)
env_layer5 <- raster(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Silt_Fraction)
env_layer6 <- raster(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Cation_Exchange_Capacity)
env_layer7 <- raster(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Organic_Carbon)
##make stacked raster object with all 7 layers########
env_vars_stacked <- stack(env_layer1, env_layer2,env_layer3, env_layer4, env_layer5,env_layer6,env_layer7)

##make points for all RMBL sites#####
gps_point_road <- data.frame(x = -106.5843, y = 38.5350)
gps_point_road <- as.data.frame(gps_point_road)

gps_point_pfeiler <- data.frame(x = -107.0153, y = 38.5738)
gps_point_pfeiler <- as.data.frame(gps_point_pfeiler)

gps_point_PBM <- data.frame(x = -107.0153, y = 38.5810)
gps_point_PBM <- as.data.frame(gps_point_PBM)

##make point in a coordinate system#####
xy_road <- gps_point_road[,c(1,2)]
spdf_road <- SpatialPointsDataFrame(coords = xy_road, data = gps_point_road,
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
xy_pfeiler <- gps_point_pfeiler[,c(1,2)]
spdf_pfeiler <- SpatialPointsDataFrame(coords = xy_pfeiler, data = gps_point_pfeiler,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

xy_PBM <- gps_point_PBM[,c(1,2)]
spdf_PBM <- SpatialPointsDataFrame(coords = xy_PBM, data = gps_point_PBM,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#function to extract env var values for each point###
env_values_function <- function(point){
  extract(env_vars_stacked, point, method = 'simple')
}

env_valuesPBM <- env_values_function(point = spdf_PBM)
env_valuespfeiler <- env_values_function(point = spdf_pfeiler)
env_valuesroad <- env_values_function(point = spdf_road)

##use model to generate suitability value at this point 
models.needed <- get_kept_models(species_ensemble_model)
formal_pred <- BIOMOD_Projection(bm.mod = species_models,
                                 new.env = env_valuesroad,
                                 proj.name = "current",
                                 new.env.xy = xy_road,
                                 models.chosen = 'all',
                                 compress = TRUE,
                                 build.clamping.mask = FALSE,
                                 do.stack = TRUE,
                                 nb.cpu = 1)
prediction_table <- get_predictions(formal_pred)
mean(prediction_table$pred)
median(prediction_table$pred)
sd(prediction_table$pred)

#range size prediction
install.packages("ecospat")
library(ecospat)
library(terra) 
library(dismo)

#convert biomod2 model to raster 
myProjDF <- get_predictions(species_models_proj_current, as.data.frame=TRUE)

as(myProjDF, "SpatRaster")

range_size_results <- ecospat.rangesize (bin.map = myProjDF,
                                         ocp = TRUE,
                                         eoo.around.model = TRUE,
                                         EOO = TRUE,
                                         Model.within.eoo = TRUE,
                                         AOO = TRUE,
                                         resol = c(2000, 2000),
                                         d = sqrt((2000 * 2)/pi),
                                         lonlat = TRUE,
                                         return.obj = TRUE)

test <- if(require(alphahull,quietly=TRUE)){ 
  rangesize2<-ecospat.rangesize(c(myProjDF), 
                                xy=species_occ_df, 
                                AOO.circles=TRUE, 
                                lonlat=FALSE)
}



  
