#tutorial from Guisan textboook - Example 1 starting on page 357

if(!require(rgbif)){}
library(rgbif)

install.packages("biomod2")
library(biomod2)
library(ggplot2)
library(gridExtra)

install.packages("raster")
install.packages("rgdal")
install.packages("sf")
install.packages("terra")
install.packages("sp")
install.packages("rgeos")
library(raster)
library(rgdal)
library(sf)
library(terra)
library(sp)
library(rgeos)
library(dplyr)


install.packages("usethis")
usethis::edit_r_environ()

#other method from RGBIF vignette (https://github.com/ropensci/rgbif/blob/dfd376e55c5a4f25f7ed5b74427b1b559ce8bc7f/vignettes/getting_occurrence_data.Rmd)

#get taxon key 
name_backbone("Agoseris aurantiaca")
key1<-name_backbone(name='Agoserisaurantiaca',rank='species')$usageKey

#download occurrence data from GBIF directly--------------------------------

occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"),
  pred_gte("year",1990),
  pred("taxonKey", 3092928),
  format = "SIMPLE_CSV",
  user="leah.veldhuisen", 
  pwd="Columbia2305", 
  email="leah.veldhuisen@gmail.com"
)

d <- occ_download_get('0008582-231002084531237') %>%
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
AgoAur_occ <- df1[df1$species == "Agoseris.aurantiaca",]

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
shapeNA <- unzip(zipfile = "politcalboundaries_shapefile.zip", exdir = "/Users/leahvedlhuisen/Downloads", overwrite = T)
list.files("bound_p", recursive = T)

mask_of_NA <- shapefile("bound_p/boundaries_p_2021_v3.shp")
bioclim_NA <- mask(bioclim_world, mask_of_NA[ mask_of_NA,]) 
bioclim_NA <- mask(bioclim_world, mask_of_NA) #alternate option if top line doesnt work 
bioclim_NA <- crop(bioclim_NA, mask_of_NA)

load("~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/distribution modeling/mask_of_ NA.rda")


###mask_of_NA code is slow to run, save object################# 

save(mask_of_NA, file="~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/distribution modeling/mask_of_ NA.rda")

writeOGR(mask_of_NA, dsn="~/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/distribution modeling/mask_of_NA.shp", layer="mask_of_NA", driver="ESRI Shapefile")

#get occurrence points to align with env data
points_agoseris <- data.frame(AgoAur_occ[1:290,c("decimalLongitude", "decimalLatitude")])
AgoAur_cell_id <- cellFromXY(subset(bioclim_NA,1), points_agoseris)

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
                                     %in% AgoAur_cell_id, levels = c("FALSE","TRUE"),
                                     labels=c("background","AgoAur")), col = c("red","blue"),
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
                                     %in% AgoAur_cell_id, levels = c("FALSE","TRUE"),
                                     labels=c("background","AgoAur")), col = c("red","blue"),
        csta = 0, cellipse = 2, cpoint = .3, pch = 16)
mtext("(a)",side = 3, line = 3, adj = 0)
s.corcircle(pca_NA_soil$co, clabel = .5)
mtext("(b)",side = 3, line = 3, adj = 0)

#try other ways to assess collinearity in soil and bioclim vars----------------
install.packages("usdm")
library(usdm)
vif(soil_NA_df)

vifcor(soil_NA_df, th =0.7)
vifcor(bioclim_NA_df, th =0.7)

#VIF with combo dataset with bioclim vars and soil vars------------------
all_NA_df <- na.omit(as.data.frame(combo_env_vars))
vif(all_NA_df)
vifcor(all_NA_df, th = 0.7)

#use chosen vars to run model -- cutoff VIF at 4-----------------
env_vars_NA_sub <- stack(subset(combo_env_vars, c("bio_4","bio_10","bio_12","Unified_NA_Soil_Map_Topsoil_Gravel_Content","Unified_NA_Soil_Map_Topsoil_Silt_Fraction","Unified_NA_Soil_Map_Topsoil_Cation_Exchange_Capacity", "Unified_NA_Soil_Map_Topsoil_Organic_Carbon")))

#run models---------------------
library(biomod2)
install.packages("devtools")
library(devtools)


#switch format of env vars to see if it runs 
env_vars_NA_sub_test <- as(env_vars_NA_sub, "SpatRaster")
expl.var <- droplevels(env_vars_NA_sub_test)

AgoAur_occ_df <- as.data.frame(AgoAur_occ)

resp_var <- rep(1,nrow(AgoAur_occ))

#format and run models--------------
library(data.table)
library(terra)
library(biomod2)

##format data for biomod2---------
setwd("/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/distribution modeling")

AgoAur_occ = fread("AgoAur_occ.csv", data.table = FALSE)
load("expl.var1.rda")


env_vars_NA_sub_test = unwrap(env_vars_NA_sub_test)
env_vars_NA_sub_test$bio_4 = as.numeric(env_vars_NA_sub_test$bio_4)
env_vars_NA_sub_test$bio_10 = as.numeric(env_vars_NA_sub_test$bio_10)
env_vars_NA_sub_test$bio_12 = as.numeric(env_vars_NA_sub_test$bio_12)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Gravel_Content = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Gravel_Content)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Silt_Fraction = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Silt_Fraction)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Cation_Exchange_Capacity = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Cation_Exchange_Capacity)
env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Organic_Carbon = as.numeric(env_vars_NA_sub_test$Unified_NA_Soil_Map_Topsoil_Organic_Carbon)

###formatting data############
AgoAur_data <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(AgoAur_occ)), 
  expl.var = env_vars_NA_sub_test,
  resp.xy = AgoAur_occ[,c('decimalLongitude','decimalLatitude')],
  resp.name = "Agoseris.aurantiaca",
  PA.nb.rep = 3,
  PA.nb.absences = 500, 
  PA.strategy = 'random',
  filter.raster = TRUE
)

AgoAur_occ[,c('decimalLongitude','decimalLatitude')]

AgoAur_data

###make table for PAs#############
PA <- bm_PseudoAbsences( resp.var = rep(1, nrow(AgoAur_occ)), 
                   expl.var = env_vars_NA_sub_test, 
                   nb.rep=3, 
                   strategy="random", 
                   nb.absences=500, 
                   sre.quant=0, 
                   dist.min=0, 
                   dist.max=NULL, 
                   user.table=NULL )

###optimizing data##############
install.packages("biomod2", dependencies = TRUE)
library(devtools)


AgoAur_opt <- BIOMOD_ModelingOptions(
  GLM = list(type ='quadratic', interaction.level=1),
  GBM = list(n.trees=1000), 
  GAM = list(algo="GAM_mgcv"))


AgoAur_opt2 <- bm_ModelingOptions(data.type = 'binary',  #######this is the code for the 4.2-5 newer version of biomod2 
                  bm.format = AgoAur_data, 
                  models = c('GLM','GBM','GAM','RF'),
                  strategy = 'default')

###running models##############
AgoAur_models <- BIOMOD_Modeling(AgoAur_data, 
                                 models = c("GLM","GBM","RF","GAM"),
  bm.options= AgoAur_opt, 
  CV.nb.rep = 4, 
  CV.perc =0.8,
  var.import = 3, 
  CV.do.full.models = F, 
  modeling.id = "ex2", 
  do.progress = T, 
  metric.eval=c('TSS','ROC'))

###get model scores#######
AgoAur_models_scores <- get_evaluations(AgoAur_models)
dim(AgoAur_models_scores)
dimnames(AgoAur_models_scores)

bm_PlotEvalMean(bm.out = AgoAur_models)

###testing different variables importance - did not do########
AgoAur_models_var_import <- get_variables_importance(AgoAur_models)
apply(AgoAur_models_var_import, c(1,2), mean)

bm_VariablesImportance( bm.model = AgoAur_models, 
                        expl.var = env_vars_NA_sub_test, 
                        method="full_rand", 
                        nb.rep=1, 
                        do.progress=TRUE)

#ensemble modeling-----------------
AgoAur_ensemble_model <- BIOMOD_EnsembleModeling(AgoAur_models, 
                                                 em.by = 'all',
                                                 metric.select = 'TSS',
                                                 metric.eval=c("KAPPA","TSS","ROC"), 
                                                 metric.select.thresh = 0.8,
                                                 em.algo = c('EMcv', 'EMca', 'EMwmean'),
                                                 var.import = 0)
AA_ensemble_model_scores <- get_evaluations(AgoAur_ensemble_model)
AA_ensemble_model_scores

#make maps using ensemble modeling------------------------------
AgoAur_models_proj_current <- BIOMOD_Projection(bm.mod = AgoAur_models, 
                                                new.env = env_vars_NA_sub_test,
                                                proj.name = 'current',
                                                metric.binary = 'TSS',
                                                output.format=".img",
                                                do.stack = FALSE,
                                                build.clamping.mask = FALSE)

AgoAur_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(bm.em = AgoAur_ensemble_model,
                                                                  bm.proj =AgoAur_models_proj_current,
                                                                  proj.name = 'test',
                                                                  metric.binary = "TSS",
                                                                  output.format = ".img",
                                                                  do.stack = FALSE,
                                                                  models.chosen = "all"
                                                                  )
plot(AgoAur_models_proj_current)
writeRaster(AgoAur_ensemble_models_proj_current, filename = "/Users/leahvedlhuisen/Library/CloudStorage/OneDrive-UniversityofArizona/Arizona PhD/Research/Chapter 2/distribution modeling/agoaur.img", format = "HFA" )

#play with plotting--------
install.packages("arulesViz")
library(arulesViz)

plot(
  AgoAur_ensemble_models_proj_current,
  coord = AgoAur_occ_df,
  plot.output = 'facet',
  do.plot = TRUE,
  std = TRUE,
  scales = 'fixed',
  size,
  maxcell = 5e+05
)

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

##generate gps point for RMBL###########
gps_point <- data.frame(x = -106.5843, y = 38.5350)
gps_point <- as.data.frame(gps_point)
extent(env_layer4)

###test point for Las Vegas#######
test_point <- data.frame(x = -75.1652, y=39.9526)

testxy <- test_point[,c(1,2)]

test_spdf <- SpatialPointsDataFrame(coords = testxy, data = test_point,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
test123 <- test_point[,c(1,2)]

##make point in a coordinate system#####
xy <- gps_point[,c(1,2)]

spdf <- SpatialPointsDataFrame(coords = xy, data = gps_point,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

plot(spdf)

##extract environmental variable values for this point#####
env_values <- extract(env_vars_stacked, spdf, method = 'simple')

test_env_values <- extract(env_vars_stacked, test_spdf, method = 'bilinear')

##use model to generate suitability value at this point 
models.needed <- get_kept_models(AgoAur_ensemble_model)
formal_pred <- BIOMOD_Projection(bm.mod = AgoAur_models,
                                 new.env = env_values,
                                 proj.name = "current",
                                 new.env.xy = xy,
                                 models.chosen = 'all',
                                 compress = TRUE,
                                 build.clamping.mask = FALSE,
                                 do.stack = TRUE,
                                 nb.cpu = 1)
prediction_table <- get_predictions(formal_pred, full.name = models.needed)
plot(prediction_table)
mean(prediction_table$pred)
median(prediction_table$pred)

###test for Walla Walla, WA give avg of 254 and for Philadelphia give 54