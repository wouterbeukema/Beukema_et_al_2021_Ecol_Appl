########################################################################
#################Bsal landscape epidemiology February 2020##############
###############Data checking, Species Distribution Modelling############
########################################################################

library(graphics) # for graphic::pairs
library(blockCV) # for blockCV::spatialBlock, to make 47-fold DataSplitTable.
library(ecospat)
library(biomod2)
library(raster) # for raster::extract
library(rgdal) # to work with shapefiles
library(sf) # for st_crs function
library(ggplot2)

# Set working directory
setwd("")

########################Check collinearity############################

# First sample environmental parameters at Bsal pres-abs sites;
# Upload pres-abs data;
Bsal <- na.exclude(read.delim("Bsal.txt", h = T, sep = "\t", dec = ","))
# Drop pres-abs column, otherwise raster::extract doesn't work;
Bsal <- Bsal[,1:2]
# Upload parameters as a stack, extract environmental data;
# The '$' prevents img.xml etc. from uploading;
img <- list.files(path = "./rasters", pattern ='\\.img$', full.names = TRUE) 
# Create the actual stack; 
env.stack <- stack(img)
# Extract;
Bsalsample <- as.data.frame(extract(env.stack, Bsal))
Bsalsample["x"] <- Bsal$x # Add x
Bsalsample["y"] <- Bsal$y # Add y
# First upload Bsal.txt again
Bsal <- na.exclude(read.delim("Bsal.txt", h = T, sep = "\t", dec = ","))
Bsalsample["Bsal"] <- Bsal$Bsal # Add pres-abs column
write.csv(Bsalsample, file="Bsalsample.csv") # export if needed

# Then check collinearity;
# Transform parameters with particularly large values;
Bsalsample$alt <- log10(Bsalsample$alt)
Bsalsample$majorroaddens <- log10(Bsalsample$majorroaddens)
Bsalsample$minorroaddens <- log10(Bsalsample$minorroaddens)
Bsalsample$prep_avg <- log10(Bsalsample$prep_avg)
Bsalsample$prep_sd <- log10(Bsalsample$prep_sd)
Bsalsample$recreation <- log10(Bsalsample$recreation)
Bsalsample$slope <- log10(Bsalsample$slope)
Bsalsample$soil_avg <- log10(Bsalsample$soil_avg)
Bsalsample$soil_sd <- log10(Bsalsample$soil_sd)
Bsalsample$tmax_avg <- log10(Bsalsample$tmax_avg)
Bsalsample$traildens <- log10(Bsalsample$traildens)
Bsalsample$waterways <- log10(Bsalsample$waterways)
# Get rid of Inf values, which appear in major- and minorroaddens;
Bsalsample <- do.call(data.frame, lapply(Bsalsample, function(x) replace(x, is.infinite(x), 0)))
# Sometimes the values become NaN, also get rid of those;
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
Bsalsample[is.nan(Bsalsample)] <- 0
# Pearson's corr;
Bsalcor <- as.data.frame(cor(Bsalsample, use = "all.obs", method = "pearson"))
write.csv(Bsalcor, file = "Bsalcor.csv")

###############################biomod2##################################
#######################Upload environmental data########################

# We initially manually set the data type of the host parameters to 'factor', which is the reason why all parameters are uploaded separately below. Including factor parameters however proved to negatively affect the models, so I changed the code to upload everything in standard format.
altitude <- raster("rasters/alt.img")
#ichalp <- raster("rasters/ichalp.img")
#lisvul <- raster("rasters/lisvul.img")
majorroaddens <- raster("rasters/majorroaddens.img")
minorroaddens <- raster("rasters/minorroaddens.img")
PDSI_avg <- raster("rasters/PDSI_avg.img")
PDSI_sd <- raster("rasters/PDSI_sd.img")
popdens <- raster("rasters/popdens.img")
recreation <- raster("rasters/recreation.img")
richness <- raster("rasters/richness.img")
#salsal <- raster("rasters/salsal.img")
shannon <- raster("rasters/shannoncorine.img")
slope <- raster("rasters/slope.img")
soil_avg <- raster("rasters/soil_avg.img")
soil_sd <- raster("rasters/soil_sd.img")
tmax_avg <- raster("rasters/tmax_avg.img")
tmax_sd <- raster("rasters/tmax_sd.img")
traildens <- raster("rasters/traildens.img")
#tricri <- raster("rasters/tricri.img")
waterways <- raster("rasters/waterways.img")

# Stack parameters selected through prior GLM model selection;
env_glm <- stack(slope, shannon, traildens, PDSI_avg, PDSI_sd, soil_avg,
                 soil_sd, tmax_avg, popdens, altitude, recreation, 
                 minorroaddens, majorroaddens)

# Use full dataset for Random Forest;
env_rf <- stack(altitude, majorroaddens, minorroaddens, PDSI_avg,  
                PDSI_sd, popdens, recreation, richness, shannon, slope, 
                soil_avg, soil_sd, tmax_avg, tmax_sd, traildens, 
                waterways)

##############Create blocks for validation (not used*)##################

# *But, I did use the spatialBlock function below to make the 47-fold DataSplitTable.
# Create blocks to restrict spatial autocorrelation and divide dataset into training and testing folds;
# I first created a 'blocks' shapefile in ArcGIS consisting of a polygon that covers Wallonia-Eifel, and one that covers Essen and other spots east of the Rhine. Attribute table only had one field, 'Blocks'. I combine this below with Bsal occurrence data, which I uploaded from a point shapefile to be sure that the blocks and points had the same coordinate system.
blocks <- readOGR(dsn = "./Raw data", layer = "blocks")
blocks <- as(blocks,"SpatialPolygons")
Bsal.SPDF <- readOGR(dsn= "./R_biomod2", layer = "Bsal18-2-2020")

# Check;
plot(shannon)      
points(Bsal.SPDF[Bsal.SPDF$bsal == "0" ,], col = 'white')
points(Bsal.SPDF[Bsal.SPDF$bsal == "1" ,], col = 'red')
plot(blocks[1,], col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3), add = TRUE)
plot(blocks[2,], col = rgb(red = 0, green = 1, blue = 0, alpha = 0.3), add = TRUE)

# Create blocks for biomod2;
folds <- spatialBlock(speciesData = Bsal.SPDF, # Bsal xy data;
                      blocks = blocks, # User-specified blocks;
                      rasterLayer = shannon, # for visualisation purposes;
                      k = 47, # Nr of folds. Cannot be more than nr of blocks
                      showBlocks = TRUE, # plot blocks on tmax_avg
                      biomod2Format = TRUE) # create DataSplitTable
# In case blocks aren't plotted;                 
plot(folds)
# Check DataSplitTable. This should contain a run containing only records in block 1, and a run containing only block-2 records;
folds$biomodTable
DataSplitTable <- folds$biomodTable
# Optional; add a run with all records for ensemble model;
DataSplitTable <- cbind(RUN3 = TRUE, DataSplitTable)

#######################Spatial autocorrelation##########################

Bsal <- na.exclude(read.delim("Bsal.txt", h = T, sep = "\t", dec = ","))
Bsalpres <- Bsal[!(Bsal$Bsal==0),] # get rid of absences
Bsalxy <- Bsalpres[,-3] # drop presence-absence
# Upload prediction parameters obtained through GLM model selection;
env <- env_glm # or env_rf
Bsal_env <- extract(env, Bsalxy) # extract environment at xy
Bsal_env_df <- as.data.frame(Bsal_env)
Bsaly <- Bsalpres[,2]
Bsalx <- Bsalpres[,1]
Bsal_autocor <- cbind(Bsal_env_df, Bsalx, Bsaly) # add x and y columns
# Remove parameters that we excluded based on collinearity;
# Bsal_sample <- Bsal_sample[, -c(12:13, 23:24, 28:29)]
# Remove species;
#Bsal_sample <- Bsal_sample[, -c(5:6, 14, 22)]

# Draw a plot with distance vs. the Mantel r value. Black circles indicate that the values are significantly different from zero. White circles indicate non significant autocorrelation. The selected distance is at the first white circle where values are non-significant different from zero;
ecospat.mantel.correlogram (dfvar = Bsal_autocor[c(1:15)], # entire data frame
                            colxy = 14:15, # x and y columns
                            colvar =  1:13, # predictor parameter columns
                            n = 100, #nr points randomly drawn among pres-abs
                            max = 4, # max distance (I guess in dec degree)
                            nclass = 50, # nr. distance classes
                            nperm = 1000) #nr. permutations

#######################Actual biomod2 modelling#########################

# As input, I ended up adding the 47 buffered leave-one-out crossvalidation folds (i.e., the DataSplitTable) to the Bsal pres-abs file, which I named "DataSplitTable.txt". I removed training records within the spatial autocorrelation range (0.4 decimal degrees) manually in ArcGIS per run. Below, I do a complete biomod2 run for each fold as buffering out records results in unequal pres-abs column lenghts, which biomod2 can't handle. As such, in each run below I first remove all NA's (the buffered-out records), after which I extract the pres-abs data and associated x-y data, and finally overwrite the DataSplitTable1 to a matrix which just contains the DataSplit for the concerning run.

DataSplitTable <- read.delim("DataSplitTable.txt", h = T, sep = "\t", dec = ",")

# Random Forest;
# Model 1
# Get rid of rows with NA;
DataSplitTable1 <- DataSplitTable[!is.na(DataSplitTable$RUN1), ]
# Define species name;
myRespName <- 'Bsal'
# Grab column from file by name;
myResp <- as.numeric(DataSplitTable1 [,myRespName])
# Grab x, y;
myRespXY <- DataSplitTable1 [,c("x","y")]
# Overwrite DataSplitTable to just contain fold 1;
DataSplitTable1 <- as.matrix(na.omit(DataSplitTable$RUN1))

# Format biomod2 input data;
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)

# There is no need to specify specific modelling options as the default settings for GLM are ok. The default 'myBiomodOption' needs to be included in the models though;
myBiomodOption <- BIOMOD_ModelingOptions()

# Random Forest in biomod2;
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable1,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal1_RF", # name of this analysis
                           VarImport = 100) 

# Save evaluation metrics;
Bsal_eval_rf <- as.data.frame(get_evaluations(Bsal_rf))

# Save variable importances;
Bsal_var_rf <- as.data.frame(get_variables_importance(Bsal_rf))

# Create ensemble. Input data is 'modeling.output', we include all models (= 1) through 'chosen models = all'. For the model assembly rules ('em.by') I chose 'PA_dataset', which in our case means that ensemble models are evaluated on the provided evaluation data (= DataSplitTable), which seems most logical to me. More complicated options are available, see function details. The eval.metric is automatically the same as in the original models, so TSS. I took the example threshold (from the function details) of 0.7, which means that only models with a TSS higher than 0.7 are included. To evaluate the ensemble that we create here I again chose the TSS (models.eval.meth). The ensemble is then made as a weighted average, in which weights are awarded for each method proportionally to their evaluation scores (prob.mean.weight.decay = 'proportional').
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')

# Project the initial GLM models. See details below. compress = FALSE means that files on the hard drive won't be compressed. build.clamping.mask = FALSE means that after doing projections we do not check where these cover areas with conditions out of their calibrating/training range. You would usually use this when projecting (transferring) to a different time or region, to check where model uncertainty might occur. build.clamping.mask = TRUE gives an 'out of bounds' error.
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)

# Project the ensemble models;
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)

# Check predictions;
plot(Bsal_rf_projensemble)

# Export predictions;
rf_predictions <- get_predictions(Bsal_rf_projensemble)

# Check which layer you want to export; in the example below I exported the second one, which is the weighted ensemble
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal1_rf.img")

#Model 2;
DataSplitTable2 <- DataSplitTable[!is.na(DataSplitTable$RUN2), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable2 [,myRespName])
myRespXY <- DataSplitTable2 [,c("x","y")]
DataSplitTable2 <- as.matrix(na.omit(DataSplitTable$RUN2))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable2,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal2_RF", # name of this analysis
                           VarImport = 100) 
Bsal2_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal2_eval)
Bsal2_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal2_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal2_rf.img")

# Model 3;
DataSplitTable3 <- DataSplitTable[!is.na(DataSplitTable$RUN3), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable3 [,myRespName])
myRespXY <- DataSplitTable3 [,c("x","y")]
DataSplitTable3 <- as.matrix(na.omit(DataSplitTable$RUN3))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable3,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal3_RF", # name of this analysis
                           VarImport = 100) 
Bsal3_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal3_eval)
Bsal3_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal3_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal3_rf.img")

# Model 4;
DataSplitTable4 <- DataSplitTable[!is.na(DataSplitTable$RUN4), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable4 [,myRespName])
myRespXY <- DataSplitTable4 [,c("x","y")]
DataSplitTable4 <- as.matrix(na.omit(DataSplitTable$RUN4))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable4,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal4_RF", # name of this analysis
                           VarImport = 100) 
Bsal4_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal4_eval)
Bsal4_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal4_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal4_rf.img")

#Model 5
DataSplitTable5 <- DataSplitTable[!is.na(DataSplitTable$RUN5), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable5 [,myRespName])
myRespXY <- DataSplitTable5 [,c("x","y")]
DataSplitTable5 <- as.matrix(na.omit(DataSplitTable$RUN5))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable5,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal5_RF", # name of this analysis
                           VarImport = 100) 
Bsal5_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal5_eval)
Bsal5_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal5_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal5_rf.img")

# Model 6;
DataSplitTable6 <- DataSplitTable[!is.na(DataSplitTable$RUN6), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable6 [,myRespName])
myRespXY <- DataSplitTable6 [,c("x","y")]
DataSplitTable6 <- as.matrix(na.omit(DataSplitTable$RUN6))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable6,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal6_RF", # name of this analysis
                           VarImport = 100) 
Bsal6_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal6_eval)
Bsal6_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal6_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal6_rf.img")

# Model 7;
DataSplitTable7 <- DataSplitTable[!is.na(DataSplitTable$RUN7), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable7 [,myRespName])
myRespXY <- DataSplitTable7 [,c("x","y")]
DataSplitTable7 <- as.matrix(na.omit(DataSplitTable$RUN7))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable7,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal7_RF", # name of this analysis
                           VarImport = 100) 
Bsal7_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal7_eval)
Bsal7_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal7_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal7_rf.img")

# Model 8;
DataSplitTable8 <- DataSplitTable[!is.na(DataSplitTable$RUN8), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable8 [,myRespName])
myRespXY <- DataSplitTable8 [,c("x","y")]
DataSplitTable8 <- as.matrix(na.omit(DataSplitTable$RUN8))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable8,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal8_RF", # name of this analysis
                           VarImport = 100) 
Bsal8_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal8_eval)
Bsal8_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal8_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal8_rf.img")

# Model 9;
DataSplitTable9 <- DataSplitTable[!is.na(DataSplitTable$RUN9), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable9 [,myRespName])
myRespXY <- DataSplitTable9 [,c("x","y")]
DataSplitTable9 <- as.matrix(na.omit(DataSplitTable$RUN9))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable9,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal9_RF", # name of this analysis
                           VarImport = 100) 
Bsal9_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal9_eval)
Bsal9_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal9_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal9_rf.img")

# Model 10;
DataSplitTable10 <- DataSplitTable[!is.na(DataSplitTable$RUN10), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable10 [,myRespName])
myRespXY <- DataSplitTable10 [,c("x","y")]
DataSplitTable10 <- as.matrix(na.omit(DataSplitTable$RUN10))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable10,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal10_RF", # name of this analysis
                           VarImport = 100) 
Bsal10_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal10_eval)
Bsal10_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal10_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal10_rf.img")

#Model 11;
DataSplitTable11 <- DataSplitTable[!is.na(DataSplitTable$RUN11), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable11 [,myRespName])
myRespXY <- DataSplitTable11 [,c("x","y")]
DataSplitTable11 <- as.matrix(na.omit(DataSplitTable$RUN11))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable11,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal11_RF", # name of this analysis
                           VarImport = 100) 
Bsal11_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal11_eval)
Bsal11_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal11_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal11_rf.img")

# Model 12;
DataSplitTable12 <- DataSplitTable[!is.na(DataSplitTable$RUN12), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable12 [,myRespName])
myRespXY <- DataSplitTable12 [,c("x","y")]
DataSplitTable12 <- as.matrix(na.omit(DataSplitTable$RUN12))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable12,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal12_RF", # name of this analysis
                           VarImport = 100) 
Bsal12_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal12_eval)
Bsal12_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal12_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal12_rf.img")

# Model 13;
DataSplitTable13 <- DataSplitTable[!is.na(DataSplitTable$RUN13), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable13 [,myRespName])
myRespXY <- DataSplitTable13 [,c("x","y")]
DataSplitTable13 <- as.matrix(na.omit(DataSplitTable$RUN13))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable13,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal13_RF", # name of this analysis
                           VarImport = 100) 
Bsal13_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal13_eval)
Bsal13_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal13_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal13_rf.img")

# Model 14;
DataSplitTable14 <- DataSplitTable[!is.na(DataSplitTable$RUN14), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable14 [,myRespName])
myRespXY <- DataSplitTable14 [,c("x","y")]
DataSplitTable14 <- as.matrix(na.omit(DataSplitTable$RUN14))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable14,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal14_RF", # name of this analysis
                           VarImport = 100) 
Bsal14_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal14_eval)
Bsal14_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal14_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal14_rf.img")

# Model 15;
DataSplitTable15 <- DataSplitTable[!is.na(DataSplitTable$RUN15), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable15 [,myRespName])
myRespXY <- DataSplitTable15 [,c("x","y")]
DataSplitTable15 <- as.matrix(na.omit(DataSplitTable$RUN15))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable15,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal15_RF", # name of this analysis
                           VarImport = 100) 
Bsal15_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal15_eval)
Bsal15_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal15_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal15_rf.img")

# Model 16;
DataSplitTable16 <- DataSplitTable[!is.na(DataSplitTable$RUN16), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable16 [,myRespName])
myRespXY <- DataSplitTable16 [,c("x","y")]
DataSplitTable16 <- as.matrix(na.omit(DataSplitTable$RUN16))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable16,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal16_RF", # name of this analysis
                           VarImport = 100) 
Bsal16_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal16_eval)
Bsal16_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal16_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal16_rf.img")

# Model 17;
DataSplitTable17 <- DataSplitTable[!is.na(DataSplitTable$RUN17), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable17 [,myRespName])
myRespXY <- DataSplitTable17 [,c("x","y")]
DataSplitTable17 <- as.matrix(na.omit(DataSplitTable$RUN17))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable17,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal17_RF", # name of this analysis
                           VarImport = 100) 
Bsal17_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal17_eval)
Bsal17_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal17_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal17_rf.img")

# Model 18;
DataSplitTable18 <- DataSplitTable[!is.na(DataSplitTable$RUN18), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable18 [,myRespName])
myRespXY <- DataSplitTable18 [,c("x","y")]
DataSplitTable18 <- as.matrix(na.omit(DataSplitTable$RUN18))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable18,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal18_RF", # name of this analysis
                           VarImport = 100) 
Bsal18_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal18_eval)
Bsal18_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal18_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal18_rf.img")

# Model 19;
DataSplitTable19 <- DataSplitTable[!is.na(DataSplitTable$RUN19), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable19 [,myRespName])
myRespXY <- DataSplitTable19 [,c("x","y")]
DataSplitTable19 <- as.matrix(na.omit(DataSplitTable$RUN19))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable19,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal19_RF", # name of this analysis
                           VarImport = 100) 
Bsal19_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal19_eval)
Bsal19_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal19_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal19_rf.img")

# Model 20;
DataSplitTable20 <- DataSplitTable[!is.na(DataSplitTable$RUN20), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable20 [,myRespName])
myRespXY <- DataSplitTable20 [,c("x","y")]
DataSplitTable20 <- as.matrix(na.omit(DataSplitTable$RUN20))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable20,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal20_RF", # name of this analysis
                           VarImport = 100) 
Bsal20_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal20_eval)
Bsal20_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal20_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal20_rf.img")

# Model 21;
DataSplitTable21 <- DataSplitTable[!is.na(DataSplitTable$RUN21), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable21 [,myRespName])
myRespXY <- DataSplitTable21 [,c("x","y")]
DataSplitTable21 <- as.matrix(na.omit(DataSplitTable$RUN21))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable21,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal21_RF", # name of this analysis
                           VarImport = 100) 
Bsal21_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal21_eval)
Bsal21_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal21_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal21_rf.img")

# Model 22;
DataSplitTable22 <- DataSplitTable[!is.na(DataSplitTable$RUN22), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable22 [,myRespName])
myRespXY <- DataSplitTable22 [,c("x","y")]
DataSplitTable22 <- as.matrix(na.omit(DataSplitTable$RUN22))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable22,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal22_RF", # name of this analysis
                           VarImport = 100) 
Bsal22_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal22_eval)
Bsal22_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal22_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal22_rf.img")

# Model 23;
DataSplitTable23 <- DataSplitTable[!is.na(DataSplitTable$RUN23), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable23 [,myRespName])
myRespXY <- DataSplitTable23 [,c("x","y")]
DataSplitTable23 <- as.matrix(na.omit(DataSplitTable$RUN23))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable23,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal23_RF", # name of this analysis
                           VarImport = 100) 
Bsal23_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal23_eval)
Bsal23_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal23_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal23_rf.img")

# Model 24;
DataSplitTable24 <- DataSplitTable[!is.na(DataSplitTable$RUN24), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable24 [,myRespName])
myRespXY <- DataSplitTable24 [,c("x","y")]
DataSplitTable24 <- as.matrix(na.omit(DataSplitTable$RUN24))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable24,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal24_RF", # name of this analysis
                           VarImport = 100) 
Bsal24_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal24_eval)
Bsal24_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal24_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal24_rf.img")

# Model 25;
DataSplitTable25 <- DataSplitTable[!is.na(DataSplitTable$RUN25), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable25 [,myRespName])
myRespXY <- DataSplitTable25 [,c("x","y")]
DataSplitTable25 <- as.matrix(na.omit(DataSplitTable$RUN25))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable25,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal25_RF", # name of this analysis
                           VarImport = 100) 
Bsal25_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal25_eval)
Bsal25_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal25_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal25_rf.img")

# Model 26;
DataSplitTable26 <- DataSplitTable[!is.na(DataSplitTable$RUN26), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable26 [,myRespName])
myRespXY <- DataSplitTable26 [,c("x","y")]
DataSplitTable26 <- as.matrix(na.omit(DataSplitTable$RUN26))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable26,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal26_RF", # name of this analysis
                           VarImport = 100) 
Bsal26_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal26_eval)
Bsal26_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal26_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal26_rf.img")

# Model 27;
DataSplitTable27 <- DataSplitTable[!is.na(DataSplitTable$RUN27), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable27 [,myRespName])
myRespXY <- DataSplitTable27 [,c("x","y")]
DataSplitTable27 <- as.matrix(na.omit(DataSplitTable$RUN27))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable27,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal27_RF", # name of this analysis
                           VarImport = 100) 
Bsal27_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal27_eval)
Bsal27_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal27_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal27_rf.img")

# Model 28;
DataSplitTable28 <- DataSplitTable[!is.na(DataSplitTable$RUN28), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable28 [,myRespName])
myRespXY <- DataSplitTable28 [,c("x","y")]
DataSplitTable28 <- as.matrix(na.omit(DataSplitTable$RUN28))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable28,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal28_RF", # name of this analysis
                           VarImport = 100) 
Bsal28_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal28_eval)
Bsal28_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal28_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal28_rf.img")

# Model 29;
DataSplitTable29 <- DataSplitTable[!is.na(DataSplitTable$RUN29), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable29 [,myRespName])
myRespXY <- DataSplitTable29 [,c("x","y")]
DataSplitTable29 <- as.matrix(na.omit(DataSplitTable$RUN29))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable29,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal29_RF", # name of this analysis
                           VarImport = 100) 
Bsal29_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal29_eval)
Bsal29_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal29_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal29_rf.img")

# Model 30;
DataSplitTable30 <- DataSplitTable[!is.na(DataSplitTable$RUN30), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable30 [,myRespName])
myRespXY <- DataSplitTable30 [,c("x","y")]
DataSplitTable30 <- as.matrix(na.omit(DataSplitTable$RUN30))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable30,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal30_RF", # name of this analysis
                           VarImport = 100) 
Bsal30_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal30_eval)
Bsal30_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal30_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal30_rf.img")

# Model 31;
DataSplitTable31 <- DataSplitTable[!is.na(DataSplitTable$RUN31), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable31 [,myRespName])
myRespXY <- DataSplitTable31 [,c("x","y")]
DataSplitTable31 <- as.matrix(na.omit(DataSplitTable$RUN31))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable31,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal31_RF", # name of this analysis
                           VarImport = 100) 
Bsal31_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal31_eval)
Bsal31_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal31_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal31_rf.img")

# Model 32;
DataSplitTable32 <- DataSplitTable[!is.na(DataSplitTable$RUN32), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable32 [,myRespName])
myRespXY <- DataSplitTable32 [,c("x","y")]
DataSplitTable32 <- as.matrix(na.omit(DataSplitTable$RUN32))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable32,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal32_RF", # name of this analysis
                           VarImport = 100) 
Bsal32_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal32_eval)
Bsal32_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal32_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal32_rf.img")

# Model 33;
DataSplitTable33 <- DataSplitTable[!is.na(DataSplitTable$RUN33), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable33 [,myRespName])
myRespXY <- DataSplitTable33 [,c("x","y")]
DataSplitTable33 <- as.matrix(na.omit(DataSplitTable$RUN33))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable33,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal33_RF", # name of this analysis
                           VarImport = 100) 
Bsal33_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal33_eval)
Bsal33_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal33_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal33_rf.img")

# Model 34;
DataSplitTable34 <- DataSplitTable[!is.na(DataSplitTable$RUN34), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable34 [,myRespName])
myRespXY <- DataSplitTable34 [,c("x","y")]
DataSplitTable34 <- as.matrix(na.omit(DataSplitTable$RUN34))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable34,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal34_RF", # name of this analysis
                           VarImport = 100) 
Bsal34_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal34_eval)
Bsal34_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal34_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal34_rf.img")

# Model 35;
DataSplitTable35 <- DataSplitTable[!is.na(DataSplitTable$RUN35), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable35 [,myRespName])
myRespXY <- DataSplitTable35 [,c("x","y")]
DataSplitTable35 <- as.matrix(na.omit(DataSplitTable$RUN35))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable35,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal35_RF", # name of this analysis
                           VarImport = 100) 
Bsal35_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal35_eval)
Bsal35_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal35_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal35_rf.img")

# Model 36;
DataSplitTable36 <- DataSplitTable[!is.na(DataSplitTable$RUN36), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable36 [,myRespName])
myRespXY <- DataSplitTable36 [,c("x","y")]
DataSplitTable36 <- as.matrix(na.omit(DataSplitTable$RUN36))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable36,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal36_RF", # name of this analysis
                           VarImport = 100) 
Bsal36_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal36_eval)
Bsal36_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal36_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal36_rf.img")

# Model 37;
DataSplitTable37 <- DataSplitTable[!is.na(DataSplitTable$RUN37), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable37 [,myRespName])
myRespXY <- DataSplitTable37 [,c("x","y")]
DataSplitTable37 <- as.matrix(na.omit(DataSplitTable$RUN37))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable37,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal37_RF", # name of this analysis
                           VarImport = 100) 
Bsal37_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal37_eval)
Bsal37_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal37_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal37_rf.img")

# Model 38;
DataSplitTable38 <- DataSplitTable[!is.na(DataSplitTable$RUN38), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable38 [,myRespName])
myRespXY <- DataSplitTable38 [,c("x","y")]
DataSplitTable38 <- as.matrix(na.omit(DataSplitTable$RUN38))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable38,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal38_RF", # name of this analysis
                           VarImport = 100) 
Bsal38_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal38_eval)
Bsal38_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal38_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal38_rf.img")

# Model 39;
DataSplitTable39 <- DataSplitTable[!is.na(DataSplitTable$RUN39), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable39 [,myRespName])
myRespXY <- DataSplitTable39 [,c("x","y")]
DataSplitTable39 <- as.matrix(na.omit(DataSplitTable$RUN39))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable39,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal39_RF", # name of this analysis
                           VarImport = 100) 
Bsal39_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal39_eval)
Bsal39_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal39_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal39_rf.img")

# Model 40;
DataSplitTable40 <- DataSplitTable[!is.na(DataSplitTable$RUN40), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable40 [,myRespName])
myRespXY <- DataSplitTable40 [,c("x","y")]
DataSplitTable40 <- as.matrix(na.omit(DataSplitTable$RUN40))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable40,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal40_RF", # name of this analysis
                           VarImport = 100) 
Bsal40_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal40_eval)
Bsal40_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal40_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal40_rf.img")

# Model 41;
DataSplitTable41 <- DataSplitTable[!is.na(DataSplitTable$RUN41), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable41 [,myRespName])
myRespXY <- DataSplitTable41 [,c("x","y")]
DataSplitTable41 <- as.matrix(na.omit(DataSplitTable$RUN41))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable41,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal41_RF", # name of this analysis
                           VarImport = 100) 
Bsal41_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal41_eval)
Bsal41_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal41_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal41_rf.img")

# Model 42;
DataSplitTable42 <- DataSplitTable[!is.na(DataSplitTable$RUN42), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable42 [,myRespName])
myRespXY <- DataSplitTable42 [,c("x","y")]
DataSplitTable42 <- as.matrix(na.omit(DataSplitTable$RUN42))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable42,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal42_RF", # name of this analysis
                           VarImport = 100) 
Bsal42_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal42_eval)
Bsal42_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal42_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal42_rf.img")

# Model 43;
DataSplitTable43 <- DataSplitTable[!is.na(DataSplitTable$RUN43), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable43 [,myRespName])
myRespXY <- DataSplitTable43 [,c("x","y")]
DataSplitTable43 <- as.matrix(na.omit(DataSplitTable$RUN43))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable43,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal43_RF", # name of this analysis
                           VarImport = 100) 
Bsal43_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal43_eval)
Bsal43_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal43_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal43_rf.img")

# Model 44;
DataSplitTable44 <- DataSplitTable[!is.na(DataSplitTable$RUN44), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable44 [,myRespName])
myRespXY <- DataSplitTable44 [,c("x","y")]
DataSplitTable44 <- as.matrix(na.omit(DataSplitTable$RUN44))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable44,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal44_RF", # name of this analysis
                           VarImport = 100) 
Bsal44_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal44_eval)
Bsal44_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal44_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal44_rf.img")

# Model 45;
DataSplitTable45 <- DataSplitTable[!is.na(DataSplitTable$RUN45), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable45 [,myRespName])
myRespXY <- DataSplitTable45 [,c("x","y")]
DataSplitTable45 <- as.matrix(na.omit(DataSplitTable$RUN45))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable45,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal45_RF", # name of this analysis
                           VarImport = 100) 
Bsal45_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal45_eval)
Bsal45_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal45_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal45_rf.img")

# Model 46;
DataSplitTable46 <- DataSplitTable[!is.na(DataSplitTable$RUN46), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable46 [,myRespName])
myRespXY <- DataSplitTable46 [,c("x","y")]
DataSplitTable46 <- as.matrix(na.omit(DataSplitTable$RUN46))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable46,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal46_RF", # name of this analysis
                           VarImport = 100) 
Bsal46_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal46_eval)
Bsal46_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal46_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal46_rf.img")

# Model 47;
DataSplitTable47 <- DataSplitTable[!is.na(DataSplitTable$RUN47), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable47 [,myRespName])
myRespXY <- DataSplitTable47 [,c("x","y")]
DataSplitTable47 <- as.matrix(na.omit(DataSplitTable$RUN47))
myBiomodData_rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = env_rf, # predictors
                                        resp.xy = myRespXY, # xy data
                                        resp.name = myRespName, # species name
                                        PA.nb.rep = 0, # no pseudoabsences
                                        na.rm = TRUE)
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_rf <- BIOMOD_Modeling(myBiomodData_rf, # input
                           models = c('RF'), # algorithm
                           models.options = myBiomodOption, # options
                           DataSplitTable = DataSplitTable47,
                           models.eval.meth = c('TSS'), # eval metric
                           modeling.id = "Bsal47_RF", # name of this analysis
                           VarImport = 100) 
Bsal47_eval <- get_evaluations(Bsal_rf)
Bsal_eval_rf <- rbind(Bsal_eval_rf, Bsal47_eval)
Bsal47_var <- get_variables_importance(Bsal_rf)
Bsal_var_rf <- cbind(Bsal_var_rf, Bsal47_var)
Bsal_rf_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_rf, 
                                       chosen.models = 'all', 
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS'),
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
Bsal_rf_proj <- BIOMOD_Projection(modeling.output = Bsal_rf, # input
                                  new.env = env_rf, # projection background
                                  proj.name = 'rf_proj', # output folder
                                  selected.models = 'all', # all models
                                  binary.meth = 'TSS',
                                  compress = FALSE,
                                  build.clamping.mask = FALSE)
Bsal_rf_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_rf_proj, EM.output = Bsal_rf_ens)
plot(Bsal_rf_projensemble)
rf_predictions <- get_predictions(Bsal_rf_projensemble)
names(rf_predictions)
writeRaster(subset(rf_predictions, 2), filename ="Bsal47_rf.img")

# Stack predictions from working directory;
rf_predictions <- list.files(path = "./R_biomod2/RF output", pattern ='\\.img$', full.names = TRUE) 
rf.stack <- stack(rf_predictions)
# Create weighted mean;
weighted.mean(rf.stack, 
              Bsal_eval_rf$Testing.data.RF.RUN1.AllData, 
              na.rm = FALSE, 
              filename = "rf_weighted.img")
# Save evaluation data;
write.csv(Bsal_eval_rf, file = "Bsal_eval_rf.csv")
# Save variables importance;
write.csv(Bsal_var_rf, file = "Bsal_var_rf.csv")


# GLM;
# Model 1
DataSplitTable1 <- DataSplitTable[!is.na(DataSplitTable$RUN1), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable1 [,myRespName])
myRespXY <- DataSplitTable1 [,c("x","y")]
DataSplitTable1 <- as.matrix(na.omit(DataSplitTable$RUN1))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable1,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal1_GLM", # name of this analysis
                            VarImport = 100) 
Bsal_eval_glm <- as.data.frame(get_evaluations(Bsal_glm))
Bsal_var_glm <- as.data.frame(get_variables_importance(Bsal_glm))
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal1_glm.img")

# Model 2;
DataSplitTable2 <- DataSplitTable[!is.na(DataSplitTable$RUN2), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable2 [,myRespName])
myRespXY <- DataSplitTable2 [,c("x","y")]
DataSplitTable2 <- as.matrix(na.omit(DataSplitTable$RUN2))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable2,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal2_GLM", # name of this analysis
                            VarImport = 100) 
Bsal2_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal2_eval)
Bsal2_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal2_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal2_glm.img")

# Model 3;
DataSplitTable3 <- DataSplitTable[!is.na(DataSplitTable$RUN3), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable3 [,myRespName])
myRespXY <- DataSplitTable3 [,c("x","y")]
DataSplitTable3 <- as.matrix(na.omit(DataSplitTable$RUN3))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable3,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal3_GLM", # name of this analysis
                            VarImport = 100) 
Bsal3_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal3_eval)
Bsal3_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal3_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal3_glm.img")

# Model 4;
DataSplitTable4 <- DataSplitTable[!is.na(DataSplitTable$RUN4), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable4 [,myRespName])
myRespXY <- DataSplitTable4 [,c("x","y")]
DataSplitTable4 <- as.matrix(na.omit(DataSplitTable$RUN4))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable4,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal4_GLM", # name of this analysis
                            VarImport = 100) 
Bsal4_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal4_eval)
Bsal4_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal4_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal4_glm.img")

# Model 5;
DataSplitTable5 <- DataSplitTable[!is.na(DataSplitTable$RUN5), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable5 [,myRespName])
myRespXY <- DataSplitTable5 [,c("x","y")]
DataSplitTable5 <- as.matrix(na.omit(DataSplitTable$RUN5))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable5,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal5_GLM", # name of this analysis
                            VarImport = 100) 
Bsal5_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal5_eval)
Bsal5_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal5_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal5_glm.img")

# Model 6;
DataSplitTable6 <- DataSplitTable[!is.na(DataSplitTable$RUN6), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable6 [,myRespName])
myRespXY <- DataSplitTable6 [,c("x","y")]
DataSplitTable6 <- as.matrix(na.omit(DataSplitTable$RUN6))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable6,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal6_GLM", # name of this analysis
                            VarImport = 100) 
Bsal6_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal6_eval)
Bsal6_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal6_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal6_glm.img")

# Model 7;
DataSplitTable7 <- DataSplitTable[!is.na(DataSplitTable$RUN7), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable7 [,myRespName])
myRespXY <- DataSplitTable7 [,c("x","y")]
DataSplitTable7 <- as.matrix(na.omit(DataSplitTable$RUN7))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable7,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal7_GLM", # name of this analysis
                            VarImport = 100) 
Bsal7_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal7_eval)
Bsal7_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal7_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal7_glm.img")

# Model 8;
DataSplitTable8 <- DataSplitTable[!is.na(DataSplitTable$RUN8), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable8 [,myRespName])
myRespXY <- DataSplitTable8 [,c("x","y")]
DataSplitTable8 <- as.matrix(na.omit(DataSplitTable$RUN8))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable8,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal8_GLM", # name of this analysis
                            VarImport = 100) 
Bsal8_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal8_eval)
Bsal8_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal8_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal8_glm.img")

# Model 9;
DataSplitTable9 <- DataSplitTable[!is.na(DataSplitTable$RUN9), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable9 [,myRespName])
myRespXY <- DataSplitTable9 [,c("x","y")]
DataSplitTable9 <- as.matrix(na.omit(DataSplitTable$RUN9))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable9,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal9_GLM", # name of this analysis
                            VarImport = 100) 
Bsal9_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal9_eval)
Bsal9_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal9_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal9_glm.img")

# Model 10;
DataSplitTable10 <- DataSplitTable[!is.na(DataSplitTable$RUN10), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable10 [,myRespName])
myRespXY <- DataSplitTable10 [,c("x","y")]
DataSplitTable10 <- as.matrix(na.omit(DataSplitTable$RUN10))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable10,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal10_GLM", # name of this analysis
                            VarImport = 100) 
Bsal10_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal10_eval)
Bsal10_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal10_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal10_glm.img")

# Model 11;
DataSplitTable11 <- DataSplitTable[!is.na(DataSplitTable$RUN11), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable11 [,myRespName])
myRespXY <- DataSplitTable11 [,c("x","y")]
DataSplitTable11 <- as.matrix(na.omit(DataSplitTable$RUN11))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable11,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal11_GLM", # name of this analysis
                            VarImport = 100) 
Bsal11_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal11_eval)
Bsal11_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal11_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal11_glm.img")

# Model 12;
DataSplitTable12 <- DataSplitTable[!is.na(DataSplitTable$RUN12), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable12 [,myRespName])
myRespXY <- DataSplitTable12 [,c("x","y")]
DataSplitTable12 <- as.matrix(na.omit(DataSplitTable$RUN12))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable12,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal12_GLM", # name of this analysis
                            VarImport = 100) 
Bsal12_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal12_eval)
Bsal12_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal12_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal12_glm.img")

# Model 13;
DataSplitTable13 <- DataSplitTable[!is.na(DataSplitTable$RUN13), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable13 [,myRespName])
myRespXY <- DataSplitTable13 [,c("x","y")]
DataSplitTable13 <- as.matrix(na.omit(DataSplitTable$RUN13))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable13,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal13_GLM", # name of this analysis
                            VarImport = 100) 
Bsal13_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal13_eval)
Bsal13_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal13_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal13_glm.img")

# Model 14;
DataSplitTable14 <- DataSplitTable[!is.na(DataSplitTable$RUN14), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable14 [,myRespName])
myRespXY <- DataSplitTable14 [,c("x","y")]
DataSplitTable14 <- as.matrix(na.omit(DataSplitTable$RUN14))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable14,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal14_GLM", # name of this analysis
                            VarImport = 100) 
Bsal14_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal14_eval)
Bsal14_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal14_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal14_glm.img")

# Model 15; 
DataSplitTable15 <- DataSplitTable[!is.na(DataSplitTable$RUN15), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable15 [,myRespName])
myRespXY <- DataSplitTable15 [,c("x","y")]
DataSplitTable15 <- as.matrix(na.omit(DataSplitTable$RUN15))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable15,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal15_GLM", # name of this analysis
                            VarImport = 100) 
Bsal15_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal15_eval)
Bsal15_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal15_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal15_glm.img")

# Model 16;
DataSplitTable16 <- DataSplitTable[!is.na(DataSplitTable$RUN16), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable16 [,myRespName])
myRespXY <- DataSplitTable16 [,c("x","y")]
DataSplitTable16 <- as.matrix(na.omit(DataSplitTable$RUN16))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable16,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal16_GLM", # name of this analysis
                            VarImport = 100) 
Bsal16_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal16_eval)
Bsal16_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal16_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal16_glm.img")

# Model 17;
DataSplitTable17 <- DataSplitTable[!is.na(DataSplitTable$RUN17), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable17 [,myRespName])
myRespXY <- DataSplitTable17 [,c("x","y")]
DataSplitTable17 <- as.matrix(na.omit(DataSplitTable$RUN17))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable17,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal17_GLM", # name of this analysis
                            VarImport = 100) 
Bsal17_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal17_eval)
Bsal17_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal17_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal17_glm.img")

# Model 18;
DataSplitTable18 <- DataSplitTable[!is.na(DataSplitTable$RUN18), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable18 [,myRespName])
myRespXY <- DataSplitTable18 [,c("x","y")]
DataSplitTable18 <- as.matrix(na.omit(DataSplitTable$RUN18))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable18,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal18_GLM", # name of this analysis
                            VarImport = 100) 
Bsal18_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal18_eval)
Bsal18_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal18_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal18_glm.img")

# Model 19;
DataSplitTable19 <- DataSplitTable[!is.na(DataSplitTable$RUN19), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable19 [,myRespName])
myRespXY <- DataSplitTable19 [,c("x","y")]
DataSplitTable19 <- as.matrix(na.omit(DataSplitTable$RUN19))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable19,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal19_GLM", # name of this analysis
                            VarImport = 100) 
Bsal19_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal19_eval)
Bsal19_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal19_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal19_glm.img")

# Model 20;
DataSplitTable20 <- DataSplitTable[!is.na(DataSplitTable$RUN20), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable20 [,myRespName])
myRespXY <- DataSplitTable20 [,c("x","y")]
DataSplitTable20 <- as.matrix(na.omit(DataSplitTable$RUN20))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable20,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal20_GLM", # name of this analysis
                            VarImport = 100) 
Bsal20_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal20_eval)
Bsal20_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal20_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal20_glm.img")

# Model 21;
DataSplitTable21 <- DataSplitTable[!is.na(DataSplitTable$RUN21), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable21 [,myRespName])
myRespXY <- DataSplitTable21 [,c("x","y")]
DataSplitTable21 <- as.matrix(na.omit(DataSplitTable$RUN21))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable21,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal21_GLM", # name of this analysis
                            VarImport = 100) 
Bsal21_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal21_eval)
Bsal21_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal21_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal21_glm.img")

# Model 22;
DataSplitTable22 <- DataSplitTable[!is.na(DataSplitTable$RUN22), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable22 [,myRespName])
myRespXY <- DataSplitTable22 [,c("x","y")]
DataSplitTable22 <- as.matrix(na.omit(DataSplitTable$RUN22))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable22,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal22_GLM", # name of this analysis
                            VarImport = 100) 
Bsal22_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal22_eval)
Bsal22_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal22_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal22_glm.img")

# Model 23;
DataSplitTable23 <- DataSplitTable[!is.na(DataSplitTable$RUN23), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable23 [,myRespName])
myRespXY <- DataSplitTable23 [,c("x","y")]
DataSplitTable23 <- as.matrix(na.omit(DataSplitTable$RUN23))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable23,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal23_GLM", # name of this analysis
                            VarImport = 100) 
Bsal23_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal23_eval)
Bsal23_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal23_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal23_glm.img")

# Model 24;
DataSplitTable24 <- DataSplitTable[!is.na(DataSplitTable$RUN24), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable24 [,myRespName])
myRespXY <- DataSplitTable24 [,c("x","y")]
DataSplitTable24 <- as.matrix(na.omit(DataSplitTable$RUN24))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable24,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal24_GLM", # name of this analysis
                            VarImport = 100) 
Bsal24_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal24_eval)
Bsal24_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal24_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal24_glm.img")

# Model 25;
DataSplitTable25 <- DataSplitTable[!is.na(DataSplitTable$RUN25), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable25 [,myRespName])
myRespXY <- DataSplitTable25 [,c("x","y")]
DataSplitTable25 <- as.matrix(na.omit(DataSplitTable$RUN25))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable25,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal25_GLM", # name of this analysis
                            VarImport = 100) 
Bsal25_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal25_eval)
Bsal25_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal25_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal25_glm.img")

# Model 26;
DataSplitTable26 <- DataSplitTable[!is.na(DataSplitTable$RUN26), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable26 [,myRespName])
myRespXY <- DataSplitTable26 [,c("x","y")]
DataSplitTable26 <- as.matrix(na.omit(DataSplitTable$RUN26))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable26,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal26_GLM", # name of this analysis
                            VarImport = 100) 
Bsal26_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal26_eval)
Bsal26_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal26_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal26_glm.img")

# Model 27;
DataSplitTable27 <- DataSplitTable[!is.na(DataSplitTable$RUN27), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable27 [,myRespName])
myRespXY <- DataSplitTable27 [,c("x","y")]
DataSplitTable27 <- as.matrix(na.omit(DataSplitTable$RUN27))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable27,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal27_GLM", # name of this analysis
                            VarImport = 100) 
Bsal27_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal27_eval)
Bsal27_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal27_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal27_glm.img")

# Model 28;
DataSplitTable28 <- DataSplitTable[!is.na(DataSplitTable$RUN28), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable28 [,myRespName])
myRespXY <- DataSplitTable28 [,c("x","y")]
DataSplitTable28 <- as.matrix(na.omit(DataSplitTable$RUN28))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable28,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal28_GLM", # name of this analysis
                            VarImport = 100) 
Bsal28_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal28_eval)
Bsal28_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal28_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal28_glm.img")

# Model 29;
DataSplitTable29 <- DataSplitTable[!is.na(DataSplitTable$RUN29), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable29 [,myRespName])
myRespXY <- DataSplitTable29 [,c("x","y")]
DataSplitTable29 <- as.matrix(na.omit(DataSplitTable$RUN29))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable29,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal29_GLM", # name of this analysis
                            VarImport = 100) 
Bsal29_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal29_eval)
Bsal29_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal29_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal29_glm.img")

# Model 30;
DataSplitTable30 <- DataSplitTable[!is.na(DataSplitTable$RUN30), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable30 [,myRespName])
myRespXY <- DataSplitTable30 [,c("x","y")]
DataSplitTable30 <- as.matrix(na.omit(DataSplitTable$RUN30))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable30,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal30_GLM", # name of this analysis
                            VarImport = 100) 
Bsal30_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal30_eval)
Bsal30_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal30_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal30_glm.img")

# Model 31;
DataSplitTable31 <- DataSplitTable[!is.na(DataSplitTable$RUN31), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable31 [,myRespName])
myRespXY <- DataSplitTable31 [,c("x","y")]
DataSplitTable31 <- as.matrix(na.omit(DataSplitTable$RUN31))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable31,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal31_GLM", # name of this analysis
                            VarImport = 100) 
Bsal31_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal31_eval)
Bsal31_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal31_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal31_glm.img")

# Model 32;
DataSplitTable32 <- DataSplitTable[!is.na(DataSplitTable$RUN32), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable32 [,myRespName])
myRespXY <- DataSplitTable32 [,c("x","y")]
DataSplitTable32 <- as.matrix(na.omit(DataSplitTable$RUN32))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable32,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal32_GLM", # name of this analysis
                            VarImport = 100) 
Bsal32_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal32_eval)
Bsal32_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal32_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal32_glm.img")

# Model 33;
DataSplitTable33 <- DataSplitTable[!is.na(DataSplitTable$RUN33), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable33 [,myRespName])
myRespXY <- DataSplitTable33 [,c("x","y")]
DataSplitTable33 <- as.matrix(na.omit(DataSplitTable$RUN33))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable33,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal33_GLM", # name of this analysis
                            VarImport = 100) 
Bsal33_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal33_eval)
Bsal33_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal33_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal33_glm.img")

# Model 34;
DataSplitTable34 <- DataSplitTable[!is.na(DataSplitTable$RUN34), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable34 [,myRespName])
myRespXY <- DataSplitTable34 [,c("x","y")]
DataSplitTable34 <- as.matrix(na.omit(DataSplitTable$RUN34))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable34,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal34_GLM", # name of this analysis
                            VarImport = 100) 
Bsal34_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal34_eval)
Bsal34_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal34_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal34_glm.img")

# Model 35;
DataSplitTable35 <- DataSplitTable[!is.na(DataSplitTable$RUN35), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable35 [,myRespName])
myRespXY <- DataSplitTable35 [,c("x","y")]
DataSplitTable35 <- as.matrix(na.omit(DataSplitTable$RUN35))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable35,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal35_GLM", # name of this analysis
                            VarImport = 100) 
Bsal35_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal35_eval)
Bsal35_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal35_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal35_glm.img")

# Model 36;
DataSplitTable36 <- DataSplitTable[!is.na(DataSplitTable$RUN36), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable36 [,myRespName])
myRespXY <- DataSplitTable36 [,c("x","y")]
DataSplitTable36 <- as.matrix(na.omit(DataSplitTable$RUN36))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable36,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal36_GLM", # name of this analysis
                            VarImport = 100) 
Bsal36_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal36_eval)
Bsal36_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal36_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal36_glm.img")

# Model 37;
DataSplitTable37 <- DataSplitTable[!is.na(DataSplitTable$RUN37), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable37 [,myRespName])
myRespXY <- DataSplitTable37 [,c("x","y")]
DataSplitTable37 <- as.matrix(na.omit(DataSplitTable$RUN37))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable37,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal37_GLM", # name of this analysis
                            VarImport = 100) 
Bsal37_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal37_eval)
Bsal37_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal37_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal37_glm.img")

# Model 38;
DataSplitTable38 <- DataSplitTable[!is.na(DataSplitTable$RUN38), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable38 [,myRespName])
myRespXY <- DataSplitTable38 [,c("x","y")]
DataSplitTable38 <- as.matrix(na.omit(DataSplitTable$RUN38))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable38,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal38_GLM", # name of this analysis
                            VarImport = 100) 
Bsal38_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal38_eval)
Bsal38_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal38_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal38_glm.img")

# Model 39
DataSplitTable39 <- DataSplitTable[!is.na(DataSplitTable$RUN39), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable39 [,myRespName])
myRespXY <- DataSplitTable39 [,c("x","y")]
DataSplitTable39 <- as.matrix(na.omit(DataSplitTable$RUN39))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable39,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal39_GLM", # name of this analysis
                            VarImport = 100) 
Bsal39_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal39_eval)
Bsal39_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal39_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal39_glm.img")

# Model 40;
DataSplitTable40 <- DataSplitTable[!is.na(DataSplitTable$RUN40), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable40 [,myRespName])
myRespXY <- DataSplitTable40 [,c("x","y")]
DataSplitTable40 <- as.matrix(na.omit(DataSplitTable$RUN40))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable40,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal40_GLM", # name of this analysis
                            VarImport = 100) 
Bsal40_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal40_eval)
Bsal40_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal40_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal40_glm.img")

# Model 41;
DataSplitTable41 <- DataSplitTable[!is.na(DataSplitTable$RUN41), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable41 [,myRespName])
myRespXY <- DataSplitTable41 [,c("x","y")]
DataSplitTable41 <- as.matrix(na.omit(DataSplitTable$RUN41))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable41,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal41_GLM", # name of this analysis
                            VarImport = 100) 
Bsal41_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal41_eval)
Bsal41_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal41_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal41_glm.img")

# Model 42;
DataSplitTable42 <- DataSplitTable[!is.na(DataSplitTable$RUN42), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable42 [,myRespName])
myRespXY <- DataSplitTable42 [,c("x","y")]
DataSplitTable42 <- as.matrix(na.omit(DataSplitTable$RUN42))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable42,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal42_GLM", # name of this analysis
                            VarImport = 100) 
Bsal42_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal42_eval)
Bsal42_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal42_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal42_glm.img")

# Model 43;
DataSplitTable43 <- DataSplitTable[!is.na(DataSplitTable$RUN43), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable43 [,myRespName])
myRespXY <- DataSplitTable43 [,c("x","y")]
DataSplitTable43 <- as.matrix(na.omit(DataSplitTable$RUN43))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable43,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal43_GLM", # name of this analysis
                            VarImport = 100) 
Bsal43_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal43_eval)
Bsal43_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal43_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal43_glm.img")

# Model 44;
DataSplitTable44 <- DataSplitTable[!is.na(DataSplitTable$RUN44), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable44 [,myRespName])
myRespXY <- DataSplitTable44 [,c("x","y")]
DataSplitTable44 <- as.matrix(na.omit(DataSplitTable$RUN44))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable44,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal44_GLM", # name of this analysis
                            VarImport = 100) 
Bsal44_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal44_eval)
Bsal44_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal44_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal44_glm.img")

# Model 45;
DataSplitTable45 <- DataSplitTable[!is.na(DataSplitTable$RUN45), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable45 [,myRespName])
myRespXY <- DataSplitTable45 [,c("x","y")]
DataSplitTable45 <- as.matrix(na.omit(DataSplitTable$RUN45))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable45,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal45_GLM", # name of this analysis
                            VarImport = 100) 
Bsal45_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal45_eval)
Bsal45_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal45_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal45_glm.img")

# Model 46;
DataSplitTable46 <- DataSplitTable[!is.na(DataSplitTable$RUN46), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable46 [,myRespName])
myRespXY <- DataSplitTable46 [,c("x","y")]
DataSplitTable46 <- as.matrix(na.omit(DataSplitTable$RUN46))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable46,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal46_GLM", # name of this analysis
                            VarImport = 100) 
Bsal46_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal46_eval)
Bsal46_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal46_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal46_glm.img")

# Model 47;
DataSplitTable47 <- DataSplitTable[!is.na(DataSplitTable$RUN47), ]
myRespName <- 'Bsal'
myResp <- as.numeric(DataSplitTable47 [,myRespName])
myRespXY <- DataSplitTable47 [,c("x","y")]
DataSplitTable47 <- as.matrix(na.omit(DataSplitTable$RUN47))
myBiomodData_glm <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = env_glm, # predictors
                                         resp.xy = myRespXY, # xy data
                                         resp.name = myRespName, #species name
                                         PA.nb.rep = 0, # no pseudoabsences
                                         na.rm = TRUE) # remove points with NA
myBiomodOption <- BIOMOD_ModelingOptions()
Bsal_glm <- BIOMOD_Modeling(myBiomodData_glm, # input
                            models = c('GLM'), # algorithm
                            models.options = myBiomodOption, # options
                            DataSplitTable = DataSplitTable47,
                            models.eval.meth = c('TSS'), # eval metric
                            modeling.id = "Bsal47_GLM", # name of this analysis
                            VarImport = 100) 
Bsal47_eval <- get_evaluations(Bsal_glm)
Bsal_eval_glm <- rbind(Bsal_eval_glm, Bsal47_eval)
Bsal47_var <- get_variables_importance(Bsal_glm)
Bsal_var_glm <- cbind(Bsal_var_glm, Bsal47_var)
Bsal_glm_ens <- BIOMOD_EnsembleModeling(modeling.output = Bsal_glm, 
                                        chosen.models = 'all', 
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = c('TSS'),
                                        eval.metric.quality.threshold = c(0.7),
                                        models.eval.meth = c('TSS'),
                                        prob.mean.weight = TRUE,
                                        prob.mean.weight.decay = 'proportional')
Bsal_glm_proj <- BIOMOD_Projection(modeling.output = Bsal_glm, # input
                                   new.env = env_glm, # projection background
                                   proj.name = 'glm_proj', # output folder
                                   selected.models = 'all', # all models
                                   binary.meth = 'TSS',
                                   compress = FALSE,
                                   build.clamping.mask = FALSE)
Bsal_glm_projensemble <- BIOMOD_EnsembleForecasting(projection.output = Bsal_glm_proj, EM.output = Bsal_glm_ens)
plot(Bsal_glm_projensemble)
glm_predictions <- get_predictions(Bsal_glm_projensemble)
writeRaster(subset(glm_predictions, 2), filename ="Bsal47_glm.img")

glm_predictions <- list.files(path = ".R_biomod2/GLM output", pattern ='\\.img$', full.names = TRUE) 
glm.stack <- stack(glm_predictions)
weighted.mean(glm.stack, 
              Bsal_eval_glm$Testing.data.GLM.RUN1.AllData, 
              na.rm = FALSE, 
              filename = "glm_weighted.img")
write.csv(Bsal_eval_glm, file = "Bsal_eval_glm.csv")
write.csv(Bsal_var_glm, file = "Bsal_var_glm.csv")
