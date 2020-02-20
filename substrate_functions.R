#----------------------------------------------------------------------------
# Script:  substrate_functions.R
# Purpose: Support building and evaluation of Random forst models for substrate paper.
# Created: January 2020. EJG
# Goal: Adapt the necessary functionality from Cole's Fit_Random_Forest.R script 
#   to use ranger() library and expand to do independent data evaluation.
#
# Notes:
#  - 
#----------------------------------------------------------------------------

#-- Load packages using neat loading code from Cole ... 
#       Currently not using: ggplot2, spatialEco, xlsx, robustbase, biomod2

# check for any required packages that aren't installed and install them
required.packages <- c("diffeR", "vegan", "randomForest", "ranger", "raster", "rgdal", "stringr", "measures", "e1071",   
                       "caret", "tidyr","dplyr", "PresenceAbsence")
uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)


#----------------------------------------------------------------------------
#-- Data sources and other constants ...  

#-- As per Cole's reporting code, set an output directory and initiate the log file.
#       Directory will be created if doesn't exist; file will be overwritten if it does.
output.dir    <- 'C:/Users/Edward/Dropbox/DFO Job/Substrate2019/Results'
model.dir     <- 'C:/Users/Edward/Dropbox/DFO Job/Substrate2019/Models'
predictor.dir <- 'C:/Data/SpaceData/Substrate2019/Predictors/'

# filepath to predictor directory (rasters must be tif format)
predictors.coastwide <- 'C:/Data/SpaceData/Substrate2019/Predictors/Coastwide'
bioregions <- c('HG','NCC','WCVI','QCS','SOG')

#-- Regression formulas
#   NOTE: Fetch not used coastwide because of deeper shoaling depth. Fetch not in the Obs database.
#   Cole has nice code to make formulas, but I am a hack here ... 
shore.formula <- "BType4 ~ bathy + broad_BPI + circulation + curvature + fine_BPI + med_BPI + rugosity + sd_slope + slope + tidal + fetch"
coast.formula <- "BType4 ~ bathy + broad_BPI + circulation + curvature + fine_BPI + med_BPI + rugosity + sd_slope + slope + tidal"


# Source of obs data. Obtained Dec 2020. 
# Working layers modifed in GIS to ensure BType4 and Region fields added. 
sourceGDB <- 'C:/Data/SpaceData/Substrate2019/IDE work.gdb'
ptSources <- list( 'obs', 'BHM_habitat', 'DropCam_reg', 'ROV_reg')
names(ptSources) <- c('train', 'dive', 'cam', 'ROV')


# proj4 string for albers projection with NAD83 datum
spat.ref <- '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

#-- Log file to be placed in output.directory 
log.file <- 'IDE_test_log.txt'

#----------------------------------------------------------------------------
#-- Set relevant RF constants. 

ntree = 1000
repl  = TRUE
imp.1 = TRUE
imp.2   = 'impurity'
test.frac <- 0.6

#----------------------------------------------------------------------------
#-- Functions.   
#   Results.Row()
#   Load.Test.Data(), Load.Predictors(), Load.Obs.Data() 
#   Test.Sample.Size()
#----------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Calculate the reporting statistics for an RF model (testModel), 
# based on a set of observations (testData). 
# Requires: 
#   testModel = a weighted RF model built with ranger()
#   testDat = dataframe with correct predictors and observed BType
#   testData must have the same predictors attached as those used to build testModel.
# Returns: a data frame with a single row. suitable for rbind().
Results.Row <- function( testModel, testData ){
  
  mant <- 3
  y <- predict( testModel, testData )
  z <- caret::confusionMatrix( y$predictions, testData$BType4 )
  ber <- BER(testData$BType4, y$predictions)
  prev <- Prev.Balance( summary(testData$BType4) / sum( summary(testData$BType4)) )
  
  out <- data.frame( 
    "N"         = nrow(testData),
    "Imbalance" = round( prev, mant ), 
    "Accuracy"  = round( as.numeric( z$overall[ 'Accuracy' ]), mant ),
    "Kappa"     = round( as.numeric( z$overall[ 'Kappa' ]), mant ), 
    "BERneg"    = round( 1-ber, mant ), 
    "OOB"       = round( testModel$prediction.error, mant ),
    round( Wtd.Stats( z ), mant ),
    round( diffr.Stats( z$table ), mant )
  )
  return(out)
}


#------------------------------------------------------------------
#-- Calculate the quantity and allocation components of the matrix.
#   Uses diffeR package and follows Pontius and Santacruz (2014).
# Takes: a caret contingency table - 4x4 hard-coded.
# Returns: A one-row data.frame of named statistics
# Notes: This will need to be updated based on detailed review of the 
#   available statistics and of Stehman and Foody (2019).
diffr.Stats <- function( cTable ){
  
  # First need to convert to something called a ctmatrix
  # Probably necessary because source is as factors.
  x <- matrix( cTable, 4,4 )
  rownames( x ) <- c(1,2,3,4); colnames( x ) <- c(1,2,3,4)
  
  y <- diffTablej( x, digits = 0, analysis = "error" )
  
#  Other stats available. Especially by class. 
#  pop <- matrix( c( 1,2,3,4, colSums( x)), 4, 2)
#  y2  <- categoryComponentsPlot( ctmatrix = x, population = pop )
  
  out <- y[ y$Category == 'Overall', -1 ]
  
  return( out )
}


#----------------------------------
#-- Weighted Sensitivity, Specificity, and TSS calculations
#   TSS calculation taken from Fit_Random_Forest.R.
#   https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2664.2006.01214.x
#   TSSwtd defined as tss/prevalence summed across classes. 
# Takes: A caret contingency table
# Returns: A one-row data.frame of named statistics
Wtd.Stats <- function( cTable ){
  
  x <- as.data.frame( cTable$byClass )
  TSS <- x$Sensitivity + x$Specificity - 1
  TSS <- sum( TSS * x$Prevalence )
  
  Sens <- sum( x$Sensitivity * x$Prevalence )
  Spec <- sum( x$Specificity * x$Prevalence )
  
  return( data.frame( 'TSSWtd' = TSS, 'SensWtd' = Sens, 'SpecWtd' = Spec ))
}


#----------------------------------
#-- Calculate the imbalance in prevalence.
# Takes: list of prevalence values
# Returns: Generalized statistic to describe deviance from prevalence.
#   0 is fully balanced. Max = N/number of categories
# Notes: Current statistic is mean global difference among classes.
#        Can diversity index serve as proxy for  prevalence?
#         diversity(x, index = "shannon", MARGIN = 1, base = exp(1))
Prev.Balance <- function( plist ){
  tot <- 0
  denom <- 0
  N <- length(plist)
  for (a in 1:N ){
    for (b in 2:N){
      if (b > a ){
      tot <- tot + abs( plist[a] - plist[b])
      denom <- denom + 1
      }
    }}
  return( as.numeric( tot/denom ))
}


#--------------------------------------------------------
# Pulls all the ID set data together for coastal tests ... 
# USES train.data as a global variable. 
Assemble.Coast.Test <- function( indat ){
  x <- rbind( indat[['NCC']], indat[['HG']], indat[['WCVI']], indat[['QCS']], indat[['SOG']] )
  return( x )
}


#------------------------------------------
# Load all the observational data, obs plus the ID sets. 
# Just a wrapper function to make a bunch of calls.
Load.Observations <- function(){
  
  #-- Load Coastwide observations (created by Cole's code)
  #   Provided here as a shape file, which I copied to the working gdb file.
  obs.in <- Load.Obs.Data( sourceGDB, ptSources$train )

    #-- Load independent data sets. 
  dive.in  <- Load.Test.Data( sourceGDB, ptSources$dive, 
                              c("SourceDB", "Key", "DepthCat", "SubCat", "SubSubCat", "BType4") )
  dive.in$BType4 <- factor(dive.in$BType4)
  
  cam.in  <- Load.Test.Data( sourceGDB, ptSources$cam,
                             c("Survey", "Transect", "RMSM_cat", "Substrate1", "Substrate2", "BType4") )
  cam.in$BType4 <- factor(cam.in$BType4)
  
  rov.in  <- Load.Test.Data( sourceGDB, ptSources$ROV,
                             c("Survey", "Transect", "RMSM_cat", "Substrate1", "Substrate2", "BType4") )
  rov.in$BType4 <- factor(rov.in$BType4)
  
  foo <- list( obs.in, dive.in, cam.in, rov.in )
  names(foo) <- names( ptSources )
  
  return(foo)
}
  
  
#-----------------------------
# Load COASTWIDE predictors. 
Add.100m.Data <- function( obsList ){
  
  rasters <- Load.Predictors( predictors.coastwide )
  preds  <- raster::extract(rasters, obsList, df = TRUE)
  merged <- cbind( as.data.frame(obsList), preds)
  return(merged)
}


#----------------------------
# Load REGIONAL predictors.
# Grabs the predictor values from the 20 m rasters. By region.
# Requires: list of points, fewer attributes is better! 
# Global var: bioregions 
Add.20m.Data <- function( obsList ){

  outList <- vector("list", length = length(bioregions) )
  names( outList ) <- bioregions
  
  for (i in bioregions ) {
    rasters <- Load.Predictors( paste0( predictor.dir, "/", i ))
    preds  <- raster::extract(rasters, obsList, df = TRUE)
    # remove awkward fields ... 
    preds <- preds[, !names(preds) %in% c('ID') ]
    merged <- cbind( as.data.frame(obsList), preds)
    outList[[ i ]] <- drop_na( merged )
    
    cat( i, ' loaded ...\n')
  }
  return( outList )
}


#----------------------------------------------------------------------------
# Loads predictors from specified subdirectory 
Load.Predictors <- function( pred.dir ) {
  
  # make a list of predictor rasters
  raster.list <- list.files(path = pred.dir, pattern = '\\.tif$', full.names = TRUE)
  
  # make list of raster names (without extension)
  raster.names <- lapply(raster.list, FUN = function(raster.layer){
    substr(basename(raster.layer), 1, nchar(basename(raster.layer)) - 4)
  } )
  
  # create a raster stack from raster_list
  raster.stack <- raster::stack(x = raster.list)
  
  return(raster.stack)
}


#-------------------------------------------------------------------
# A simple function to report the differences in sample sizes before 
# and after adding predictor data. 
Diff.Sets <- function( x, y){
  a <- dim( x )[[1]]
  b <- sum( unlist( lapply( y, function(z) dim(z)[[1]] ) ))
  cat(  a, '-', b, "=", a-b, '(', round((a-b)/a, 4)*100, '% loss)' )
}


#----------------------------------------------------------------------------
# Load test data from specified geodatabase and feature class. Preserves the columns. 
# BType4 is required field, containing the classifed idendependent test data. 
Load.Test.Data <- function( gdb, fc, columns ) {
  
  # for reading into ESRI File Geodatabase
  subset(ogrDrivers(), grepl('GDB', name))
  
  # read the feature class specified by user
  GIS.data <- readOGR(dsn = gdb, layer = fc)
  
  # Keep only COLUMNS specified from the GIS POINTS DATA FRAME ... 
  observations <- GIS.data[, (names(GIS.data) %in% columns)] 
  
  # make X, Y columns (even if already exist to recreate spatial points dataframe later)
  observations$X <- observations@coords[, 1]
  observations$Y <- observations@coords[, 2]
  
  # number of observations before dropping invalid BType4 values
  observations.before <- nrow(observations)
 
  # make sure that BType4 column has only values 1:4
  observations <- observations[observations$BType4 %in% c(1, 2, 3, 4), ]
  
  # number of observations after dropping any records with invalid BType4 values
  observations.after <- nrow(observations)
  
    # write a summary to the log file ... a BIT UNORTHODOX for a function perhaps. 
  sink( file = file.path( output.dir, log.file ), append=TRUE)
  
  cat('\n\n\n***********************************************\n\n\n')
  cat(paste0('\nNumber of observations before omitting records with invalid BType4 values:\n', observations.before))
  cat(paste0('\n\nNumber of observations omitted due to records with invalid BType4 values:\n', observations.before - observations.after))
  # percent of observations dropped
  cat(paste0('\n\n', round(observations.after / observations.before * 100, 2), '% of observations have been kept after omitting records with invalid BType4 values'))
  
  # stop sink to log file
  sink()
  
  return( observations )
}


#----------------------------------------------------------------------------
# Loads observations from specified geodatabase and feature class
Load.Obs.Data <- function( gdb, fc ) {

  # for reading into ESRI File Geodatabase
  subset(ogrDrivers(), grepl('GDB', name))
  
  # read the feature class specified by user
  obs <- readOGR(dsn = gdb, layer = fc)
  obs$BType4 <- factor( obs$BType4 )
  
  # Columns to keep ... 
  # columns <- c("BType4", "bathy", "brd_BPI", "circltn", "curvatr", "fetch", "fin_BPI", "med_BPI", 
  #              "rugosty", "sd_slop", "slope", "tidal") 
  # 
  # # convert to data frame with no spatial attributes ... 
  # obs <- obs[, (names(obs) %in% columns)] 
  
  return(obs)
}


#---------------------------------------------------------------------------------------------
# For the provided  sample: 1) partition the data; 2) build train model, 3) build a full model
# build a RF model with partitioned subset data
# evaluate model with withheld data (should be same sample size ans i-1 model )
Test.Sample.Size <- function( obs, howMany, part, theForm ){
  
  N <- nrow( obs )
  results <- data.frame()
  
  for (i in 1:howMany) {
    
    # pull the sample ... 
    RF.sample <- obs[ sample( 1:N, (i/howMany * N) ), ]
    
    # calculate the wts ...  
    props <- RF.sample %>% group_by(BType4) %>% count(BType4) 
    #prop.Btype4 <- obs@data %>% group_by(BType4) %>% count(BType4) 
    pWts <- 1 - ( props$n / sum( props$n ))
    
    # Build the  model ... 
    anRF <- ranger( theForm,
                    data = RF.sample,
                    num.trees = ntree, replace = repl, importance = imp.2,
                    case.weights = pWts[ RF.sample$BType4 ])    #Define case.weights
    
    #-- Report on RF.Train vs. itself.
    x <- kappa(anRF$confusion.matrix)$coef
    y <- BER(RF.sample$BType4, anRF$predictions)
    
    r1 <- data.frame("test" = i, "sample" = nrow(RF.sample), "evaluation" = "full", 
                     "kappa" = round(x, 4), "BER" = round(y, 4)  )
    
    results <- rbind( results, r1) 
    
    #--------------------------------------
    # Now partition the data, and report statistics for testing partition.
    x <- sample( 1:nrow(RF.sample), part * nrow(RF.sample) )
    train <- RF.sample[ x, ]
    test  <- RF.sample[ -x, ]
    
    cat("train N", nrow(train), "\n") 
    # Build the weighting vector ...
    prop.Btype4 <- train %>% group_by(BType4) %>% count(BType4) 
    wts <- 1 - ( prop.Btype4$n / sum( prop.Btype4$n ))
    
    # Build the coastwide model ... 
    anRF <- ranger( theForm,
                    data = train,
                    num.trees = ntree, replace = repl, importance = imp.2,
                    case.weights = wts[ train$BType4 ])    #Define case.weights
    
    #-- Report on Train vs. Testing.
    test.preds <- predict( anRF, test )
    
    x <- kappa( confusionMatrix(test.preds$predictions, test$BType4)$table )
    y <- BER(test$BType4, test.preds$predictions)
    
    cat( x$coef, '\n')
    cat( y, '\n')
    
    r1 <- data.frame("test" = i, "sample" = nrow(RF.sample), "evaluation" = "partition", 
                     "kappa" = round(x$coef, 4), "BER" = round(y, 4)  )
    
    results <- rbind( results, r1) 
    
    #save( anRF, file = file.path(output.directory, paste0( i, "_", 'msm_model.RData')) )  
  }
  return( results ) 
}


### fin.
