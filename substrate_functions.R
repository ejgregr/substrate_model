#----------------------------------------------------------------------------
# Script:  substrate_functions.R
# Created: January 2020. EJG
#
# Purpose: Support building and evaluation of Random forst models for substrate paper.
# Goal: Adapt the necessary functionality from Cole's Fit_Random_Forest.R script 
#   to use ranger() library and expand to do independent data evaluation.
#
# Notes:
#  - 
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
#-- Load packages using neat loading code from Cole ... 

# check for any required packages that aren't installed and install them
required.packages <- c("ggplot2", "reshape2", "tidyr","dplyr", 
                       "diffeR", "vegan", "randomForest", "ranger", "raster", "rgdal", "stringr",
                       "measures", "e1071", "caret", "PresenceAbsence", 
                       "superheat", "PNWColors")
# Currently not using:  "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr",
#                       "tinytex", "rmarkdown", "binr"

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

#============================================================================
#-- Functions.   


#--------------------------------------------------------------
# Builds a random forest model for the specified region
# Returns a list with the final model fitted to all the data, 
# and a summary table of the re-sampled performance statistics.
# Requires - list of model names
# Calls    - Results.Row()
Make.Ranger.Model <- function( x.data, x.formula, partition = 0.7, iterations = 1  ){
  
  y <- nrow( x.data )
  out.table <- NULL
  
  for (i in 1:iterations ){
    
    #-- build the train/test partition ... 
    idx <- sample( 1:y, round(partition * y), replace = FALSE )
    z <- replace( rep(0, y), idx, 1)
    xx <- cbind(x, 'partition' = z )
    
    x.train <- xx[ xx$partition == 1, ]
    x.test  <- xx[ xx$partition == 0, ]
    
    #-- build a weighted model ... 
    props <- x.train %>% group_by(BType4) %>% count(BType4) 
    wts <- 1 - ( props$n / sum( props$n ))
    x.model <- ranger( x.formula,
                       data = x.train,
                       num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                       case.weights = wts[ x.train$BType4 ])    #Define case.weights
    
    #-- Evaluate with testing partition  ...
    out.table <- rbind( out.table, Results.Row( x.model, x.test )) 
    print(i)
  }
  
  #-- Build the model result entry.
  if (iterations > 1)
    z <- Summary.Row( out.table[, -1] )
  else {
    z <- out.table
  }
  #-- Build the model with all the training data.
  props <- x.data %>% group_by(BType4) %>% count(BType4) 
  wts <- 1 - ( props$n / sum( props$n ))
  x.model <- ranger( x.formula,
                     data = x.data,
                     num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                     case.weights = wts[ x.data$BType4 ])    #Define case.weights
  
    return( list( 'Stats' = z, 'Model' = x.model ))
}

# Function to predict a raster surface, write to disk, and export a png of the surface
# env.predictors: raster stack of environmental predictors, ranger.model: model object, output.directory: path to output
Predict.Surface <- function(env.predictors, ranger.model, output.directory){
  
  # create output directory - this creates it in the working directory. If it exists, do not show warning
  dir.create(output.directory, showWarnings = FALSE)
  
  # save ranger model to disk as RData file
  save(ranger.model, file = file.path(output.directory, 'ranger_model.RData'))
  
  # predict substrate using raster stack and model file
  raster.obj <- raster::predict(object   = env.predictors,  
                                model    = ranger.model,               
                                progress = 'text',
                                fun = function(model, ...) predict(model, ...)$predictions)
  
  # write raster file to disk
  writeRaster(raster.obj, file.path(output.directory, 'classified_substrate.tif'), format = 'GTiff', datatype = 'INT2S')
  
  # generate table of proportions for each class in predicted raster
  raster.prop <- round(100 * prop.table(table(na.omit(as.data.frame(getValues(raster.obj))))), 2)
  
  # colour palette for map
  pal <- c("#999999", "#33bbff", "#ffff99", "#ffb84d")
  
  # labels for legend
  labels <- c("Rock", "Mixed", "Sand", "Mud")
  
  # Map (up to 5,000,000 pixels)
  png(file=file.path(output.directory, "substrate_raster.png"),
       height = 7, width = 6, units = "in", res = 400)
  plot(raster.obj, maxpixels=5000000, col=pal, legend=FALSE,
       xlab = "Easting", ylab = "Northing", cex.axis = .5, cex.lab = .75)
  legend(x = "bottomleft",
         legend = labels, fill = pal, title=NA, bg = NA, box.col = NA)
  
  dev.off()
  
  return(list(raster.obj, raster.prop))
  
}

       
#------------------------------------------------------------------------------------
# Calculate the reporting statistics for an RF model (testModel), 
# based on a set of observations (testData). 
# Requires: 
#   testModel: a weighted RF model built with ranger()
#   testData : dataframe with correct predictors and observed BType
#    testData must have the same predictors attached as those used to build testModel.
# Calls:  diffr.Stats()
#         Prev.Balance()
#         Wtd.Stats()
# Returns: a data frame with a single row. suitable for rbind().
Results.Row <- function( testModel, testData ){
  
  mant <- 3
  y <- predict( testModel, testData )
  z <- caret::confusionMatrix( y$predictions, testData$BType4 )
  ber <- BER(testData$BType4, y$predictions)
  prev <- Prev.Balance( summary(testData$BType4) / sum( summary(testData$BType4)) )
  
  out1 <- data.frame( 
    "N"         = nrow( testData),
    "Imbalance" = round( prev, mant ), 
    "OOB"       = round( testModel$prediction.error, mant ),
    "Accuracy"  = round( as.numeric( z$overall[ 'Accuracy' ]), mant ),
    "BERneg"    = round( 1-ber, mant ), 
    "TSS"       = TSS.Calc( z$table ),
    round( Wtd.Stats( z ), mant ),
    round( diffr.Stats( z$table )/nrow(testData), mant )
  )
  out2 <- data.frame( 
    "User"  = round( diag(z$table) / rowSums(z$table), mant ), 
    "Prod"  = round( diag(z$table) / colSums(z$table), mant ),
    "PrevObs"  = round( as.vector( table( testData$BType4 )) ),
    "PrevPred" = round( as.vector( table( y$predictions )) ))
  
  return( list( 'Integrated' = out1, 'PerClass' = t(out2) ))
}


# #-- Testing Results.Row ...
# Results.Row( foo$Model, x )
# y <- predict( foo$Model, x )
# rows(x)
# y$predictions
# z <- caret::confusionMatrix( y$predictions, x$BType4 )
# z$table
# dim(x)[[1]]
# round( diag(z$table/dim(x)[[1]]) / rowSums(z$table/dim(x)[[1]]), 3)


#---------------------------------------
# Summarizes a data.frame of statistics
# RETURNS mean and standard deviation. 
Summary.Row <- function( inRes ){
  
  x <- round( apply( inRes, 2, mean ), 3)
  y <- round( apply( inRes, 2, sd ), 4)
  
  z <- rbind(x,y)
  row.names( z ) <- c('Mean', 'StdDev')
  return( z )
  
  #Nice but returns a string
  #formatC(y, format = "e", digits = 2)
}


#----------------------------------------------------------------
#-- Calculate the True Skill Statistic for a 4x4 matrix following
#   Allouche et al. 2006. Applies 2 steps - first calculates the 
#   Expected correct classifications expected by chance; 
#   then substracts this from the diagonal values, and reports the 
#   sum of the diagonal values. 
# Assumes: 4x4 matrix
# Takes: a caret contingency table
# Returns: A single TSS statistic for the entire table.
# NOTES: Correcting for chance success reduces the sample size on 
#   the diagonal ... 
TSS.Calc <- function( cTable ){
  x <- matrix( cTable, 4,4 )

  # Expected correct predictions by chance (diagonal)
  E <- rowSums( x ) * colSums( x ) / sum( x )
  # Correct predictions corrected for success due to chance 
  R <- diag(x) - E
  
  # The diagonal values from a perfect forecast 
  n.star <- rowSums( x )
  
  # Perfect forecast corrected for success due to chance 
  R.star <- n.star - ( n.star^2/sum(x) )
    
  return( round( sum(R) / sum(R.star), 3) )
    
}


# #-- Testing TSS Calculation ...
# test <- matrix( c(10,0,0,0, 0,10,0,0, 0,0,10, 0, 0, 0, 0, 10), 4, 4 )
# test <- matrix( c(.25,0,0,0, 0,.25,0,0, 0,0,.25, 0, 0, 0, 0, .25), 4, 4 )
# 
# #-- A simplee bad model ...
# test <- matrix( c(0,10,0,0, 0,0,0,10, 10,0,0,0, 0,0,10,0), 4, 4 )
# 
# #-- Shouldn't these have greater badness?
# test <- matrix( c(0,100,0,0, 0,0,0,10, 10,0,0,0, 0,0,10,0), 4, 4 )
# test <- matrix( c(0,33,33,34, 10,0,10,10, 10,10,0,10, 10,10,10,0), 4, 4 )
# 
# test
# TSS.Calc(test)
# 
# #-- Random matrix has values around 0 which is good ...
# test <- matrix( sample(1:100, 16), 4, 4)
# TSS.Calc(test)
# 


#------------------------------------------------------------------
#-- Calculate the quantity and allocation components of the matrix.
#   Uses diffeR package and follows Pontius and Santacruz (2014).
# Assumes: 4x4 matrix
# Takes: a caret contingency table
# Returns: A one-row data.frame of named statistics
# Notes: This will need to be updated based on detailed review of the 
#   available statistics and Stehman and Foody (2019).
#   diffTablej() returns Omission, Agreement, Comission, Quantity, Exchange, and Shift metrics
#     by classes and overall.
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
  return( out[ c('Quantity', 'Exchange', 'Shift') ] ) #return desired stats and order
}

# #-- Testing diffeR stats ... 
# z <- caret::confusionMatrix( rf.region.HG$predictions, train.data.20m$HG$BType4)
# diffr.Stats( z$table )
# 
# x <- matrix( z$table, 4,4 )
# rownames( x ) <- c(1,2,3,4); colnames( x ) <- c(1,2,3,4)
# y <- diffTablej( x, digits = 0, analysis = "error" )



#----------------------------------
#-- Weighted Specificity (true negative rate), and TSS calculations
#   TSS calculation taken from Fit_Random_Forest.R and based on Allouche et al. 2006.
#   TSSwtd defined as tss/prevalence summed across classes. 
# Takes:   A caret contingency table
# Returns: A one-row data.frame of named statistics
Wtd.Stats <- function( cTable ){
  
  x <- as.data.frame( cTable$byClass )
  TSS <- x$Sensitivity + x$Specificity - 1
  TSS <- sum( TSS * x$Prevalence )
  
  #-- Dropped cuz this looks like the accuracy calculation - 2020/02/27.
  #Sens <- sum( x$Sensitivity * x$Prevalence )
  Spec <- sum( x$Specificity * x$Prevalence )
  
  return( data.frame( 'TSSWtd' = TSS, 'TNRWtd' = Spec ))
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
  
  ptSources <- list( 'obs', 'BHM_habitat', 'DropCam_reg', 'ROV_reg')
  names(ptSources) <- c('train', 'dive', 'cam', 'ROV')
  
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
# For the provided sample: 1) partition the data; 2) build train model, 3) build a full model
# build a RF model with partitioned subset data
# evaluate model with withheld data (should be same sample size ans i-1 model )
Test.Sample.Size <- function( obs, howMany, part, theForm ){
  
  N <- nrow( obs )
  results <- NULL
  
  for (i in 1:howMany) {
    
    # pull the sample ... 
    RF.sample <- obs[ sample( 1:N, (i/howMany * N) ), ]
    
    # partition the sample ... 
    x <- sample( 1:nrow(RF.sample), part * nrow(RF.sample) )
    x.train <- RF.sample[ x, ]
    x.test  <- RF.sample[ -x, ]
    
    cat("train N", nrow(x.train), "\n") 
    # Build the weighting vector ...
    prop.Btype4 <- x.train %>% group_by(BType4) %>% count(BType4) 
    wts <- 1 - ( prop.Btype4$n / sum( prop.Btype4$n ))
    
    # Build the coastwide model ... 
    anRF <- ranger( theForm,
                    data = x.train,
                    num.trees = ntree, replace = repl, importance = imp.2,
                    case.weights = wts[ x.train$BType4 ])    #Define case.weights
    
    #-- Evaluate the model using the test partition.
    results <- rbind( results, Results.Row( anRF, x.test )) 
    print(i)
    
    }
  return( results ) 
}


#---------------------------------------------------------------------------------------------
# For the provided sample: 1) partition the data ranging from even prevalence to imbalance = NN;
# 2) build an model with training data; 3) evaluate the model with testing partition.
# REQUIRES: obs class attribute = BType; Needs to have 4 classes.
Test.Prevalence <- function( obs, N, part, theForm ){
  
  results <- NULL

  #-- Build the partion sizes first  ...
  #   NEED to 1) fix sample size and 
  #           2) change prevlance of one category only to avoid conflating things ... 

  #-- increase prevalence of rock from 0.25 to 0.85 ... 
  w <- list()
#  i = 0
  for (i in seq(0, 0.6, by=0.05) ) {

    w[1] <- 0.25 - i/3
    w[2] <- 0.25 + i
    w[3] <- 0.25 - i/3
    w[4] <- 0.25 - i/3
    
    cat( paste( w[1], w[2], w[3], w[4] ), '\n')

    RF.sample <- rbind( 
      b1 <- sample_n( obs[ obs$BType4 == 1, ], round( w[[1]] * N )),
      b2 <- sample_n( obs[ obs$BType4 == 2, ], round( w[[2]] * N )),
      b3 <- sample_n( obs[ obs$BType4 == 3, ], round( w[[3]] * N )),
      b4 <- sample_n( obs[ obs$BType4 == 4, ], round( w[[4]] * N ))
    )  

    # partition the sample ... 
    x <- sample( 1:nrow(RF.sample), part * nrow(RF.sample) )
    x.train <- RF.sample[ x, ]
    x.test  <- RF.sample[ -x, ]
    
    cat("train N", nrow(x.train), "\n") 
    # Build the weighting vector ...
    prop.Btype4 <- x.train %>% group_by(BType4) %>% count(BType4) 
    wts <- 1 - ( prop.Btype4$n / sum( prop.Btype4$n ))
    
        # Build the coastwide model ... 
    anRF <- ranger( theForm,
                    data = x.train,
                    num.trees = ntree, replace = repl, importance = imp.2,
                    case.weights = wts[ x.train$BType4 ])    #Define case.weights
    
    #-- Evaluate the model using the test partition.
    results <- rbind( results, Results.Row( anRF, x.test )) 
    print(i)
    
  }
  return( results ) 
}


### fin.
