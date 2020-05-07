#----------------------------------------------------------------------------
# Script:  substrate_functions.R
# Created: January 2020. EJG
#
# Purpose: Support building and evaluation of Random forst models for substrate paper.
# Goal: Adapt the necessary functionality from Cole's Fit_Random_Forest.R script 
#   to use ranger() library and expand to do independent data evaluation.
#
# Notes:
#  - 2020/05/01: Some serious pruning of the data load functions. Combined with 
# consolidating of data cleaning. 
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
results.dir   <- 'C:/Users/Edward/Dropbox/DFO Job/Substrate2019/Results'
model.dir     <- 'C:/Users/Edward/Dropbox/DFO Job/Substrate2019/Models'
predictor.dir <- 'C:/Data/SpaceData/Substrate2019/Predictors'
coastPred.dir <- 'C:/Data/SpaceData/Substrate2019/Predictors/Coastwide'
source.dir    <- 'C:/Data/SpaceData/Substrate2019'

# Sources of all the necessary spatial data including regions and observations.
# Data last updated May 30 2020. 
source.files  <- list( 'Obs'     = 'obs.shp',
                       'Dive'    = 'IndependentData/All__BHM_SpatPts.shp',
                       'Cam'     = 'IndependentData/All_DropCam_SpatPts.shp',
                       'ROV'     = 'IndependentData/All_ROV_SpatPts.shp',
                       'Regions' = 'regions/BC_coast_regions.shp' )

# filepath to predictor directory (rasters must be tif format)

bioregions <- c('HG','NCC','WCVI','QCS','SOG')

#-- Regression formulas
#   NOTE: Fetch not used coastwide because of deeper shoaling depth. Fetch not in the Obs database.
#   Cole has nice code to make formulas, but I am a hack here ... 
shore.formula <- "BType4 ~ bathy + broad_BPI + circulation + curvature + fine_BPI + med_BPI + rugosity + sd_slope + slope + tidal + fetch"
coast.formula <- "BType4 ~ bathy + broad_BPI + circulation + curvature + fine_BPI + med_BPI + rugosity + sd_slope + slope + tidal"

# proj4 string for albers projection with NAD83 datum
spat.ref <- '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

#-- Log file to be placed in output.directory 
log.file <- 'IDE_test_log.txt'

#-- Relevant RF constants. (Had 2 imps but one was for Random Forest)
RFC <- list( 'ntree'=1000, 'repl'=TRUE, 'imp'='impurity', 'test.frac'=0.6 )


#============================================================================
#-- Functions.   

#=============== ANALYSIS ===================


#---------------------------------------------------------------------------------------
# Builds a random forest model using ranger() using the provided train and test samples.
# 2020/04/13: Created to surface train and testing subsaamples. 
# Returns a list with the final fitted model
# and a summary table of performance statistics.
# Requires - list of model names
# Calls    - Results.Row()
Fixed.Ranger.Model <- function( x.train, x.test, x.formula ){
  
  out.table <- NULL
  
  #-- build a weighted model ... 
  props <- x.train %>% group_by(BType4) %>% count(BType4) 
  wts <- 1 - ( props$n / sum( props$n ))
  x.model <- ranger( x.formula,
                     data = x.train,
                     num.trees = RFC$ntree, replace = RFC$repl, importance = RFC$imp, oob.error = T,
                     case.weights = wts[ x.train$BType4 ])    #Define case.weights
  
  #-- Evaluate with testing partition  ...
  out.table <- rbind( out.table, Results.Row( x.model, x.test )) 

  return( list( 'Stats' = out.table, 'Model' = x.model ))
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
    "TSS"       = TSS.Calc( z$table ),
    "Accuracy"  = round( as.numeric( z$overall[ 'Accuracy' ]), mant ),
#    "BERneg"    = round( 1-ber, mant ), 
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
# z <- predict( foo$Model, y )
# a <- caret::confusionMatrix( z$predictions, y$BType4 )
# a$byClass
# 
# a$overall
# 
# round( diag(a$table/dim(y)[[1]]) / rowSums(z$table/dim(x)[[1]]), 3)
# 
# sum( diag(a$table/dim(y)[[1]]) ) # Accuracy calculation. NOTE dividing by total N. 

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
  
#  return( data.frame( 'TSSWtd' = TSS, 'TNRWtd' = Spec ))
return( data.frame( 'TNRWtd' = Spec ))
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
                    num.trees = RFC$ntree, replace = RFC$repl, importance = RFC$imp,
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
                    num.trees = RFC$ntree, replace = RFC$repl, importance = RFC$imp,
                    case.weights = wts[ x.train$BType4 ])    #Define case.weights
    
    #-- Evaluate the model using the test partition.
    results <- rbind( results, Results.Row( anRF, x.test )) 
    print(i)
    
  }
  return( results ) 
}



#============================================ DATA PREPARATION  =================================

#---------------------------------------------
# Partition input points according to regions.
# Region names and source features hard-coded in function
# Returns: List of point dataframes, one for each region
Partition.Test.Data <- function( pts ){
  
  out <- list()
  
  #-- Load the regions shape file containing the region pgons ... 
  pgons <- readOGR( file.path(source.dir, "/regions/BC_coast_regions.shp") )
  
  drop.list <- c('BT_Sorc','Rock','X','Y','ID')
  
  x <- pts@data
  x$BType4 <- as.factor( x$BType4 )
  #-- Return corrected coast  first ... 
  y <- x[ , !colnames(x) %in% drop.list ]
  names(y)[3:12] <- names.100m 
  out <- c( out, list( 'Coast' = y) )
  
  x <- pts[ pgons[ pgons$Region == "Haida Gwaii",], ]@data
  y <- x[ , !colnames(x) %in% drop.list ]
  names(y)[3:12] <- names.100m 
  out <- c( out, list( 'HG' = y) )
  
  x = pts[ pgons[ pgons$Region == "North Central Coast",], ]@data
  y <- x[ , !colnames(x) %in% drop.list ]
  names(y)[3:12] <- names.100m 
  out <- c( out, list( 'NCC' = y) )
  
  x = pts[ pgons[ pgons$Region == "West Coast Vancouver Island",], ]@data
  y <- x[ , !colnames(x) %in% drop.list ]
  names(y)[3:12] <- names.100m 
  out <- c( out, list( 'WCVI' = y) )
  
  x = pts[ pgons[ pgons$Region == "Queen Charlotte Strait",], ]@data
  y <- x[ , !colnames(x) %in% drop.list ]
  names(y)[3:12] <- names.100m 
  out <- c( out, list( 'QCS' = y) )
  
  x = pts[ pgons[ pgons$Region == "Strait of Georgia",], ]@data
  y <- x[ , !colnames(x) %in% drop.list ]
  names(y)[3:12] <- names.100m 
  out <- c( out, list( 'SOG' = y) )
  
  return( out )
}


#----------------------------
# Load 20 m (regional) predictor variables.
# Grabs the predictor values from the 20 m rasters. By region.
# Requires: list of points, fewer attributes is better! 
# Global var: bioregions 
Add.20m.Preds <- function( obsList ){
  
  outList <- vector("list", length = length(bioregions) )
  names( outList ) <- bioregions
  
  for (i in bioregions ) {
    rasters <- Load.Predictors( paste0( predictor.dir, "/", i ))
    preds  <- raster::extract(rasters, obsList, df = TRUE)
    # Drop ID as it duplicates source points  ... 
    preds <- preds[, !names(preds) %in% c('ID' ) ]
    merged <- cbind( as.data.frame(obsList), preds)
    
    merged <- merged[, !names(merged) %in% c('coords.x1', 'coords.x2' ) ]
    
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


#-----------------------------
# Load COASTWIDE predictors. 
Add.100m.Data <- function( obsList ){
  rasters <- Load.Predictors( predictors.coastwide )
  preds  <- raster::extract(rasters, obsList, df = TRUE)
  merged <- cbind( as.data.frame(obsList), preds)
  return(merged)
}


#----------------------------------------------------------------------------
# Straight up load of the independent point data sets from source shapefiles. 
# REQUIRES: source.dir to be set.  
# RETURNS: list of PointFeatures
Load.Point.Data <- function( shplist ) {
  
  a <- list()

  for (i in shplist) {
    a <- c( a, 
            list( 'dat' = readOGR( file.path(source.dir, i )) ))
  }
  
  names( a ) <- names( shplist )
  return( a )
}


#---------------------------------------------------------------------------
# Partition observations into training and testing when resampling the data. 
# train is assigned 1, testing remains 0. 
split.Obs.Data <- function( obs, seed, train.size ){
  
  set.seed( seed  )  # ensure repeatable results
  y <- nrow( obs ) # num of rows in the data.frame
  
  b <- sample( 1:y, round( train.size * y ), replace = FALSE )
  c <- replace( rep(0, y), b, 1)
  
  return( c )
}


#--------------------------------------------------------
# Pulls all the ID set data together for coastal tests ... 
# USES train.data as a global variable. 
Assemble.Coast.Test <- function( indat ){
  x <- rbind( indat[['NCC']], indat[['HG']], indat[['WCVI']], indat[['QCS']], indat[['SOG']] )
  return( x )
}


#--------------------------------------------------------------
# Builds a random forest model with internal partitioning of train/test data. 
# 2020/04/13: No longer used. Intent was to re-sample partitions but dropped from the paper. 
Rand.Ranger.Model <- function( x.data, x.formula, partition = 0.7, iterations = 1  ){
  
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
                       num.trees = RFC$ntree, replace = RFC$repl, importance = RFC$imp, oob.error = T,
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
                     num.trees = RFC$ntree, replace = RFC$repl, importance = RFC$imp, oob.error = T,
                     case.weights = wts[ x.data$BType4 ])    #Define case.weights
  
  return( list( 'Stats' = z, 'Model' = x.model ))
}

# # Older number-checking code for Load.Test.Data(). 
#   # for reading into ESRI File Geodatabase
#   subset(ogrDrivers(), grepl('GDB', name))
#   
#   # read the feature class specified by user
#   GIS.data <- readOGR(dsn = gdb, layer = fc)
#   
#   # Keep only COLUMNS specified from the GIS POINTS DATA FRAME ... 
#   observations <- GIS.data[, (names(GIS.data) %in% columns)] 
#   
#   # make X, Y columns (even if already exist to recreate spatial points dataframe later)
#   observations$X <- observations@coords[, 1]
#   observations$Y <- observations@coords[, 2]
#   
#   # number of observations before dropping invalid BType4 values
#   observations.before <- nrow(observations)
#  
#   # make sure that BType4 column has only values 1:4
#   observations <- observations[observations$BType4 %in% c(1, 2, 3, 4), ]
#   
#   # number of observations after dropping any records with invalid BType4 values
#   observations.after <- nrow(observations)
#   
#     # write a summary to the log file ... a BIT UNORTHODOX for a function perhaps. 
#   sink( file = file.path( output.dir, log.file ), append=TRUE)
#   
#   cat('\n\n\n***********************************************\n\n\n')
#   cat(paste0('\nNumber of observations before omitting records with invalid BType4 values:\n', observations.before))
#   cat(paste0('\n\nNumber of observations omitted due to records with invalid BType4 values:\n', observations.before - observations.after))
#   # percent of observations dropped
#   cat(paste0('\n\n', round(observations.after / observations.before * 100, 2), '% of observations have been kept after omitting records with invalid BType4 values'))
#   
#   # stop sink to log file
#   sink()

### fin.
