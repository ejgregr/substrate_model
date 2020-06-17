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
#  - 2020/05/10: All unused code removed. Includes Rand.Ranger.Model and its associated functions
# like Summary.Row(), Wtd.Stats, ...

#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
#-- Load packages using neat loading code from Cole ... 

# check for any required packages that aren't installed and install them
required.packages <- c("superheat", "ggplot2", "reshape2", "tidyr","dplyr", 
                       "diffeR", "vegan", "ranger", "rgdal", "raster", "stringr",
                       "measures", "caret", "PresenceAbsence" )
# Currently not using:  "e1071", "randomForest", "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr",
#                       "tinytex", "rmarkdown", "binr"

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)

#lapply(required.packages, library, character.only = TRUE)

#-- devtools required above to properly install colours ... the first time?
#devtools::install_github("jakelawlor/PNWColors") 
library(PNWColors)

#----------------------------------------------------------------------------
#-- Data sources and other constants ...  

#-- As per Cole's reporting code, set an output directory and initiate the log file.
#       Directory will be created if doesn't exist; file will be overwritten if it does.
results.dir   <- 'C:/Users/Edward/Dropbox/DFO Job/Substrate2019/Results'
model.dir     <- 'C:/Users/Edward/Dropbox/DFO Job/Substrate2019/Models'
predictor.dir <- 'C:/Data/SpaceData/Substrate2019/Predictors'
source.dir    <- 'C:/Data/SpaceData/Substrate2019'
raster.dir    <- 'C:/Users/Edward/Dropbox/DFO Job/Substrate2019/Rasters'

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
# 2020/04/13: Updated to surface train and testing subsamples. 
# 2020/05/07: Dropped the randomization code earlier (hence the name Fixed). Code for the 
#             Randomize version is in the pre- May 10 commit. 
# Returns:  a list with the final fitted model and a summary table of performance statistics.
# Requires: training and testing data sets, a formula. 
# Calls:    Results.Row()
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


#-----------------------------------------------------------------------------------------------
# Calculate statistics for a given RF model (testModel) and test set of observations (testData). 
# 2020/05/10: Updated to include a subset of potential statistics, as approved by project team. 
# Returns: a data frame with a single row. suitable for rbind().
# Requires: testModel: a weighted RF model built with ranger()
#           testData : dataframe with correct predictors and observed BType
#           testData must have the same predictors attached as those used to build testModel.
# Calls:  diffr.Stats(); Prev.Balance()
Results.Row <- function( testModel, testData, mant = 3 ){

  x <- predict( testModel, testData )
  y <- caret::confusionMatrix( x$predictions, testData$BType4 )
  z <- data.frame(y$byClass)
  prev <- Prev.Balance( summary(testData$BType4) / sum( summary(testData$BType4)) )
  
  out1 <- data.frame( 
    "N"         = nrow( testData),
    "Imbalance" = round( prev, mant ), 
    "OOB"       = round( testModel$prediction.error, mant ),
    "TSS"       = TSS.Calc( y$table ),
#    "Kappa"     = round( as.numeric( y$overall[ 'Kappa' ]), mant ),
    "Accuracy"  = round( as.numeric( y$overall[ 'Accuracy' ]), mant ),
    'TNRWtd'    = round( sum( z$Specificity * z$Prevalence ), mant ),
    round( diffr.Stats( y$table )/nrow(testData), mant ) #returns some pre-named, Pontius stuff.
  )
  out2 <- data.frame( 
    "User"  = round( diag(y$table) / rowSums(y$table), mant ), 
    "Prod"  = round( diag(y$table) / colSums(y$table), mant ),
    "PrevObs"  = round( as.vector( table( testData$BType4 )) ),
    "PrevPred" = round( as.vector( table( x$predictions )) ))
  
  return( list( 'Integrated' = out1, 'PerClass' = t(out2) ))
}


# w <- obs.20mIV$HG[ obs.20mIV$HG$TestDat ==1, ]
# Results.Row( rf.region.HG, w )
# 
# x <- predict( rf.region.HG, w )
# y <- caret::confusionMatrix( x$predictions, w$BType4 )
# y$byClass
# y$overall
# 
# z <- diffr.Stats( y$table )
# names(z)


#----------------------------------------------------------------
#-- Calculate the True Skill Statistic for a 4x4 matrix following
#   Allouche et al. 2006. Applies 2 steps - first calculates the 
#   Expected correct classifications expected by chance; 
#   then substracts this from the diagonal values, and reports the 
#   sum of the diagonal values. 
# Returns: A single TSS statistic for the entire table.
# Takes:   a caret contingency table
# Assumes: 4x4 matrix
# Notes:   Correcting for chance success reduces the sample size on 
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
# Returns: A one-row data.frame of named statistics
# Takes:   a caret contingency table
# Assumes: 4x4 matrix
# Notes:   This will need to be updated based on detailed review of the 
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
#-- Calculate the imbalance in prevalence.
# Returns: Generalized statistic to describe deviance from prevalence.
#   0 is fully balanced. Max = N/number of categories
# Takes: list of prevalence values
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


#----------------------------------
# Function to predict a raster surface, write to disk, and export a png of the surface
# env.predictors: raster stack of environmental predictors, ranger.model: model object, output.directory: path to output
Predict.Surface <- function(env.predictors, ranger.model, output.directory, nm, pngpal){
  
  # create output directory - this creates it in the working directory. If it exists, do not show warning
#  dir.create(output.directory, showWarnings = FALSE)
  
  # save ranger model to disk as RData file
#  save(ranger.model, file = file.path(output.directory, 'ranger_model.RData'))
  
  # predict substrate using raster stack and model file
  raster.obj <- raster::predict(object   = env.predictors,  
                                model    = ranger.model,               
                                progress = 'text',
                                fun = function(model, ...) predict(model, ...)$predictions)
  
  # write raster file to disk
  writeRaster(raster.obj, file.path(output.directory, 
                                    paste0(nm, '_classified_substrate.tif')), 
                                    format = 'GTiff', datatype = 'INT2S', overwrite = TRUE)
  
  # generate table of proportions for each class in predicted raster
  raster.prop <- round(100 * prop.table(table(na.omit(as.data.frame(getValues(raster.obj))))), 2)
  
  # colour palette for map
  #pal <- c("#999999", "#33bbff", "#ffff99", "#ffb84d")
  
  # labels for legend
  labels <- c("Rock", "Mixed", "Sand", "Mud")
  
  # Map (up to 5,000,000 pixels)
  png(file=file.path(output.directory, paste0(nm, "substrate_raster.png")) ,
      height = 7, width = 6, units = "in", res = 400)
  plot(raster.obj, maxpixels=5000000, col=pngpal, legend=FALSE,
       xlab = "Easting", ylab = "Northing", cex.axis = .5, cex.lab = .75)
  legend(x = "bottomleft",
         legend = labels, fill = pngpal, title=NA, bg = NA, box.col = NA)
  
  dev.off()
  
  return(list(raster.obj, raster.prop))
  
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
# Assumes: obs class attribute = BType4 with 4 classes.
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


#===============================================================================
#-- Independent Data Evaluation
#
# 2020/05/10: Needs updating based on objectives.
#
# Comparisons include:
#   1. Coastwide model vs. all 3 ID sets.
#     ID merged from the regions to test coastwide.
#     summarizes the merged IDS
#   2. Regional models vs. all 3 ID sets. 
#     summarizes the regional IDS
# Build a table to hold the performance scores for the various comparisons.
# Requires: loaded independent data (dive, cam, ROV)
#           all RF models built.
Do.Independent.Evaluation <- function( results=NULL ){
  #-- Part 1: Coast model vs. each ID set. Regional IDE sets need to be assembled.  
  results.int <- NULL
  results.class <- NULL
  
  # Dive Data - pull the regional data together  ... 
  x.test <- Assemble.Coast.Test( dive.20mIV )
  
  #-- Build the table piece ... 
  compare.what <- data.frame( 'Region' = 'Coast', 'IDS' = 'Dive' )
  w <- Results.Row( rf.region.Coast, x.test )
  x <- cbind( compare.what, w$Integrated )
  results.int <- rbind( results.int, x ) 

  y <- cbind( compare.what, 'Stat' = row.names( w$PerClass ), w$PerClass )
  results.class <- rbind( results.class, y )
  print( 'Dive test done ...')
  
  # Repeat for the Cam IDS  ... 
  x.test <- Assemble.Coast.Test( cam.20mIV )
  
  #-- Build the table piece ... 
  compare.what <- data.frame( 'Region' = 'Coast', 'IDS' = 'Cam' )
  w <- Results.Row( rf.region.Coast, x.test )
  x <- cbind( compare.what, w$Integrated )
  results.int <- rbind( results.int, x ) 
  
  y <- cbind( compare.what, 'Stat' = row.names( w$PerClass ), w$PerClass )
  results.class <- rbind( results.class, y )
  print( 'Cam test done ...')
  
    # Repeat for the ROV IDS  ... 
  x.test <- Assemble.Coast.Test( ROV.20mIV )
  
  #-- Build the table piece ... 
  compare.what <- data.frame( 'Region' = 'Coast', 'IDS' = 'ROV' )
  w <- Results.Row( rf.region.Coast, x.test )
  x <- cbind( compare.what, w$Integrated )
  results.int <- rbind( results.int, x ) 
  
  y <- cbind( compare.what, 'Stat' = row.names( w$PerClass ), w$PerClass )
  results.class <- rbind( results.class, y )
  print( 'ROV test done ...')
  
  #------------------------------------------------------------
  #-- Part 2: Regional runs, each with all 3 ID sets.
  #   THIS IS 5 x 3 = 15 comparisions ... with ROV exception.
  #   Requres: Models to be loaded (rf.region.XX)
  
  IDS <- c("dive", "cam", "ROV")
  
  for (i in bioregions) {
    # select the regional model
    m.model <- eval(parse( text=paste0( "rf.region.", i ) ))
    cat(i, "\n")
    
    for (j in IDS ){
      cat("  ",j, "\n")  
      
      if (j == "ROV" && i %in% c('WCVI', 'QCS', 'SOG') ){
        # Bail ... 
        print( "skipping ... ")
      } else {
        # good to go ... except need to upper case first letter of IDS ...  
        jj <- sapply(j, function(x) 
          paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="") )  # toupper for first character
        compare.what <- data.frame( 'Region' = i, 'IDS' = jj )
        
        # grab the IDS data ... 
        x <- eval(parse( text=paste0( j, ".20mIV" ) ))
        x.test <- x[[ i ]]
        # make the results ... 
        w <- Results.Row( m.model, x.test )
        x <- cbind( compare.what, w$Integrated )
        results.int <- rbind( results.int, x ) 
        
        y <- cbind( compare.what, 'Stat' = row.names( w$PerClass ), w$PerClass )
        results.class <- rbind( results.class, y )
        
      }
    }
  }
  rownames( results ) <- NULL
  return( list( 'Integrated' = results.int, 'PerClass' = results.class ))
}


#============================================ DATA PREPARATION  =================================

#-------------------------------------------------------------------------------
# Construct summary tables of build results. 
# Takes: The build.results table from building all the RF models. 
# Returns: Nothing. Side effect is to write results to CSV files:
#     Build_results_Integrated.csv
#     Build_results_byClassStats.csv
#     Build_results_test_ClassPrevalence.csv
#     Build_results_varImportance.csv
Summarize.Build <- function( build.df ){
  # Bunch of hacking here to pull the tables together ... 
  # Requires: build.results data.frame as input.
  #   Includes pairs of lists for each run. 
  #   First is list of Integrated and perClass stats; Second is list of variable importance
  # ASSUMES 6 regions done.
  
  out <- list()
  # Build a list of names ... 
  a <- names( build.df )[c(1,3,5,7,9,11)]
  nm <- unlist( lapply( a, strsplit, split = '\\.' ))[c(1,3,5,7,9,11)]
  
  
  # Pull Integrated Stats ... number are the rows for each region-specific stat
  x <- do.call(rbind.data.frame, 
               build.df[ c(1,3,5,7,9,11) ] )
  
  # x has 2 components, the Integrated bit, .... 
  y <- do.call( rbind.data.frame, x$Integrated )
  
  z <- cbind( 'Region' = nm, y )
  row.names(z) <- NULL
  
  out.file <- 'Build_results_Integrated.csv'
  write.csv( z, file = file.path(results.dir, out.file) )
  out[[ 'build.results.Integrated' ]] <- z
  
  # and the PerClass bit which includes:
  #   1) A table with both User and Producer. A df with 2 sets of rows.
  y <- do.call( rbind,  x$PerClass )
  z <- y[ (row.names(y) == 'User') | (row.names(y) == 'Prod') , ]
  
  # User first ... 
  y.usr <- data.frame( 'Region' = nm, z[ c(1,3,5,7,9,11), ])
  row.names(y.usr) <- NULL
  
  # Now producer  
  y.prd <- data.frame( 'Region' = nm, z[ c(2,4,6,8,10,12), ])
  row.names(y.prd) <- NULL
  
  zz <- rbind(
    cbind( 'Stat' = 'User', y.usr ),
    cbind( 'Stat' = 'Prod', y.prd ) )
  colnames(zz) <- c('Stat','Region','Hard','Mixed','Sand','Mud')
  out.file <- 'Build_results_byClassStats.csv'
  write.csv( zz, file = file.path(results.dir, out.file) )
  out[[ 'build.results.ByClass' ]] <- zz
  
  #   2) The prevalence of the testing obs vs. predicted.
  z <- y[ (row.names(y) == 'PrevObs') | (row.names(y) == 'PrevPred') , ]
  
  # Obs first ... 
  y.obs <- data.frame( 'Region' = nm, z[ c(1,3,5,7,9,11), ])
  row.names(y.obs) <- NULL
  
  # Now predicted  
  y.prd <- data.frame( 'Region' = nm, z[ c(2,4,6,8,10,12), ])
  row.names(y.prd) <- NULL
  
  zz <- rbind(
    cbind( 'Stat' = 'Obs', y.obs ),
    cbind( 'Stat' = 'Pred', y.prd ) )
  colnames(zz) <- c('Stat','Region','Hard','Mixed','Sand','Mud')
  out.file <- 'Build_results_test_ClassPrevalence.csv'
  write.csv( zz, file = file.path(results.dir, out.file) )
  out[[ 'build.results.ClassPrev' ]] <- zz
  
  #-- And finally variable importance ... 
  # Integrated stats first ... 
  # 2020/04/09: Moved the ranking to the plot function so values can go thru and be displayed 
  
  y <- data.frame(
    rbind( as.vector( c( build.df$Coast.Import, 0 )),
           as.vector( build.df$HG.Import    ),
           as.vector( build.df$NCC.Import   ),
           as.vector( build.df$WCVI.Import  ),
           as.vector( build.df$QCS.Import   ),
           as.vector( build.df$SOG.Import   )
    ))
  
  row.names(y) <- nm
  colnames(y)  <- names( build.df$HG.Import )
  
  out.file <- 'Build_results_varImportance.csv'
  write.csv( y, file = file.path(results.dir, out.file) )
  out[[ 'build.results.VarImport' ]] <- y
  
  return( out )
}


#---------------------------------------------
# Partition input points according to regions.
# Returns: List of point dataframes, one for each region
# Requires: names.100m as a global (defined during data load)
# Notes:   Region names and source features hard-coded in function
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


#---------------------------------------------
# Partition input points into low and high density piles
# Returns: List of 2 dataframes, one for each of low/hi density regions
# Requires: names.100m as a global (defined during data load)
#           Obs points, and the density shape file. 
Partition.By.Density <- function( pts ){
  
  out <- list()
  
  #-- Load the regions shape file containing the High Density pgons ... 
  pgons <- readOGR( file.path(source.dir, "/regions/hi_density_area.shp") )
  
  drop.list <- c('BT_Sorc','Rock','X','Y','ID')
  
  
  # Split the points into those that fall in the High density area ... 
  dens.pts <- x[ pgons[ pgons$name == "High_Density_Area",], ]@data
  # and those that don't ... (confirmed IDs are unique)
  sparse.pts <- x[ !(x$ID %in% dens.pts$ID), ]@data
  
  # Now fix the attribute lists 
  dens.pts$BType4   <- as.factor( dens.pts$BType4 )
  sparse.pts$BType4 <- as.factor( sparse.pts$BType4 )
  
  x <- dens.pts[ , !colnames(dens.pts) %in% drop.list ]
  names(x)[3:12] <- names.100m 
  out <- c( out, list( 'Dense' = x) )
  
  x <- sparse.pts[ , !colnames(sparse.pts) %in% drop.list ]
  names(x)[3:12] <- names.100m 
  out <- c( out, list( 'Sparse' = x) )
  
  
  return( out )
}


#----------------------------
# Loads 20 m (regional) predictor data onto observation points.
# Returns:  List of regions ea containing list of points with predictors attached.
# Takes:    List of points, fewer attributes is better! 
# Requires: Global var: bioregions 
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
    merged$BType4 <- as.factor( merged$BType4 )  #Ensure this is stored as a factor ... 
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
# Straight up load of the independent point data sets from source shapefiles. 
# Returns: list of PointFeatures
# Requires: source.dir to be set.  
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


### fin.
