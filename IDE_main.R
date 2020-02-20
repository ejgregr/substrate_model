#*******************************************************************************
# Script:  IDE_main.R
# Purpose: Evaluate wtd RF models with Independent data (ID) across various data partitions
#   Part 1: Load all the observational data
#   Part 2: Build the necessary RF models
#   Part 3: Compare the various models and ID sets 
# Created: January 2020. EJG
# Updated: 2020/20/02. EJG
#  - Moved and renamed fully operational version to GitHub. 
#  - Main objective here is bulding the necessary loops to do all the tests. 
#           Feb 2020
#  - All data inspection and comparison now in IDE_data_compare_V1.R script.
#  - Entire script runs to completion, yielding data table and report (assuming it can find the data).
#  - First upload to GitLab
# Notes:
#  - Question as to whether CI is required, or if kappa SE is sufficient
#  - Use ConfusionMatrix() for consistency of statistics. 
#     Always order parameters as (predicted, observed)
#  - ranger has inbag.counts - do we need to calculate our own OOB error?
#  - Adjust SIGNS on the 100 m and 20 m bathymetries if comparing models.
#       Flipped sign on 20 m during data loading.
#*******************************************************************************

rm(list=ls(all=T))  #Erase all previously saved objects

#-- Load necessary packages and functions ... 
source( "substrate_functions.R" )


#-------------------------------------------------------------------------
#-- PART 1: Load, prep, and inspect all the data. 
#-- Wrapper functions used to call other functions and keep it clean here. 

#-- Load existing observational data .. 
load( 'C:\\Users\\Edward\\Dropbox\\DFO Job\\Substrate2019\\Data\\loaded_data_2020-02-10-1246.RData')
ls()

#-- OR Build new ... can take 20-30 min.
starttime <- Sys.time()

#-- Load all point observations from ArcGIS GDB. 
# Returns a list of the 4 point files.
obs.data <- Load.Observations()

#-- Attach predictor data to the observations. Both resolutions
#   NOTE For the obs data already have the predictors attached, so just take the data. 
train.data.100m <- obs.data$train@data

# Now use the spatial points to pull predictor data ... 
#   First part splits the Train Obs into the bioregions and replaces the predictors with 20 m.
#   Remove the 100 m data beforehand ... 
#   takes ~10 minutes because of the extents/size of the coastwide Obs training dataset.
#   Necessary to compare models across resolutions.
foo <- obs.data$train[ , c('BType4', 'TestDat', 'ID') ]
train.data.20m <- Add.20m.Data( foo )

#-- NOTE there are 977 non-unique points in train.data.20m because of overlap between 
# some of the regions (mainly HG and NCC)

#-- Load predictors on the remaining observations (the ID sets). 
#   NB: Different set of obs attributes retained here. 
dive.data  <- Add.20m.Data( obs.data$dive )
cam.data   <- Add.20m.Data( obs.data$cam )
ROV.data   <- Add.20m.Data( obs.data$ROV )

#-- Minor adustments to the 100 m predictor data.
#   First change predictor names to match the 20 m so we can use a single RF formula
#   Cross-feference the names you need, removing fetch (not used in coastal model).
#   REVIEW if source files change. 
x <- names( dive.data$NCC )
x <- x[ !x %in% c("fetch")]
x <- x[12:20]
names(train.data.100m)[9:17] <- x

#-- And flip bathy so that all are positive ... 
train.data.100m$bathy <- train.data.100m$bathy * -1


#-- Compare lengths to see how many points were lost ... 
Diff.Sets( obs.data$train, train.data.20m ); Diff.Sets( obs.data$dive, dive.data )
Diff.Sets( obs.data$cam, cam.data );         Diff.Sets( obs.data$ROV, ROV.data )


#-----------------------------
#-- Done loadings source data. 
endtime <- Sys.time()
cat('---------------------------\n')
cat( 'Data build time: ', endtime - starttime, '\n\n' )

#-- Save out the populated observational data.
# Build a time stamp  ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( obs.data, dive.data, cam.data, ROV.data, train.data.100m, train.data.20m, 
      file = file.path(model.dir, paste0('loaded_data_', x, '.RData') ))


#==================================================
#-- PART 2: Build and save the required RF models. 
# All models use wtd ranger() unless otherwise noted. 
# Models include:
#   1. Coastwide, 100 m, weighted and unweighted.
#   2. Regional, 20 m, using training data

#-- Load the latest RF models: Last update: 2020/02/02.
load( 'C:\\Users\\Edward\\Dropbox\\DFO Job\\Substrate2019\\Models\\rf_allModels_2020-02-11-1059.RData')
ls()

#-- Or build them ... 
starttime <- Sys.time()

#-----------------------------------------
# Coastwide, weighted model using ranger()

# The weighting vector ...
props <- train.data.100m %>% group_by(BType4) %>% count(BType4) 
#prop.Btype4 <- obs@data %>% group_by(BType4) %>% count(BType4) 
wts <- 1 - ( props$n / sum( props$n ))

rf.coast.wtd  <- ranger( coast.formula,
                       data = train.data.100m,
                       num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                       case.weights = wts[ train.data.100m$BType4 ])    #Define case.weights

#---------------------------------------
# Regional models, weighted, w 20 m data.
# Using all the Training data with the 20 m predictors mapped onto it. 
names(train.data.20m)
str(train.data.20m$NCC)

#-- HG ..
x <- train.data.20m$HG
props <- x %>% group_by(BType4) %>% count(BType4) 
wts <- 1 - ( props$n / sum( props$n ))
rf.region.HG <- ranger(  shore.formula,
                         data = x,
                         num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                         case.weights = wts[ x$BType4 ])    #Define case.weights

#-- NCC ..
x <- train.data.20m$NCC
props <- x %>% group_by(BType4) %>% count(BType4) 
wts <- 1 - ( props$n / sum( props$n ))
rf.region.NCC <- ranger( shore.formula,
                        data = x,
                        num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                        case.weights = wts[ x$BType4 ])    #Define case.weights

#-- WCVI ..
x <- train.data.20m$WCVI
props <- x %>% group_by(BType4) %>% count(BType4) 
wts <- 1 - ( props$n / sum( props$n ))
rf.region.WCVI <- ranger( shore.formula,
                        data = x,
                        num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                        case.weights = wts[ x$BType4 ])    #Define case.weights

#-- QCS ..
x <- train.data.20m$QCS
props <- x %>% group_by(BType4) %>% count(BType4) 
wts <- 1 - ( props$n / sum( props$n ))
rf.region.QCS <- ranger( shore.formula,
                        data = x,
                        num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                        case.weights = wts[ x$BType4 ])    #Define case.weights

#-- SOG ..
x <- train.data.20m$SOG
props <- x %>% group_by(BType4) %>% count(BType4) 
wts <- 1 - ( props$n / sum( props$n ))
rf.region.SOG <- ranger( shore.formula,
                        data = x,
                        num.trees = ntree, replace = repl, importance = imp.2, oob.error = T,
                        case.weights = wts[ x$BType4 ])    #Define case.weights

#-----------------------------
#-- Done loadings source data. 
endtime <- Sys.time()
cat('---------------------------\n')
cat( 'Model build time: ', endtime - starttime, '\n\n' )


#-- SAVE The resulting models. Takes minutes - its a big file. 

# Build a time stamp ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( rf.coast.wtd, rf.region.HG, rf.region.NCC, 
      rf.region.WCVI, rf.region.QCS, rf.region.SOG, 
      file = file.path( model.dir, paste0('rf_allModels_', x, '.RData')) )


#=====================================================================
#-- Do some Independent Data Evaluation --
# Comparisions include:
#   1. Coastwide model vs. all 3 ID sets.
#     ID merged from the regions to test coastwide.
#     summarizes the merged IDS
#   2. Regional models vs. all 3 ID sets. 
#     summarizes the regional IDS

log.file <- 'independent_data_evaluation_log.txt'

sink( file = file.path(output.dir, log.file))
cat("Independent Data Evaluation \n")
cat("Run Date: ", date(), "\n")
cat("---------------------------\n\n")

#-- Summarize the different observational data sets.
cat("Summary of the observational data\n")
cat("---------------------------------\n\n")

cat( '100 m train data - Summary\n')
x.test <- train.data.100m
x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
cat( 'N =         ', nrow(x.test), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')


#- Build a table to hold the performance scores for the various comparisons.
# All models built with full set of training observations.

results.table <- NULL

#-- Coast model vs. 100 m train ...
compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = '100 m Train' )
x <- cbind( compare.what, Results.Row( rf.coast.wtd, x.test ))
results.table <- rbind( results.table, x ) 


#-- Coast model vs. each ID set ... 

# Dive Data - pull the regional data together  ... 
rm( 'x', 'x.test', 'x.pvec')
x.test <- Assemble.Coast.Test( dive.data )

cat("\n---------------------------------\n\n")
cat("20 m Dive IDS - Summary\n")
x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
cat( 'N =         ', nrow(x.test), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')

#-- Build the table piece ... 
compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = '20 m Dive' )
x <- cbind( compare.what, Results.Row( rf.coast.wtd, x.test ))
results.table <- rbind( results.table, x ) 


# Repeat for the Cam IDS  ... 
rm( 'x', 'x.test', 'x.pvec')
x.test <- Assemble.Coast.Test( cam.data )

cat("\n---------------------------------\n\n")
cat("20 m Camera IDS - Summary\n")
x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
cat( 'N =         ', nrow(x.test), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')

#-- Build the table piece ... 
compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = '20 m Cam' )
x <- cbind( compare.what, Results.Row( rf.coast.wtd, x.test ))
results.table <- rbind( results.table, x ) 


# Repeat for the ROV IDS  ... 
rm( 'x', 'x.test', 'x.pvec')
x.test <- Assemble.Coast.Test( ROV.data )

cat("\n---------------------------------\n\n")
cat("20 m ROV IDS - Summary\n")
x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
cat( 'N =         ', nrow(x.test), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')

#-- Build the table piece ... 
compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = '20 m ROV' )
x <- cbind( compare.what, Results.Row( rf.coast.wtd, x.test ))
results.table <- rbind( results.table, x ) 

# cat( results.table )

cat("\n\n---------------------------------\n")
cat("---- IDS evaluation Coastwide ---------\n\n")

results.table

sink()


#------------------------------------------------------------
#-- Regional, each by all 3 ID sets.
#   THIS IS 5 x 3 = 15 comparisions ... with ROV exception.

#-- Requres: Models to be loaded (rf.region.XX)
  
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
      # good to go ... 
        compare.what <- data.frame( 'Model' = i, 'Test Data' = j )
      
        # grab the IDS data ... 
        x <- eval(parse( text=paste0( j, ".data" ) ))
        x.test <- x[[ i ]]
        # make the results ... 
        x <- cbind( compare.what, Results.Row( m.model, x.test ))
        results.table <- rbind( results.table, x ) 
        }
  }
}

sink( file = file.path(output.dir, log.file), append=TRUE)
cat("\n\n----------------------------------------------\n")
cat("---------- IDS evaluation by Region ----------\n")
cat("    Regional models built w 20 m data \n\n")
results.table
sink()




#----------------------------
# Now check for depth effect. 








#-- FIN.
