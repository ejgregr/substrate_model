#*******************************************************************************
# Script:  IDE_main.R
# Created: January 2020. EJG
# Updated: 2020/04/01 EJG (now detailed in README.md)
#
# This first script in the set: 
#   1: Loads all the observational data
#   2: Load predictors onto the observations
#   3: Build the necessary RF models
#   4: Compare the various models and IDE sets 

# Subsequent scripts: depth_effect.R, sample_effect.R, model_summaries.R.
# Supporting scripts: substrate_functions.R, Plot_Functions.R
#*******************************************************************************

rm(list=ls(all=T))  # Erase environment.

#-- Load necessary packages and functions ... 
source( "substrate_functions.R" )
source( "model_summaries.R" )
source( "Plot_Functions.R" )


#-------------------------------------------------------------------------
#-- PART 1: Load, prep, and inspect all the data. 
#-- Wrapper functions used to call other functions and keep it clean here. 

#-- Load existing observational data .. 
load( 'C:\\Users\\Edward\\Dropbox\\DFO Job\\Substrate2019\\Data\\loaded_data_2020-02-10-1246.RData')
# obs.data - training data
# dive.data, cam.data, ROV.data = Independent data sets
ls()

#-- OR RE-LOAD observational data  ... can take 20-30 min.
starttime <- Sys.time()

#-- Load all point observations from ArcGIS GDB. 
# Returns a list of the 4 point files.
obs.data <- Load.Observations()

names( obs.data )
head( obs.data$train )
table( obs.data$train$BT_Sorc )

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

dim( train.data.20m$HG )

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
x <- x[12:21]
names(train.data.100m)[8:17] <- x

#-- And flip bathy so that all are positive ... 
train.data.100m$bathy <- train.data.100m$bathy * -1

#-- Compare lengths to see how many points were lost ... 
Diff.Sets( obs.data$train, train.data.20m ); Diff.Sets( obs.data$dive, dive.data )
Diff.Sets( obs.data$cam, cam.data );         Diff.Sets( obs.data$ROV, ROV.data )


#-----------------------------
#-- Done loading source data. 
endtime <- Sys.time()
cat('---------------------------\n')
cat( 'Data build time: ', endtime - starttime, '\n\n' )


#-----------------------------
#-- Save out the populated observational data.
# Build a time stamp  ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( obs.data, dive.data, cam.data, ROV.data, train.data.100m, train.data.20m, 
      file = file.path(model.dir, paste0('loaded_data_', x, '.RData') ))


#==================================================
#-- PART 2: Build and evaluate the required RF models. 
# All models use wtd ranger(). 
# Models include:
#   1. Coastwide, 100 m
#   2. Regional, 20 m

#-- Load the latest RF models: Last update: 2020/03/21.
#   This is a BIG file ... 
load( 'C:\\Users\\Edward\\Dropbox\\DFO Job\\Substrate2019\\Models\\rf_allModels_2020-04-06-1409.RData' )

#  Also want the build statistics and results ... 
load( 'C:\\Users\\Edward\\Dropbox\\DFO Job\\Substrate2019\\Models\\buildResults_2020-04-06-1409.RData' )

ls()

#-- Or build them ... 
#=======================================================================================
#=== BUILD AND SUMMARIZE 6 RF MODELS 

starttime <- Sys.time()
set.seed( 42 )  # ensure repeatable results

tpart <- 0.7
iter  <- 1

#-- Need to save the performance stats in something ... 
#   Use a list of 3 results for each model.
#   Summary of variable importance currently only reports those with higher than median values

#-- Where the results live ... 
build.results <- list()
prev.table <- data.frame()

#----- Coastwide -----
x <- train.data.100m

foo <- Make.Ranger.Model( x, coast.formula, tpart, iter )
imp <- foo$Model$variable.importance

# Store results  ... 
rf.coast <- foo$Model
build.results <- c( build.results, 'Coast' = list( 'stats' = foo$Stats, 'Import' = imp) )
prev.table <- rbind( prev.table, 
                     cbind( 'Region' = 'Coast', foo$Prev))


#----- Region HG  -----
x <- train.data.20m$HG

foo <- Make.Ranger.Model( x, shore.formula, tpart, iter)
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.HG <- foo$Model
build.results <- c( build.results, 'HG' = list( 'stats' = foo$Stats, 'Import' = imp) )
prev.table <- rbind( prev.table, 
                      cbind( 'Region' = 'HG', foo$Prev))


#----- Region NCC  -----
x <- train.data.20m$NCC

foo <- Make.Ranger.Model( x, shore.formula, tpart, iter)
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.NCC <- foo$Model
build.results <- c( build.results, 'NCC' = list( 'stats' = foo$Stats, 'Import' = imp) )
prev.table <- rbind( prev.table, 
                     cbind( 'Region' = 'NCC', foo$Prev))


#----- Region WCVI  -----
x <- train.data.20m$WCVI

foo <- Make.Ranger.Model( x, shore.formula, tpart, iter)
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.WCVI <- foo$Model
build.results <- c( build.results, 'WCVI' = list( 'stats' = foo$Stats, 'Import' = imp) )
prev.table <- rbind( prev.table, 
                     cbind( 'Region' = 'WCVI', foo$Prev))


#----- Region QCS  -----
x <- train.data.20m$QCS

foo <- Make.Ranger.Model( x, shore.formula, tpart, iter)
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.QCS <- foo$Model
build.results <- c( build.results, 'QCS' = list( 'stats' = foo$Stats, 'Import' = imp) )
prev.table <- rbind( prev.table, 
                     cbind( 'Region' = 'QCS', foo$Prev))


#----- Region SOG  -----
x <- train.data.20m$SOG

foo <- Make.Ranger.Model( x, shore.formula, tpart, iter)
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.SOG <- foo$Model
build.results <- c( build.results, 'SOG' = list( 'stats' = foo$Stats, 'Import' = imp) )
prev.table <- rbind( prev.table, 
                     cbind( 'Region' = 'SOG', foo$Prev))

names(build.results)

#-----------------------------
#-- Done building ranger RF models. 

endtime <- Sys.time()
cat( 'Model build time: ', endtime - starttime, '\n\n' )
#last run for 1 iteration was 7.4 minutes.

#--------------------------------------------
#-- SAVE the resulting models. Takes minutes - its a big file. 

# Build a time stamp ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( rf.coast, rf.region.HG, rf.region.NCC, 
      rf.region.WCVI, rf.region.QCS, rf.region.SOG, 
      file = file.path( model.dir, paste0('rf_allModels_', x, '.RData')) )

save( build.results, prev.table,
      file = file.path( model.dir, paste0('buildResults_', x, '.RData')) )

#------------------------------------------------

#--- Create summaries of model builds in csv files ... for figure production.
# Build_results_Integrated.csv
# Build_results_byClassStats.csv
# Build_results_varImportance.csv
Summarize.Build( build.results )

#------- Produce and save some plots ----------------

#--- Heatmaps (tigures!) of build stats as png files from above csv files.
Build.Plots()

#--- Plot some facets of the Build results ... 

#--- Prevalence plots for obs and predicted for build testing partition faceted by Region
a <- Plot.Obs.Pred.Prevalence.Build( prev.table, pal.cole ) 
ggsave("Class Prevalence for Obs&Pred (Build) for Regions.png", a, dpi = 300, width = 16, height = 10, path = output.dir)

#-- User and producer accuracies by class ... 
a <- Plot.Class.Stats.For.Regions( 'Build_results_byClassStats.csv', pal.cole )
ggsave( 'Class Stats (build) for Regions.png', a, dpi = 300, width = 16, height = 10, path = output.dir)



#====================================================================================
#-- Independent Data Evaluation --
# Comparisons include:
#   1. Coastwide model vs. all 3 ID sets.
#     ID merged from the regions to test coastwide.
#     summarizes the merged IDS
#   2. Regional models vs. all 3 ID sets. 
#     summarizes the regional IDS

#-- Load it ... 
results.table <- read.csv( file.path(output.dir, 'IDS_evaluation_by_region.csv'))

#-- or build it ...

results.table <- Do.Independent.Evaluation()

#-- Fix the above build function BUT naming depends on existing data structures. ugh. 
results.table$Test.Data <- c('Dive','Cam','ROV','Dive', 'Cam','ROV','Dive','Cam','ROV',
                             'Dive','Cam','Dive','Cam','Dive','Cam')

out.file <- 'IDS_evaluation_by_region.csv'
write.csv( results.table, file = file.path(output.dir, out.file) )


#------------------------------------------------
#--- Plot some IDE results! ---

#-- User and producer accuracies by class ---
#     NEEDS A NEW DATA STRUCTURE, SO UPDATE TO THE SUMMARY CODE ... 
a <- Plot.Class.Stats.For.Regions( 'IDS_evaluation_by_region.csv', pal.cole )
ggsave( 'Class Stats (IDE) for Regions.png', a, dpi = 300, width = 16, height = 10, path = output.dir)


#--- Integrated statistics for each IDS by Region ---
# Use a predetermined collection of stats (see function)
a <- Plot.Stat.By.IDS.For.Regions( 'IDS_evaluation_by_region.csv', pal.cole )
ggsave( paste0( 'Integrated Stats by IDE for Regions.png'), a, dpi = 300, width = 16, height = 10, path = output.dir)






#-- FIN.
