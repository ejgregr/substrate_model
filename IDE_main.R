#*******************************************************************************
# Script:  IDE_main.R
# Created: January 2020. EJG
# Updated: 2020/04/09 EJG (updates now detailed in README.md)
#
# This first script in the set:  
#   1: Loads all the observational data
#   2: Load predictors onto the observations
#   3: Build the necessary RF models
#   4: Compare the various models and IDE sets 

# Subsequent scripts: depth_effect.R, sample_effect.R, model_summaries.R.
# Supporting scripts: substrate_functions.R, Plot_Functions.R

# NOTES:
#   2020/04/09: sample_effect.R DEFERRED UNTIL SUBSTRATE PAPER COMPLETE
#
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
load( file.path( model.dir, 'loaded_data_2020-05-04-1426.RData' ))
# point.data - training data
# dive.data, cam.data, ROV.data = Independent data sets
ls()

#-- OR RE-LOAD observational data  ... can take 20-30 min.
starttime <- Sys.time()

#-- Load all point observations from ArcGIS GDB. 
# Returns a list of the 4 point files, with attributes unchanged.
point.data <- Load.Point.Data( source.files )

names( point.data )
str( point.data$Dive )

head( point.data$Obs )
dim( point.data$Obs )

#-- Check the train/test partition. Modify if desired ... 
#obs.100m@data$TestDat <- split.Obs.Data( obs.100m, 42, 0.7 )
sum( point.data$Obs$TestDat == 0 ) / nrow( point.data$Obs )


#-- Adjust a couple things in the Obs data ... 
# Ensure BType4 is a factor ... 
point.data$Obs$BType4 <- factor( point.data$Obs$BType4  )

# Flip bathy so that all are positive; and ensure BType4 is a factor. 
point.data$Obs$bathy <- point.data$Obs$bathy * -1

# check ... 
str( point.data$Obs )

#-------- Now associate the predictor data (IVs) with the points ------

#-- Part 1) Prepare the 20 m training/testing data (obs). 
# These are divided into bioregions, cuz the 20 m predictors depend on the 
# 20 m regional bathymetries. 

#-- Drop NRCan data since these were added to improve performance of the 100 m model at depth, 
#   and were determined to not be of sufficient accuracy at 20 m. 

table( point.data$Obs$BT_Sorc )

#-- For the 20 m models, keep only the NRCan-free points ... 
a <- point.data$Obs[ point.data$Obs[['BT_Sorc']] != 'NRCan_Observations',] 
dim(a)
names( a )
str( a )

#-- use the points in 'a' to pull the 20 m predictor data. 
#   Divides the observation points into bioregions, then pulls predictors from 20 m raster stack
#   Keep only the relevant attributes (i.e., no 100 m IVs or unnecessary fields.
#   Takes ~10 minutes because of the extents/size of the coastwide Obs training dataset.

a <- a[ , c('ID', 'BType4', 'TestDat') ]
obs.20mIV <- Add.20m.Preds( a )
rm( a )

# -- NOTE there are 977 non-unique points in obs.20mIV - these are duplicates, 
# created because of of overlap between the regions (mainly HG and NCC).

#-- Part 2) Load predictors onto the remaining observations (the ID sets). 
#   Standardize attributes on the IDS, and also divide by bioregions. 
#   BType4 turned into a factor in Add.20m.Preds() by data.frame()

a <- point.data$Dive[ , c('ID', 'RMSM', 'bathy_m' )]
names( a ) <- c('ID', 'BType4', 'bathy_m')
dive.20mIV  <- Add.20m.Preds( a )

a <- point.data$Cam[ , c('cellID', 'RMSM', 'maxDepth_m' )]
names( a ) <- c('ID', 'BType4', 'bathy_m')
cam.20mIV   <- Add.20m.Preds( a )

a <- point.data$ROV[ , c('cellID', 'RMSM', 'maxDepth_m' )]
names( a ) <- c('ID', 'BType4', 'bathy_m')
ROV.20mIV   <- Add.20m.Preds( a )
rm(a)

names( ROV.20mIV$HG )
str( ROV.20mIV$HG )

# Assess the success of assigning 20 m IVs to the point observations.
# Compare number of records with and w/out the 20 m IVs attached to check how many points
# could not be assigned 20 m IVs: 

Diff.Sets( point.data$Obs, obs.20mIV ) # All the Obs data 
Diff.Sets( point.data$Dive, dive.20mIV )    # dive ... etc. 
Diff.Sets( point.data$Cam, cam.20mIV )     
Diff.Sets( point.data$ROV, ROV.20mIV )

# Obs losses likely due to circulation and tidal models missing at 20 m (not checked)
# Extra point picked up by Dive and Cam data because of overlapping regional bathymetries. Not a problem for analysis. 

#-- Part 3) Prepare the train/testing Obs data with 100 m predictors. 
obs.100mIV <- point.data$Obs@data

#-- Make some necessary adustments ...

# Drop unused attributes.
obs.100mIV <- obs.100mIV[ , !names(obs.100mIV) %in% c('Rock', 'X', 'Y') ]

# Rename 100 m predictors to match names used for 20 m predictors 
x <- names( dive.20mIV$NCC )
x <- x[ !x %in% c("fetch")]
names.100m <- x[4:13]
names(obs.100mIV)[5:14] <- names.100m


#--------- Done loading source data --------------
endtime <- Sys.time()
cat('---------------------------\n')
cat( 'Data build time: ', endtime - starttime, '\n\n' )

#-----------------------------
#-- Save out the populated observational data.
# Build a time stamp  ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( point.data, dive.20mIV, cam.20mIV, ROV.20mIV, 
      obs.100mIV, obs.20mIV, 
      file = file.path(model.dir, paste0('loaded_data_', x, '.RData') ))


#=================================================================================
#-- PART 2: Build and evaluate the required RF models. 
# All models use wtd ranger(). 
# Models include:
#   1. Coastwide, 100 m
#   2. Regional, 20 m

#-- Load the latest RF models (This is a BIG file) ... 
load( file.path( model.dir, 'rf_allModels_2020-04-06-1701.RData' ))

#  Also want the build statistics and results ... 
load( file.path( model.dir, 'buildResults_2020-04-06-1701.RData' ))

ls()

#-- Or build them ... 
#===================================
#=== BUILD AND SUMMARIZE 6 RF MODELS 
#   Previoulsly built using Rand.Ranger.Model but simplified now to surface train/test data.
#   Fixed.Ranger.Model() also using only the train data to parameterize. 

starttime <- Sys.time()

#-- Need to save the performance stats in something ... 
#   Use a list of 3 results for each model.
#   Summary of variable importance currently only reports those with higher than median values

#-- Where the results live ... 
build.results <- list()

#----- Coastwide -----
x <- obs.100mIV[ obs.100mIV$TestDat == 0, ] #train
y <- obs.100mIV[ obs.100mIV$TestDat == 1, ] #test

foo <- Fixed.Ranger.Model( x, y, coast.formula)
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.Coast <- foo$Model
build.results <- c( build.results, 'Coast' = list( 'stats' = foo$Stats, 'Import' = imp) )


#----- Region HG  -----
a <- obs.20mIV$HG
x <- a[ a$TestDat == 0, ]
y <- a[ a$TestDat == 1, ]

foo <- Fixed.Ranger.Model( x, y, shore.formula )
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.HG <- foo$Model
build.results <- c( build.results, 'HG' = list( 'stats' = foo$Stats, 'Import' = imp) )


#----- Region NCC  -----
a <- obs.20mIV$NCC
x <- a[ a$TestDat == 0, ]
y <- a[ a$TestDat == 1, ]

foo <- Fixed.Ranger.Model( x, y, shore.formula )
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.NCC <- foo$Model
build.results <- c( build.results, 'NCC' = list( 'stats' = foo$Stats, 'Import' = imp) )


#----- Region WCVI  -----
a <- obs.20mIV$WCVI
x <- a[ a$TestDat == 0, ]
y <- a[ a$TestDat == 1, ]

foo <- Fixed.Ranger.Model( x, y, shore.formula )
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.WCVI <- foo$Model
build.results <- c( build.results, 'WCVI' = list( 'stats' = foo$Stats, 'Import' = imp) )


#----- Region QCS  -----
a <- obs.20mIV$QCS
x <- a[ a$TestDat == 0, ] 
y <- a[ a$TestDat == 1, ] 

foo <- Fixed.Ranger.Model( x, y, shore.formula )
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.QCS <- foo$Model
build.results <- c( build.results, 'QCS' = list( 'stats' = foo$Stats, 'Import' = imp) )


#----- Region SOG  -----
a <- obs.20mIV$SOG
x <- a[ a$TestDat == 0, ] 
y <- a[ a$TestDat == 1, ] 

foo <- Fixed.Ranger.Model( x, y, shore.formula )
imp <- foo$Model$variable.importance

# Store results  ... 
rf.region.SOG <- foo$Model
build.results <- c( build.results, 'SOG' = list( 'stats' = foo$Stats, 'Import' = imp) )

#-----------------------------
#-- Done building ranger RF models. 

endtime <- Sys.time()
cat( 'Model build time: ', endtime - starttime, '\n\n' )
#last run for 1 iteration was 7.4 minutes.

#--------------------------------------------
#-- SAVE the resulting models. Takes minutes - its a big file. 

# Build a time stamp ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( rf.region.Coast, rf.region.HG, rf.region.NCC, 
      rf.region.WCVI, rf.region.QCS, rf.region.SOG, 
      file = file.path( model.dir, paste0('rf_allModels_', x, '.RData')) )

save( build.results, file = file.path( model.dir, paste0('buildResults_', x, '.RData')) )

build.results$Coast.stats[2]
build.results$WCVI.stats[1]

#--------------------------------------------------------------------------
#--- Create summaries of model builds for figure production
build.sum <- Summarize.Build( build.results )
# BOTH as csv files and data frames.
#   Build_results_Integrated.csv - build.results.Integrated
#   Build_results_byClassStats.csv - build.results.ByClass
#   Build_results_testPrevalence.csv - build.results.ClassPrev
#   Build_results_varImportance.csv - build.results.VarImport
#   CSVs are good for tables but loose useful DF bits ... 

#------- Produce and save some plots ----------------

#--- Heatmaps (tigures!) of build stats as png files from above csv files.

Plot.Build.Class.Stats( 'Build_results_byClassStats.csv', pal.10, 800, 600 )

Plot.Build.Var.Import( 'Build_results_varImportance.csv', pal.11, 1000, 600 )
  

#--- Facets of the Build results ... 

#--- Prevalence plots (obs & predicted) for build (testing partition), faceted by Region
a <- Plot.Obs.Pred.Prevalence.Build( build.sum$build.results.ClassPrev, pal.cole ) 
ggsave("Class Prevalence for Obs&Pred (Build) for Regions.png", a, dpi = 300, width = 16, height = 10, path = output.dir)

#--- Prevalence plots (obs vs. predicted space) faceted by Region.
#-- Build a prevalence table for the regional predictions ... 
# DEFERRED cuz raster scaping. 
# a <- Plot.Obs.RegionPred.Prevalence( build.sum$build.results.ClassPrev, pal.cole ) 
# ggsave("Class Prevalence for Obs vs Study area for Regions.png", a, dpi = 300, width = 16, height = 10, path = output.dir)


#-- User and producer accuracies by class ... 
a <- Plot.Class.Stats.For.Regions( build.sum$build.results.ByClass, pal.cole )
ggsave( 'Class Stats (Build) for Regions.png', a, dpi = 300, width = 16, height = 10, path = output.dir)

#-- Main statistics, faceted ... doesn't make a very nice bar plot ... table is better. 


#====================================================================================
#-- How does the 100m model perform regionally?

#-- Feed the SpatialPointsDataFrame containing the testing partition
# To a function to split it into regions ... 
#-- NOTE: Obs data already partitioned, but has the 20 m data assigned ... 

test.regions <- Partition.Test.Data( point.data$Obs[ point.data$Obs$TestDat == 1, ] )
length( test.regions )
names( test.regions )
names( test.regions$NCC )

str( test.regions$NCC )

x <- rbind( 
  cbind( 'Region' = 'Coast', Results.Row( rf.region.Coast, test.regions$Coast )$Integrated), 
  cbind( 'Region' = 'HG', Results.Row( rf.region.Coast, test.regions$HG )$Integrated), 
  cbind( 'Region' = 'NCC', Results.Row( rf.region.Coast, test.regions$NCC )$Integrated),
  cbind( 'Region' = 'WVVI', Results.Row( rf.region.Coast, test.regions$WCVI )$Integrated), 
  cbind( 'Region' = 'QCS', Results.Row( rf.region.Coast, test.regions$QCS )$Integrated),
  cbind( 'Region' = 'SOG', Results.Row( rf.region.Coast, test.regions$SOG )$Integrated)
)

out.file <- '100m_Performance_Regionally.csv'
write.csv( x, file = file.path(results.dir, out.file) )

#-- Point of above was to see if 100m model performance was skewed in regions. Not really. 
#   Some question re: differences in regional sample sizes compared to 20 m. Defer unless needed.


#====================================================================================
#-- Very meaningful section title ... 

# Test 100 m model by depth class for each region ... 

x <- rbind(
  cbind( 'Region' = 'HG', Coast.Fit.By.Region.By.Depth( test.regions$HG )),
  cbind( 'Region' = 'NCC', Coast.Fit.By.Region.By.Depth( test.regions$NCC )),
  cbind( 'Region' = 'WCVI', Coast.Fit.By.Region.By.Depth( test.regions$WCVI )),
  cbind( 'Region' = 'QCS', Coast.Fit.By.Region.By.Depth( test.regions$QCS )),
  cbind( 'Region' = 'SOG', Coast.Fit.By.Region.By.Depth( test.regions$SOG )) )

#====================================================================================
#-- Within each region, how do the 20 m models perform across depth classes
#   Each model uses its own withheld Obs data.

y <- rbind(
#  cbind( 'Region' = 'Coast',Model.Fit.Obs.By.Depth( 'Coast', 3 )),
  cbind( 'Region' = 'HG',   Model.Fit.Obs.By.Depth( 'HG', 3 )),
  cbind( 'Region' = 'NCC',  Model.Fit.Obs.By.Depth( 'NCC', 3 )),
  cbind( 'Region' = 'WCVI', Model.Fit.Obs.By.Depth( 'WCVI', 3 )),
  cbind( 'Region' = 'QCS',  Model.Fit.Obs.By.Depth( 'QCS', 3 )),
  cbind( 'Region' = 'SOG',  Model.Fit.Obs.By.Depth( 'SOG', 3 )) )


#-- Single stat. No legend required. 
#   Looked at Accuracy and TNR and they contribute little.
Plot.TSS.By.Depth.For.Regions( x[ , c('Region', 'Ribbon', 'TSS')], pal.cole[3] )


#-- Above is fine, but to be fair really want the 100 m model evaluated for each region ... 

z <- rbind( 
  cbind( Model = '100 m', x ),
  cbind( Model = '20 m', y ) )

Plot.TSS.By.Depth.For.Regions2( z[ , c('Model', 'Region', 'Ribbon', 'TSS')], pal.cole[-3] )


#-- If you want to know what's driving TSS ... 
names(z)
cor( z[, -c(1:3)] )
a <- lm( TSS ~ Model + Region + Ribbon + N + Imbalance, data = z )
a <- lm( TSS ~ Model + Region + Ribbon, data = z )
summary(a)

#-- What about other metrics?
a <- lm( Accuracy ~ Model + Region + Ribbon, data = z )
anova( a )



#====================================================================================
#-- Create study area-wide predictions so that can compare prevalence of 
# study area prediction to training data ... 

#- coastwide 
coast.stack <- Load.Predictors( predictors.coastwide )
# standardize var names and bathymetry sign
names(coast.stack)[1] <- 'bathy'
coast.stack$bathy <- coast.stack$bathy * -1
#hist(coast.stack$bathy)

#- HG
HG.stack <- Load.Predictors( paste0( predictor.dir, 'HG' ) )
names(HG.stack)


names(rf.region.Coast)

#-- BOTH fail with 'not a multple' error ... 
foo <- raster::predict( coast.stack, rf.region.Coast )
foo <- raster::predict( HG.stack, rf.region.HG, 'text' )


print(rf.region.Coast$variable.importance)
print(rf.region.Coast$predictions)

treeInfo( rf.region.Coast, 1)


# coles...
#writeRaster(foo, file.path(single.step.dir, single.step), format = 'GTiff', datatype = 'INT2S')



#====================================================================================
#-- Independent Data Evaluation --
# Comparisons include:
#   1. Coastwide model vs. all 3 ID sets.
#     ID merged from the regions to test coastwide.
#     summarizes the merged IDS
#   2. Regional models vs. all 3 ID sets. 
#     summarizes the regional IDS


#--- 1. Test each model (n=6) vs. each of the IDS.

#-- Load it ... 
results.table <- read.csv( file.path(output.dir, 'IDS_evaluation_by_region.csv'))

#-- or build it ...
results.table <- Do.Independent.Evaluation()

out.file <- 'IDS_evaluation_by_region.csv'
write.csv( results.table, file = file.path(output.dir, out.file) )

#--------------------------------------
#-- Present the initial IDE results ... 

x <- build.sum$build.results.Integrated
x[1,1]

#--- Integrated statistics for each IDS by Region ---
# Use a predetermined collection of stats (see function)
a <- Plot.Stat.By.IDS.For.Regions( results.table, pal.3 )
ggsave( paste0( 'Integrated Stats by IDE for Regions.png'), a, dpi = 300, width = 16, height = 10, path = output.dir)


#-- UNDER CONSTRUCTION - need some by-class results ... 

#-- User and producer accuracies by class ---
#####   NEEDS A NEW DATA STRUCTURE, SO requires different/updated IDs SUMMARY CODE ####
a <- IDS.Class.Stats.For.Regions( results.table, pal )
ggsave( 'Class Stats (IDE) for Regions.png', a, dpi = 300, width = 16, height = 10, path = output.dir)






#-- FIN.
