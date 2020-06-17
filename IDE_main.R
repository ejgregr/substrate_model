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
#   2020/06/15: re-structure with all plotting happening in RMD file. 
#     This script now makes data structures and write CSV summaries, picked up by RMD.
#     The process is documented in the RMD.
#
#*******************************************************************************

rm(list=ls(all=T))  # Erase environment.

#-- Load necessary packages and functions ... 
source( "substrate_functions.R" )
source( "Plot_Functions.R" )

#-------------------------------------------------------------------------
#-- PART 1: Load, prep, and inspect all the data. 
#-- Wrapper functions used to call other functions and keep it clean here. 

#-- Load existing observational data .. 
load( file.path( model.dir, 'loaded_data_2020-06-01-1257.RData' ))
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
      obs.100mIV, obs.20mIV, names.100m,
      file = file.path(model.dir, paste0('loaded_data_', x, '.RData') ))


#=================================================================================
#-- PART 2: Build and evaluate the required RF models. 
# All models use wtd ranger(). 
# Models include:
#   1. Coastwide, 100 m
#   2. Regional, 20 m

#-- Load the latest RF models (This is a BIG file) ... 
load( file.path( model.dir, 'rf_allModels_2020-05-21-1019.RData' ))

#  Also want the build statistics and results ... 
load( file.path( model.dir, 'buildResults_2020-05-21-1019.RData' ))

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
#last run for 1 iteration was 2.7 minutes.

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
#   CSVs are good for tables but lose useful DF bits ... 


#---Part 3: Produce and save some build results plots ------------

#--- Heatmaps (tigures!) of build stats as png files from above csv files.

Plot.Build.Class.Stats( build.sum$build.results.ByClass, pal.10, 800, 600 )

Plot.Build.Var.Import( build.sum$build.results.VarImport , pal.11, 1000, 600 )
  

#-- Facets of the Build results ... 

#-- Prevalence plots (obs & predicted) for build (testing partition), faceted by Region
a <- Plot.Obs.Pred.Prevalence.Build( build.sum$build.results.ClassPrev, pal.cole ) 
ggsave("Class Prevalence for Obs&Pred (Build) for Regions.png", a, dpi = 300, width = 16, height = 10, path = output.dir)

#-- Prevalence plots (obs vs. predicted space) faceted by Region.
#   PENDING - waiting for raster scaping. 
# a <- Plot.Obs.RegionPred.Prevalence( build.sum$build.results.ClassPrev, pal.cole ) 
# ggsave("Class Prevalence for Obs vs Study area for Regions.png", a, dpi = 300, width = 16, height = 10, path = output.dir)

#-- User and producer accuracies by class ... 
a <- Plot.Class.Stats.For.Regions( build.sum$build.results.ByClass, pal.cole )
ggsave( 'Class Stats (Build) for Regions.png', a, dpi = 300, width = 16, height = 10, path = output.dir)

#-- Main statistics (Build_results_Integrated.csv) don't make a very nice bar plot. Report as table. 


#======================================================================
#---Part 4: Model resolution tests 

#-- 4a: How does the 100m model perform regionally?
#   Requires patitioning of the Obs testing data with 100 m predictors. 
#   Use the spatial info on the points to separate using the Region shape file.
#   NOTE: This is not the same as the Obs partitioned and loaded w the 20 m predictors.

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

#-- Above shows that 100m model performance is NOT very skewed by regions.


#------------------------------------------------------------------
#-- Part 4b: How do the two model resolutions compare across depths?
#     Assess ribbons within each region

# First test the 100 m model by ribbon within each region. 
#   test.regions are the 100 m test data partitioned above. 

x <- rbind(
  cbind( 'Region' = 'HG', Coast.Fit.By.Region.By.Depth( test.regions$HG )),
  cbind( 'Region' = 'NCC', Coast.Fit.By.Region.By.Depth( test.regions$NCC )),
  cbind( 'Region' = 'WCVI', Coast.Fit.By.Region.By.Depth( test.regions$WCVI )),
  cbind( 'Region' = 'QCS', Coast.Fit.By.Region.By.Depth( test.regions$QCS )),
  cbind( 'Region' = 'SOG', Coast.Fit.By.Region.By.Depth( test.regions$SOG )) )

# Now test the regional models; each model uses its own, named 20 m Obs testing data.

y <- rbind(
#  cbind( 'Region' = 'Coast',Model.Fit.Obs.By.Depth( 'Coast', 3 )),
  cbind( 'Region' = 'HG',   Model.Fit.Obs.By.Depth( 'HG', 3 )),
  cbind( 'Region' = 'NCC',  Model.Fit.Obs.By.Depth( 'NCC', 3 )),
  cbind( 'Region' = 'WCVI', Model.Fit.Obs.By.Depth( 'WCVI', 3 )),
  cbind( 'Region' = 'QCS',  Model.Fit.Obs.By.Depth( 'QCS', 3 )),
  cbind( 'Region' = 'SOG',  Model.Fit.Obs.By.Depth( 'SOG', 3 )) )

# Combine the two results ... 

depth.results <- rbind( 
  cbind( Model = '100 m', x ),
  cbind( Model = '20 m', y ) )
row.names( depth.results ) <- NULL

str( depth.results )

out.file <- 'Stats_byRibbon_byRegion_byModel.csv'
write.csv( depth.results, file = file.path(results.dir, out.file) )

Plot.TSS.By.Depth.For.Regions( depth.results[ , c('Model', 'Region', 'Ribbon', 'TSS')], pal.cole[-3] )


#-- vs. What does a tigure look like? Not so good ... 
out.file <- 'TSS_byRibbon_byRegion_byModel.csv'
write.csv( depth.results[ , c("Model", "Region", "Ribbon", "TSS")], file = file.path(results.dir, out.file) )

Plot.Region.Class.Stats( 'TSS_byRibbon_byRegion_byModel.csv', '20 m', pal.10 )
Plot.Region.Class.Stats( 'TSS_byRibbon_byRegion_byModel.csv', '100 m', pal.10 )




#-------------------------------------------------------------------------
#-------- Do all the above a second time, with the longer list of ribbons.

x <- rbind(
  cbind( 'Region' = 'HG', Coast.Fit.By.Region.By.Depth2( test.regions$HG )),
  cbind( 'Region' = 'NCC', Coast.Fit.By.Region.By.Depth2( test.regions$NCC )),
  cbind( 'Region' = 'WCVI', Coast.Fit.By.Region.By.Depth2( test.regions$WCVI )),
  cbind( 'Region' = 'QCS', Coast.Fit.By.Region.By.Depth2( test.regions$QCS )),
  cbind( 'Region' = 'SOG', Coast.Fit.By.Region.By.Depth2( test.regions$SOG )) )

# Now test the regional models; each model uses its own, named 20 m Obs testing data.

y <- rbind(
  #  cbind( 'Region' = 'Coast',Model.Fit.Obs.By.Depth( 'Coast', 3 )),
  cbind( 'Region' = 'HG',   Model.Fit.Obs.By.Depth2( 'HG', 3 )),
  cbind( 'Region' = 'NCC',  Model.Fit.Obs.By.Depth2( 'NCC', 3 )),
  cbind( 'Region' = 'WCVI', Model.Fit.Obs.By.Depth2( 'WCVI', 3 )),
  cbind( 'Region' = 'QCS',  Model.Fit.Obs.By.Depth2( 'QCS', 3 )),
  cbind( 'Region' = 'SOG',  Model.Fit.Obs.By.Depth2( 'SOG', 3 )) )

# Combine the two results ... 

depth.results2 <- rbind( 
  cbind( Model = '100 m', x ),
  cbind( Model = '20 m', y ) )
row.names( depth.results2 ) <- NULL


out.file <- 'Stats_byRibbon2_byRegion_byModel.csv'
write.csv( depth.results, file = file.path(results.dir, out.file) )

Plot.TSS.By.Depth.For.Regions( depth.results2[ , c('Model', 'Region', 'Ribbon', 'TSS')], pal.cole[-3] )


#-- vs. What does a tigure look like? Not so good ... 
out.file <- 'TSS_byRibbon2_byRegion_byModel.csv'
write.csv( depth.results[ , c("Model", "Region", "Ribbon", "TSS")], file = file.path(results.dir, out.file) )

Plot.Region.Class.Stats( 'TSS_byRibbon2_byRegion_byModel.csv', '20 m', pal.10 )
Plot.Region.Class.Stats( 'TSS_byRibbon2_byRegion_byModel.csv', '100 m', pal.10 )








#--------------------------------------------------------------------
#-- Some regression to see what's driving the different metrics  ... 
#   Have tested with/without N and Imbalance. No change. 
#   REPLACE TSS with other metrics to create the results 
names(depth.results)
cor( depth.results[, -c(1:3)] )
a <- lm( TSS ~ Model + Region + Ribbon + N + Imbalance, data = depth.results )
summary(a)
anova( a )



#=======================================
#-- Part 5 - Independent Data Evaluation
# Comparisons include:
#   1. Coastwide model vs. all 3 ID sets.
#     ID merged from the regions to test coastwide.
#     summarizes the merged IDS
#   2. Regional models vs. all 3 ID sets. 
#     summarizes the regional IDS


#--- 1. Test each model (n=6) vs. each of the IDS.
# Build new every time to ensure its fresh. Doesn't take long.
# Question(2020/05/18): Since WCVI, QCS, and SOG regions now dropped because of low sample size,
#   is it necessary to drop them  out of the Coast evaluation here? Prob not. Just don't get used.

IDE.results <- Do.Independent.Evaluation()

#--------------------------------------
#-- Present some Initial IDE results ... 


#--- 1) Begin by looking at IDS sample sizes using the perClass prevalences 

### Code section copied to RMD file (Figure 7). 

#----- Above results suggest IDS focus needs to be on Coast, HG, and NCC ... 
#     This should apply to the data density test as well then. See below.


#--- 2) Look at class statistics by IDS ... User/producer accuracies by regions by IDS. 

### Challenge here is there is an extra dimension: Stat, Class, Region, IDS ...
#     Could stack a bar within a facet but ... :\
#    Plot for each region? For each IDS? Not clear. :(

#--- Assemble the data ... 
x <- IDE.results$PerClass
y <- x[ x$Stat %in% c( 'User', 'Prod'), ]
y <- y[ y$Region %in% c( 'Coast', 'HG', 'NCC'), ]

# Tidy the data a bit  ... 
colnames(y) <- c( 'Region', 'IDS', 'Stat', 'Rock', 'Mixed', 'Sand', 'Mud')
rownames(y) <- NULL

a <- IDS.Class.Stats.For.Regions( y, pal.3.win, sz = 25, lx=0.83, ly=0.75 )


#---------------------
# 3) Look at Integrated statistics for each IDS by Region.
# Use a predetermined collection of stats (see function)

# Again, start with a little data prep ... 

y <- IDE.results$Integrated
z <- y[ y$Region %in% c('Coast', 'HG', 'NCC'), ]
rownames(z) <- NULL
z
  
# This function can be used to generate various views of the data. 
# 2020/05/25 DH preferred IDS as the grouping variable ... 
Plot.Stats.By.IDS.For.Regions( z, pal.3.win, sz = 25, lx=0.76, ly=0.87 )
ggsave( paste0( 'Integrated Stats by IDE for Regions.png'), a, dpi = 300, width = 16, height = 10, path = output.dir)


# What about TSS vs. Quantity and Allocation?
# Remember: Accuracy = 1 – (Quantity + Exchange + Shift).
zz <- z$Shift + z$Exchange
zz <- zz[, c( 'Region', 'IDS', 'Accuracy', 'Quantity', 'Exchange','Shift' )]
head(zz)

Plot.Pontius.By.IDS.For.Regions( zz, rev(pal.4), sz = 25 )


#=========================================== 
#--- Part 6 - Test data density effect -----

#-- Use function to return lists of Obs testing data in low/high density regions.
#   Calls shape file directly. 

dens.lists <- Partition.By.Density( point.data$Obs[ point.data$Obs$TestDat == 1, ])

# Some double checking: lengths and names ... 
dim( point.data$Obs[ point.data$Obs$TestDat == 1, ] )
dim( dens.lists$Dense  )
dim( dens.lists$Sparse )
names( dens.lists$Sparse )

# So ... do the two data sets evaluate differently?

x <- rbind(
  cbind( 'Density' = 'Low',  Results.Row( rf.region.Coast, dens.lists$Sparse )$Integrated ), 
  cbind( 'Density' = 'High', Results.Row( rf.region.Coast, dens.lists$Dense )$Integrated )
)


y <- Partition.By.Density( point.data$Obs[ point.data$Obs$TestDat == 1, ])





#-- Build the table piece ... 
compare.what <- data.frame( 'Region' = 'Coast', 'IDS' = 'Dive' )
w <- Results.Row( rf.region.Coast, x.test )
x <- cbind( compare.what, w$Integrated )
results.int <- rbind( results.int, x ) 



#--- Potential plots... 
#     Compare the IDs results to some of the Obs testing results if desired ...Options = 
#     1) compare model ranks across test data; 
#     2) test how sample size affects results? What about empty classes?



#=============================================
#-- Part 7 - Spatial substrate predictions ... 
#-- Create study area-wide predictions so that can compare prevalence of 
# study area prediction to training data ... 


#-------------------------------------------
#  Load the predicted rasters and prevalences  ... 
load( file.path( model.dir, 'rasterMapObjects_2020-05-21-1019.RData' ))

#-------------------------------------------
#  Or build some new ones here. TAKES HOURS!

#- Somewhere to put the data ... 
map.prev <- list()

#- coastwide 
a <- 'Coast'
b <- rf.region.Coast
a.stack <- Load.Predictors( paste0( predictor.dir, '/Coastwide' ) )

# standardize var names and bathymetry sign
names(a.stack)[1] <- 'bathy'
a.stack$bathy <- a.stack$bathy * -1

y <- Predict.Surface( a.stack, b, raster.dir, a, pal.map )
map.prev <- c(map.prev, 'Coast' = y )

#- HG
a <- 'HG'
b <- rf.region.HG
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.map )
map.prev <- c(map.prev, 'HG' = y )

#- NCC
a <- 'NCC'
b <- rf.region.NCC
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.4 )
map.prev <- c(map.prev, 'NCC' = y )

#- WCVI
a <- 'WCVI'
b <- rf.region.NCC
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.4 )
map.prev <- c(map.prev, 'WCVI' = y )

#- QCS
a <- 'QCS'
b <- rf.region.NCC
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.4 )
map.prev <- c(map.prev, 'QCS' = y )

#- SOG
a <- 'SOG'
b <- rf.region.SOG
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.4 )
map.prev <- c(map.prev, 'SOG' = y )


#-- SAVE the prevalence and the predicted objects ... 

#Build a time stamp ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( map.prev,
      file = file.path( model.dir, paste0('rasterMapObjects_', x, '.RData')) )


#----- Plot the Study Area prevalence results ------

# Prepare the prevalence tables for each region, knock into shape, and ggplot() ... 
# Called by RMD: 
#   Plot.Pred.Map.Prevalence( map.prev, build.sum)


#---- Manual plotting of predicted tifs ----
#  UNDER DEVELOPMENT - Saving and plotting predictions has become a bit of a pain. 

x <- pal.map

a <- 'QCS'
#y <- unlist( map.prev$HG1 )

b <- paste0( raster.dir, '/',a, '_classified_substrate.tif')
y <- raster::raster( b )
str(y)

png(file=b,
    height = 7, width = 6, units = "in", res = 400)
raster::plot(y, maxpixels=5000000, col = x, legend=FALSE,
             xlab = "Easting", ylab = "Northing", cex.axis = .5, cex.lab = .75)
legend(x = "bottomleft",
       legend =  c("Rock", "Mixed", "Sand", "Mud"), fill = x, title=NA, bg = NA, box.col = NA)

dev.off()


# How to use the string in the list naming? ugh. i.e., 
map.prev <- c(map.prev, eval( parse(paste0( a, '=9876' )) ))

x[[2]]
map.prev






#-- FIN.
