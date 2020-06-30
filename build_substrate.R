#======= READ source and BUILD data structures ===========================

#-------------------------------------------------------------------------
#-- PART 1: Load, prep, and inspect all the data. 
#-- Wrapper functions used to call other functions and keep it clean here. 
#-- RE-BUILDING observational data can take 20-30 min.
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
obs.100mIV <- obs.100mIV[ , !names(obs.100mIV) %in% drop.list ]

obs.100mIV <- Rename.100m.Preds( obs.100mIV )


#--------- Done loading source data --------------
endtime <- Sys.time()
cat('---------------------------\n')
cat( 'Data build time: ', endtime - starttime, '\n\n' )

#-----------------------------
#-- Save out the populated observational data.
# Build a date/time stamp  ... 
# x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

# Build a date stamp  ... 
x <- substr( endtime, 1, 10)


save( point.data, dive.20mIV, cam.20mIV, ROV.20mIV, 
      obs.100mIV, obs.20mIV,
      file = file.path(model.dir, paste0('loaded_data_', x, '.RData') ))


#=====================================================================================
#-- PART 2: Build and evaluate 6 RF models. (one Coastwide, 100m; 5 regional, 20m)
# All models use wtd ranger(). 
# Previoulsly built using Rand.Ranger.Model but simplified now to surface train/test data.
# Fixed.Ranger.Model() uses only the train data to parameterize. 

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

#----------------------------------
#-- Done building ranger RF models. 

endtime <- Sys.time()
cat( 'Model build time: ', endtime - starttime, '\n\n' )
#last run for 1 iteration was 2.7 minutes.

#--------------------------------------------
#-- SAVE the resulting models. Takes minutes - its a big file. 

# Build a time stamp ... 
#x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

# Build a date stamp ... 
x <- substr( endtime, 1, 10)

save( rf.region.Coast, rf.region.HG, rf.region.NCC, 
      rf.region.WCVI, rf.region.QCS, rf.region.SOG, 
      file = file.path( model.dir, paste0('rf_allModels_', x, '.RData')) )

save( build.results, file = file.path( model.dir, paste0('buildResults_', x, '.RData')) )



#-------------------------------------------
#  Build model predictions. TAKES HOURS!

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
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.RMSM )
map.prev <- c(map.prev, 'HG' = y )

#- NCC
a <- 'NCC'
b <- rf.region.NCC
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.RMSM )
map.prev <- c(map.prev, 'NCC' = y )

#- WCVI
a <- 'WCVI'
b <- rf.region.WCVI
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.RMSM )
map.prev <- c(map.prev, 'WCVI' = y )

#- QCS
a <- 'QCS'
b <- rf.region.QCS
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.RMSM )
map.prev <- c(map.prev, 'QCS' = y )

#- SOG
a <- 'SOG'
b <- rf.region.SOG
a.stack <- Load.Predictors( paste0( predictor.dir, '/', a ) )
y <- Predict.Surface( a.stack, b, raster.dir, a, pal.RMSM )
map.prev <- c(map.prev, 'SOG' = y )


#-- SAVE the prevalence and the predicted objects ... 

#Build a date stamp ... 
endtime <- Sys.time()
x <- substr( endtime, 1, 10)

save( map.prev,
      file = file.path( model.dir, paste0('rasterMapObjects_', x, '.RData')) )


### Fin.