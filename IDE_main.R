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
source( "depth_effect.R" )
source( "Plot_Functions.R" )

#-- This will re-load all the data and re-build the RF models ... 
source( "build_substrate.R" )


#======= Load Data and Results =================

#----------------------
#-- Part 1: Source data

#-- Load existing observational data .. 
load( file.path( model.dir, 'loaded_data_2020-06-25.RData' ))
# point.data - All the point observations from ArcGIS. 
# obs.100mIV
# obs.20mIV, dive.20mIV, cam.20mIV, ROV.20mIV
# names.100m

#---------------------------------------
#-- Part 2: RF Models and Build results. 

#-- Load the latest RF models (This is a BIG file) ... 
load( file.path( model.dir, 'rf_allModels_2020-06-25.RData' ))
# rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI, rf.region.QCS, rf.region.SOG

#  Also want the build statistics and results ... 
load( file.path( model.dir, 'buildResults_2020-06-25.RData' ))
# build.results

#build.results$Coast.stats[2]
#build.results$WCVI.stats[1]


#--------------------------------------------------------------------------
#-- Part 3: Create data summaries and CSVs for Heat maps.

build.sum <- Summarize.Build( build.results )
# BOTH as csv files and data frames.
#   Build_results_Integrated.csv - build.results.Integrated
#   Build_results_byClassStats.csv - build.results.ByClass
#   Build_results_testPrevalence.csv - build.results.ClassPrev
#   Build_results_varImportance.csv - build.results.VarImport
#   CSVs are good for tables but lose useful DF bits ... 


#--------------------------------------------------------------------------
#---Part 4: Model resolution tests 

#-- 4a: Test how the 100m model performs at the regions.

test.regions <- Regional.Test.100m( point.data$Obs[ point.data$Obs$TestDat == 1, ] )
#-- Write CSV of results showing that 100m model performance is NOT skewed by regions.


#------------------------------------------------------------------
#-- 4b: Test how the two model resolutions compare across depths.
#     Assess ribbons within each region

# First test the 100 m model by ribbon within each region. 
# Then test the regional models; each model uses its own 20 m Obs testing data.
# test.regions are the 100 m test data partitioned above. 
#-- Also writes CSV of results.

depth.results <- Models.Across.Ribbons( test.regions )

# Do it gain with the longer list of ribbons.
depth.results2 <- Models.Across.Ribbons2( test.regions )


#--------------------------------------------------------------------
#-- Some regression to see what's driving the different metrics  ... 
#   Have tested with/without N and Imbalance. No change. 
#   REPLACE TSS with other metrics to create the results 
names(depth.results)
cor( depth.results[, -c(1:3)] )
a <- lm( TSS ~ Model + Region + Ribbon + N + Imbalance, data = depth.results )
summary(a)
anova( a )


#--------------------------------------------------------------------------
#-- Part 5 - Independent Data Evaluation
# Comparisons include:
#   a. Coastwide model vs. all 3 ID sets.
#     ID merged from the regions to test coastwide.
#     summarizes the merged IDS
#   b. Regional models vs. all 3 ID sets. 
#     summarizes the regional IDS

#--- 5a. Test each model (n=6) vs. each of the IDS.
# Build new every time to ensure its fresh. Doesn't take long.
# Question(2020/05/18): Since WCVI, QCS, and SOG regions now dropped because of low sample size,
#   is it necessary to drop them out of the Coast evaluation here? Prob not. Just don't get used.

IDE.results <- Do.Independent.Evaluation()

#--------------------------------------
#-- Present some Initial IDE results ... 

#- A) Fig. 6: IDS sample sizes using the perClass prevalences 
# Results suggest IDS focus needs to be on Coast, HG, and NCC ... 
# This should apply to the data density test as well then. See below.

#---------------------
# B) Fig. 7: Look at Integrated statistics for each IDS by Region.
# Use a predetermined collection of stats (see function)



#-- 2) Look at class statistics by IDS ... User/producer accuracies by regions by IDS. 

#-- Assemble the data ... 
x <- IDE.results$PerClass
y <- x[ x$Region %in% c( 'Coast', 'HG', 'NCC'), ]
y <- y[ y$Stat %in% c( 'User', 'Prod'), ]

#-- FIX Stat levels to match earlier figure.  
y <- mutate(y, Stat = factor(Stat, levels=c('User', 'Prod', 'PrevObs', 'PrevPred' )))
            

# Tidy the data a bit  ... 
colnames(y) <- c( 'Region', 'IDS', 'Stat', 'Rock', 'Mixed', 'Sand', 'Mud')
rownames(y) <- NULL

a <- IDS.Class.Stats.For.Regions( y, pal.cb2, sz = 25, lx=0.83, ly=0.75 )
a



#--- Pontius stats for Obs data: ribbon by region. 


#--- Pontius stats for IDS data: ribbon by region. 






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
