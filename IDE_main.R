#*******************************************************************************
# Script:  IDE_main.R
# Created: January 2020. EJG
# Updated: 2020/10/06 EJG (updates now intermittently detailed in README.md)
#
# This script sources the necessary libraries and functions, then uses flags to 
# either loads results from saved Rdata files, or import and analyse raw data. 
# The data build/loaded here gets picked up by th R markdown file.

#*******************************************************************************

rm(list=ls(all=T))  # Erase environment.

#-- Load necessary packages and functions ... 
#source( "depth_effect.R" )
source( "substrate_functions.R" )
source( "Plot_Functions.R" )


#--------------------------------------------------------------------------
#-- Part 1. Load Source Points and Results 

# Processsing FLAGS. When set to true, the data structure will be re-built from imported data. 
# Otherwise, it will be loaded from rData. 
# THere are dependencies, e.g., if the RF models are re-built, the analysis  needs to be redone.

# Parameters to the re-build ... 
reloadpts  <- F
wtdrngr    <- F
nowtrngr   <- F
makeheat   <- F
mapplots   <- F

#---- Step 1 of 2: Load the things that don't need to be built ... 

#-- RData file names:
data.files <- list(
  pts        = 'loaded_data_2020-11-19.RData',
  modelsWtd  = 'rf_allModels_2020-11-19.RData',
  buildRes   = 'buildResults_2020-11-19.RData',
  modelsNoWt = 'nwrf_allModels_2020-11-19.RData',
  buildResNW = 'nwbuildResults_2020-11-19.RData',
  predMaps   = 'rasterMapObjects_2020-11-20.RData'
)

# 1) Source data
if (reloadpts == F){
#-- Load existing observational data .. 
  load( file.path( model.dir, data.files[["pts"]] ))
# point.data - All the point observations from ArcGIS. 
# obs.100mIV
# obs.20mIV, dive.20mIV, cam.20mIV, ROV.20mIV
# names.100m
}

# 2): RF Models and Build results. 
if (wtdrngr == F){
  #-- Load the latest RF models (This is a BIG file) ... 
  load( file.path( model.dir, data.files[["modelsWtd"]] ))
  # rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI, rf.region.QCS, rf.region.SOG

  #  Also want the associated build statistics and results ... 
  load( file.path( model.dir, data.files[["buildRes"]] ))
  # build.results
}

if (nowtrngr == F){
  #-- Load the latest RF models (This is a BIG file) ... 
  load( file.path( model.dir, data.files[["modelsNoWt"]] ))
  # rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI, rf.region.QCS, rf.region.SOG
  
  #  Also want the associated build statistics and results ... 
  load( file.path( model.dir, data.files[["buildResNW"]] ))
  # build.results
}

#-------------------------------------------
#  Load the predicted rasters and prevalences (map.prev) ... 
if (mapplots == F){
  load( file.path( model.dir, data.files[["predMaps"]] ))
}

# end of data/results load

#---------------------------------------------------------------------------------
#-- Step 2 of 2: Source the build script to build any flagged results/objects.
source( "build_substrate.R" )


#--------------------------------------------------------------------------
#-- Part 2. Summarize Results 

build.sum <- Summarize.Build( build.results )

#----
# Duplicate the above to summarize the no-weight build. 
build.sum.nw <- Summarize.Build( build.results.nw )


#--------------------------------------------------------------------------
#-- Part 3. Model resolution tests

#-- A) How does the 100m model perform at the regional scale?
 
# First separate the Obs testing data by region. 
test.regions <- Partition.Test.Data( point.data$Obs[ point.data$Obs$TestDat == 1, ] )

#-- B) How do the two model resolutions compare across depths?
#   Assess ribbons within each region

# 1) Test the 100 m model by depth  within each region. 
# 2) Test the regional models; each using its own 20 m Obs testing data.
# test.regions contain the 100 m test data partitioned above. 

depth.results <- Models.Across.Depths( test.regions )


#-------------------------------------------------------------------
#-- Part 4 - Build the Heat maps 
#  Needs to be done here cuz behaviour is odd in the RMD file ... 

if (makeheat == T){
  
  Heat.Build.Class.Stats( build.sum$build.results.ByClass, T, 'TPR', rev( pal.heat.10 ), 800, 600, 'black' )
  Heat.Build.Class.Stats( build.sum$build.results.ByClass, T, 'TNR', rev( pal.heat.10 ), 800, 600, 'black' )
  Heat.Build.Class.Stats( build.sum$build.results.ByClass, T, 'User', rev( pal.heat.10 ), 800, 600, 'black' )
  
  Heat.Build.Class.Stats( build.sum.nw$build.results.ByClass, F, 'TPR', rev( pal.heat.10 ), 800, 600, 'black' )
  Heat.Build.Class.Stats( build.sum.nw$build.results.ByClass, F, 'TNR', rev( pal.heat.10 ), 800, 600, 'black' )
  Heat.Build.Class.Stats( build.sum.nw$build.results.ByClass, F, 'User', rev( pal.heat.10 ), 800, 600, 'black' )
  
  Heat.Build.Var.Import( build.sum$build.results.VarImport, pal.heat.11, 1000, 600, 'black' )
}

#--------------------------------------------------------------
#-- Regression to see what's driving the different metrics  ... 
#   Have tested with/without N and Imbalance. No change. 
#   REPLACE TSS with other metrics to create the results 

# names(depth.results)
# cor( depth.results[, -c(1:3)] )
# a <- lm( TSS ~ Model + Region + Ribbon + N + Imbalance, data = depth.results )
# summary(a)
# anova( a )


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

IDE.results.wtd  <- IDS.Evaluation( 'rf' )
row.names( IDE.results.wtd ) <- NULL

# IDE.results.trm  <- IDS.Evaluation( 'trm' )
IDE.results.nowt <- IDS.Evaluation( 'nwrf' )
row.names( IDE.results.nowt ) <- NULL


#-------------------------------------------------------------
#---- Make some plots with the above ... but so many plots ... 
#---- USES the data stucture currently built for Table S2 in rmd file.

# y <- x[ x$Stat == 'TPR', c('Model', 'Region', 'Hard', 'Mixed', 'Sand', 'Mud')]
# a <- Plot.ClassStats.IDE( y, 'Accuracy', pal.cb3b, sz=30, lx=0, ly=0 )
#   
# y <- x[ x$Stat == 'TNR', c('Model', 'Region', 'Hard', 'Mixed', 'Sand', 'Mud')]
# b <- Plot.ClassStats.IDE( y, 'Specificity', pal.cb3b, sz=30, lx=0, ly=0 )
# 
# y <- x[ x$Stat == 'User', c('Model', 'Region', 'Hard', 'Mixed', 'Sand', 'Mud')]
# c <- Plot.ClassStats.IDE( y, 'Reliability', pal.cb3b, sz=30, lx=0, ly=0 )


#-------------------------------------------------------------
#-- IDS test results across depth ribbons by region. 
# Used to plot of Pontius stats 
# 2020/07/22: Modified the Fit function to use the extended depth zones. 
IDE.depths <- rbind(
  cbind( 'IDS' = 'Dive', rbind( 
    cbind( 'Region' = 'Coast', Model.Fit.IDS.By.Depth( 'dive', 'Coast' )),
    cbind( 'Region' = 'HG',    Model.Fit.IDS.By.Depth( 'dive', 'HG' )),
    cbind( 'Region' = 'NCC',   Model.Fit.IDS.By.Depth( 'dive', 'NCC' )) )
  ),
  cbind( 'IDS' = 'Cam', rbind( 
    cbind( 'Region' = 'Coast', Model.Fit.IDS.By.Depth( 'cam', 'Coast' )),
    cbind( 'Region' = 'HG',    Model.Fit.IDS.By.Depth( 'cam', 'HG' )),
    cbind( 'Region' = 'NCC',   Model.Fit.IDS.By.Depth( 'cam', 'NCC' )) )
  ),
  cbind( 'IDS' = 'ROV', rbind( 
    cbind( 'Region' = 'Coast', Model.Fit.IDS.By.Depth( 'ROV', 'Coast' )),
    cbind( 'Region' = 'HG',    Model.Fit.IDS.By.Depth( 'ROV', 'HG' )),
    cbind( 'Region' = 'NCC',   Model.Fit.IDS.By.Depth( 'ROV', 'NCC' )) )
  )
)

row.names( IDE.depths ) <- NULL


#-- FIN.
