#*******************************************************************************
# Script:  IDE_main.R
# Created: January 2020. EJG
# Updated: 2020/09/04 EJG (updates now intermittently detailed in README.md)
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
#source( "depth_effect.R" )
source( "Plot_Functions.R" )

#-- This will re-load all the data and re-build the RF models ... 
#  **** FIX this so it doesn't necessarily re-do all the plots. ****

# Parameters for the re-build ... 
reloadpts  <- F
wtdrngr    <- F
nowtrngr   <- F
trimrngr   <- F
mapplots   <- F
boptests   <- F


#======= Load Data and Results =================

#----------------------
#-- Part 1: Source data

if (reloadpts == F){
#-- Load existing observational data .. 
  load( file.path( model.dir, 'loaded_data_2020-08-10.RData' ))
# point.data - All the point observations from ArcGIS. 
# obs.100mIV
# obs.20mIV, dive.20mIV, cam.20mIV, ROV.20mIV
# names.100m
}
#---------------------------------------
#-- Part 2: RF Models and Build results. 

if (wtdrngr == F){
  #-- Load the latest RF models (This is a BIG file) ... 
  load( file.path( model.dir, 'rf_allModels_2020-08-24.RData' ))
  # rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI, rf.region.QCS, rf.region.SOG

  #  Also want the associated build statistics and results ... 
  load( file.path( model.dir, 'buildResults_2020-08-24.RData' ))
  # build.results
}

if (nowtrngr == F){
  #-- Load the latest RF models (This is a BIG file) ... 
  load( file.path( model.dir, 'nwrf_allModels_2020-09-02.RData' ))
  # rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI, rf.region.QCS, rf.region.SOG
  
  #  Also want the associated build statistics and results ... 
  load( file.path( model.dir, 'nwbuildResults_2020-09-02.RData' ))
  # build.results
}

if (trimrngr == F){
  #-- Load the latest RF models (This is a BIG file) ... 
  load( file.path( model.dir, 'trm_allModels_2020-08-31.RData' ))
  # rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI, rf.region.QCS, rf.region.SOG
  
  #  Also want the associated build statistics and results ... 
  load( file.path( model.dir, 'trmbuildResults_2020-08-31.RData' ))
  # build.results
}

#build.results$Coast.stats[2]
#build.results$WCVI.stats[1]


#-------------------------------------------
#  Load the predicted rasters and prevalences (map.pred) ... 
if (mapplots == F){
  load( file.path( model.dir, 'rasterMapObjects_2020-06-29.RData' ))
}


#-- Build any data/results not loaded above (uses the same flags) ... 
source( "build_substrate.R" )


#------------------------------------------------------------
#-- Part 3: Create data summaries for tables and figures.
build.sum <- Summarize.Build( build.results )
# BOTH as csv files and data frames.
#   Build_results_Integrated.csv - build.results.Integrated
#   Build_results_byClassStats.csv - build.results.ByClass
#   Build_results_testPrevalence.csv - build.results.ClassPrev
#   Build_results_varImportance.csv - build.results.VarImport
#   CSVs are good for tables but lose useful DF bits ... 

# REPRODUCE the above to create the no-weight build summaries. 
build.sum.nw <- Summarize.Build( build.results.nw )

# Just take the integrated stats for the trimmed build summaries for now. 
trm.build.sum <- rbind( 
  cbind( 'Region' = 'Coast', data.frame( trm.build.results$Coast.stats[ 1 ] ) ),
  cbind( 'Region' = 'HG',    data.frame( trm.build.results$HG.stats[ 1 ] ) ),
  cbind( 'Region' = 'NCC',   data.frame( trm.build.results$NCC.stats[ 1 ] ) ),
  cbind( 'Region' = 'WCVI',  data.frame( trm.build.results$WCVI.stats[ 1 ] ) ),
  cbind( 'Region' = 'QCS',   data.frame( trm.build.results$QCS.stats[ 1 ] ) ),
  cbind( 'Region' = 'SOG',   data.frame( trm.build.results$SOG.stats[ 1 ] ) )
)
    
  
#--------------------------------------------------------------------------
#---Part 4: Model resolution tests 

#-- 4a: Test how the 100m model performs at the regional scale
 
# First separate the Obs testing data by region. 
test.regions <- Partition.Test.Data( point.data$Obs[ point.data$Obs$TestDat == 1, ] )

# Test results tabulated in the RMD file.


#------------------------------------------------------------------
#-- 4b: Test how the two model resolutions compare across depths.
#     Assess ribbons within each region

# First test the 100 m model by ribbon within each region. 
# Then test the regional models; each model uses its own 20 m Obs testing data.
# test.regions are the 100 m test data partitioned above. 
#-- Also writes CSV of results.

depth.results <- Models.Across.Ribbons( test.regions )


#--------------------------------------------
# Build the Heat maps 
#  Needs to be done here cuz behaviour is odd in the RMD file ... 


Heat.Build.Class.Stats( build.sum$build.results.ByClass, T, 'TPR', rev( pal.heat.10 ), 800, 600, 'black' )
Heat.Build.Class.Stats( build.sum$build.results.ByClass, T, 'TNR', rev( pal.heat.10 ), 800, 600, 'black' )
Heat.Build.Class.Stats( build.sum$build.results.ByClass, T, 'User', rev( pal.heat.10 ), 800, 600, 'black' )

Heat.Build.Class.Stats( build.sum.nw$build.results.ByClass, F, 'TPR', rev( pal.heat.10 ), 800, 600, 'black' )
Heat.Build.Class.Stats( build.sum.nw$build.results.ByClass, F, 'TNR', rev( pal.heat.10 ), 800, 600, 'black' )
Heat.Build.Class.Stats( build.sum.nw$build.results.ByClass, F, 'User', rev( pal.heat.10 ), 800, 600, 'black' )


Heat.Build.Var.Import( build.sum$build.results.VarImport, pal.heat.11, 1000, 600, 'black' )


#--------------------------------------------------------------------
#-- Some regression to see what's driving the different metrics  ... 
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
IDE.results.trm  <- IDS.Evaluation( 'trm' )
IDE.results.nowt <- IDS.Evaluation( 'nwrf' )


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



# This function can be used to generate various views of the data. 
# 2020/05/25 DH preferred IDS as the grouping variable ... 
# 2020/08/23 Standardized them all as the difference from random baseline
Plot.Stats.By.IDS.For.Regions( z, pal.cb3b, sz = 25, lx=0.85, ly=0.87 )


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


#-- Assess BoP performance on IDS data by region.
# 2020/07/22: Some ugliness in this function including
#   Hard-coding of the region, needed to overcome the factoring of the names
#   Regions w low IDE sample size have 0 values in some classes
#   BoPs are drawn from an unrelated directory

if ( boptests == F ){
  
# Pre-build QCS and SOG ... 
  a <-  Build.IDE.BoP.Results('QCSSOG_BoPs_v1.0.gdb', 'BoP18_merged', 'Queen Charlotte Strait')
  a$Region <- substr( a$Region, 1, 3)  
  
  b <-  Build.IDE.BoP.Results('QCSSOG_BoPs_v1.0.gdb', 'BoP18_merged', 'Strait of Georgia')
  b$Region <- substr( b$Region, 4, 6)  

# Now combine with the other regions ...   
  IDE.BoP <- rbind( 
    Build.IDE.BoP.Results('HG_BoPs_v2.gdb',    'BoP18_merged_again', 'Haida Gwaii'),
    Build.IDE.BoP.Results('NCC_BoPs_v1.1.gdb', 'BoP18_merged', 'North Central Coast'),
    Build.IDE.BoP.Results('WCVI_BoPs_v1.gdb', 'BoPs', 'West Coast Vancouver Island'),
    a, b
  )  

  # Sort rows by IDS first then Region ... 
  IDE.BoP <- IDE.BoP[ order( IDE.BoP$IDS ), ]
  rownames( IDE.BoP ) <- NULL

  # NOTE: Regional name corrections for 'QCSSOG_BoPs_v1.0.gdb' done manually. 
  
  x <- substr( Sys.time(), 1, 10)
  save( IDE.BoP, file = file.path( model.dir, paste0('BoP_test_', x, '.RData')) )
    
} else {
  
  load( file.path( model.dir, 'BoP_test_2020-09-01.RData' ))
}
  





##################################################################
#------------ rMarkdown separation complete to here ------------


  
  
  
  
  

#=============================================
#-- Part 7 - Spatial substrate predictions ... 
#-- Create study area-wide predictions so that can compare prevalence of 
# study area prediction to training data ... 


#----- Plot the Study Area prevalence results ------

# Prepare the prevalence tables for each region, knock into shape, and ggplot() ... 
# Called by RMD: 
#   Plot.Pred.Map.Prevalence( map.prev, build.sum)


#---- Manual plotting of predicted tifs ----
plotme<-F

if (plotme == TRUE) {
  a <- raster( file.path( raster.dir, 'WCVI_classified_substrate.tif' ))

  #-- This works ok, but image has much poorer resolution than TIF loaded into ArcGIS ... 
  #   Try bumping up the bits ... 
  
  # original parameters ... 
  #   height = 7, width = 6, units = "in", res = 400)
  #   plot(a, maxpixels=5000000, col=pal.RMSM, legend=FALSE,
       
  # Map (up to 5,000,000 pixels)
  png( file=file.path( raster.dir, "TIF test.png"),
      height = 14, width = 12, units = "in", res = 600)
  plot(a, maxpixels=10000000, col=pal.RMSM, legend=FALSE,
       xlab = "Easting", ylab = "Northing", cex.axis = .5, cex.lab = .75)
  legend(x = "bottomleft",
         legend = c("Rock", "Mixed", "Sand", "Mud"), 
         fill = pal, title= NA, bg = NA, box.col = NA)
  
  dev.off()
}
#-- FIN.
