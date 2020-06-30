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
#source( "depth_effect.R" )
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


#-------------------------------------------
#  Load the predicted rasters and prevalences (map.pred) ... 
load( file.path( model.dir, 'rasterMapObjects_2020-06-29.RData' ))


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

# Do it again with the longer list of ribbons.
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


#-- IDS test results across depth ribbons by region. 
# Used to plot of Pontius stats 
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


#-- Assess BoP performance on IDS for 3(4) regions w sufficient data
IDE.BoP <- rbind( 
  Build.IDE.BoP.Results('HG_BoPs_v2.gdb',    'BoP18_merged_again'),
  Build.IDE.BoP.Results('NCC_BoPs_v1.1.gdb', 'BoP18_merged')
)  
rownames( IDE.BoP ) <- NULL

# Build.IDE.BoP.Results('QCSSOG_BoPs_v1.0.gdb', 'BoP18_merged')
#--> Need to either separate the BoPs or select BoP pgons for each region. Ugh.



#------------ rMarkdown separation complete to here ------------

# IDE still a bit of a work in progress. 


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




#----- Plot the Study Area prevalence results ------

# Prepare the prevalence tables for each region, knock into shape, and ggplot() ... 
# Called by RMD: 
#   Plot.Pred.Map.Prevalence( map.prev, build.sum)


#---- Manual plotting of predicted tifs ----
#  UNDER DEVELOPMENT - Saving and plotting predictions has become a bit of a pain. 

#-- FIN.
