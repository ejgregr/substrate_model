#*******************************************************************************
# Script:  IDE_main.R
# Purpose: Evaluate wtd RF models with Independent data (ID) across various data partitions
#   Part 1: Load all the observational data
#   Part 2: Build the necessary RF models
#   Part 3: Compare the various models and ID sets 
# Created: January 2020. EJG
# Updated: 2020/03/10. EJG
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
# obs.data - training data
# dive.data, cam.data, ROV.data = Independent data sets
ls()

#-- OR Build new ... can take 20-30 min.
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

#-- Load the latest RF models: Last update: 2020/03/18.
load( 'C:\\Users\\Edward\\Dropbox\\DFO Job\\Substrate2019\\Models\\rf_allModels_2020-03-18-1558.RData')
ls()

#-- Or build them ... 
starttime <- Sys.time()

#=======================================================================================
#=== BUILD AND SUMMARIZE 6 RF MODELS 

tpart <- 0.7
iter  <- 1

#-- Need to save the performance stats in something ... 
#   Use a list of 3 results for each model.
#   Summary of variable importance currently only reports those with higher than median values

build.results <- list()


log.file <- 'build_ranger_models_log.txt'
sink( file = file.path(output.dir, log.file))

cat('Ranger Model Building Results\n')
cat('Run Date: ', date(), '\n')
cat('-----------------------------\n')
cat('Train partition: ', tpart, '\n')
cat('Iterations:      ', iter,  '\n')

cat('\n------------------------------\n')
cat('100 m Coastwide - Data Summary\n')
x <- train.data.100m

x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
cat( 'N =         ', nrow(x), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')

sink()
foo <- Make.Ranger.Model( x, coast.formula, tpart, iter )
imp <- foo$Model$variable.importance

#ASSIGN and REPORT on MODEL ... 
rf.coast <- foo$Model
build.results <- c( build.results, 'Coast' = list( 'stats' = foo$Stats, 'Import' = imp) )

sink( file = file.path(output.dir, log.file), append=TRUE)
cat("-- Model fit --\n")
foo$Stats[[1]]  #-- Integrated
cat("\n--- by class ---\n")
foo$Stats[[2]]
cat("\n-- Variable importance --\n")
sort( imp[ imp >= median( imp ) ], decreasing = TRUE)


cat('\n-------------------- -\n')
cat('20 m HG - Data Summary\n')
x <- train.data.20m$HG
x$rugosity <- x$rugosity / max( x$rugosity )

x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
cat( 'N =         ', nrow(x), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')

sink()
foo <- Make.Ranger.Model( x, shore.formula, tpart, iter)
imp <- foo$Model$variable.importance

#ASSIGN MODEL ... 
rf.region.HG <- foo$Model
build.results <- c( build.results, 'HG' = list( 'stats' = foo$Stats, 'Import' = imp) )


sink( file = file.path(output.dir, log.file), append=TRUE)
cat("-- Model fit --\n")
foo$Stats[[1]]  #-- Integrated
cat("\n--- by class ---\n")
foo$Stats[[2]]
cat("\n-- Variable importance --\n")
sort( imp[ imp >= median( imp ) ], decreasing = TRUE)


cat('\n-----------------------\n')
cat('20 m NCC - Data Summary\n')
x <- train.data.20m$NCC
x$rugosity <- x$rugosity / max( x$rugosity )

x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
cat( 'N =         ', nrow(x), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')

sink()
foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
imp <- foo$Model$variable.importance

#ASSIGN MODEL ... 
rf.region.NCC <- foo$Model
build.results <- c( build.results, 'NCC' = list( 'stats' = foo$Stats, 'Import' = imp) )

sink( file = file.path(output.dir, log.file), append=TRUE)
cat("-- Model fit --\n")
foo$Stats[[1]]  #-- Integrated
cat("\n--- by class ---\n")
foo$Stats[[2]]
cat("\n-- Variable importance --\n")
sort( imp[ imp >= median( imp ) ], decreasing = TRUE)


cat('\n------------------------\n')
cat('20 m WCVI - Data Summary\n')
x <- train.data.20m$WCVI
x$rugosity <- x$rugosity / max( x$rugosity )

x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
cat( 'N =         ', nrow(x), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')

sink()
foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
imp <- foo$Model$variable.importance

#ASSIGN MODEL ... 
rf.region.WCVI <- foo$Model
build.results <- c( build.results, 'WCVI' = list( 'stats' = foo$Stats, 'Import' = imp) )

sink( file = file.path(output.dir, log.file), append=TRUE)
cat("-- Model fit --\n")
foo$Stats[[1]]  #-- Integrated
cat("\n--- by class ---\n")
foo$Stats[[2]]
cat("\n-- Variable importance --\n")
sort( imp[ imp >= median( imp ) ], decreasing = TRUE)


cat('\n-----------------------\n')
cat('20 m QCS - Data Summary\n')
x <- train.data.20m$QCS
x$rugosity <- x$rugosity / max( x$rugosity )

x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
cat( 'N =         ', nrow(x), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')

sink()
foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
imp <- foo$Model$variable.importance

#ASSIGN MODEL ... 
rf.region.QCS <- foo$Model
build.results <- c( build.results, 'QCS' = list( 'stats' = foo$Stats, 'Import' = imp) )

sink( file = file.path(output.dir, log.file), append=TRUE)
cat("-- Model fit --\n")
foo$Stats[[1]]  #-- Integrated
cat("\n--- by class ---\n")
foo$Stats[[2]]
cat("\n-- Variable importance --\n")
sort( imp[ imp >= median( imp ) ], decreasing = TRUE)


cat('\n-----------------------\n')
cat('20 m SOG - Data Summary\n')
x <- train.data.20m$SOG
x$rugosity <- x$rugosity / max( x$rugosity )

x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
cat( 'N =         ', nrow(x), '\n')
cat( 'Prevalence =', round( x.pvec, 4), '\n')
cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')

sink()
foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
imp <- foo$Model$variable.importance

#ASSIGN MODEL ... 
rf.region.SOG <- foo$Model
build.results <- c( build.results, 'SOG' = list( 'stats' = foo$Stats, 'Import' = imp) )

sink( file = file.path(output.dir, log.file), append=TRUE)
cat("-- Model fit --\n")
foo$Stats[[1]]  #-- Integrated
cat("\n--- by class ---\n")
foo$Stats[[2]]
cat("\n-- Variable importance --\n")
sort( imp[ imp >= median( imp ) ], decreasing = TRUE)

sink()


#-----------------------------
#-- Done building ranger RF models. 

endtime <- Sys.time()
cat('---------------------------\n')
cat( 'Model build time: ', endtime - starttime, '\n\n' )
#last run for 1 iteration was 7.4 minutes.

#-- SAVE The resulting models. Takes minutes - its a big file. 

# Build a time stamp ... 
x <- substr( endtime, 1, 16); x <- gsub(":", "", x); x <- gsub(" ", "-", x)

save( rf.coast, rf.region.HG, rf.region.NCC, 
      rf.region.WCVI, rf.region.QCS, rf.region.SOG, 
      file = file.path( model.dir, paste0('rf_allModels_', x, '.RData')) )

save( build.results, 
      file = file.path( model.dir, paste0('modelResults', x, '.RData')) )

#-- Can't just save this to file cuz its a complicated structure. So ... 


#-------------------------------------------------------------------------------
#-- Build some summary tables ... 
#   build.results includes a pair of list entries for each run. 
#   First of pair includes Integrated and perClass values;
#   Second of pair is a list of variable importance

#-- Bunch of hacking here to pull the tables together ... 
#   ASSUMES 6 models done, as above.

# Build a list of names ... 
a <- names( build.results )[c(1,3,5,7,9,11)]
nm <- unlist( lapply( a, strsplit, split = '\\.' ))[c(1,3,5,7,9,11)]
nm

# Integrated stats first ... 
x <- do.call(rbind.data.frame, 
  build.results[ c(1,3,5,7,9,11) ] )

y <- do.call( rbind.data.frame, x$Integrated )

z <- cbind( 'Model' = nm, y )
row.names(z) <- NULL
z

out.file <- 'Build_results_Integrated.csv'
write.csv( z, file = file.path(output.dir, out.file) )


# Now Class-based stats  ... Have a data.frame with 2 rows .... 
# User first 
x <- do.call(rbind.data.frame, 
             build.results[ c(1,3,5,7,9,11) ] )

y <- do.call( rbind.data.frame,  x$PerClass )

y.usr <- cbind( 'Model' = nm, y[ c(1,3,5,7,9,11), ])
row.names(y.usr) <- NULL

# Now producer  
y.prd <- cbind( 'Model' = nm, y[ c(2,4,6,8,10,12), ])
row.names(y.prd) <- NULL

zz <- rbind(
cbind( 'Stat' = 'User', y.usr ),
cbind( 'Stat' = 'Prod', y.prd ) )


out.file <- 'Build_results_byClassStats.csv'
write.csv( zz, file = file.path(output.dir, out.file) )


#-- And finally variable importance ... 
# Integrated stats first ... 

y <- data.frame(
rbind( c( rank( -as.vector( build.results$Coast.Import )), NA),
       rank( -as.vector( build.results$HG.Import    )),
       rank( -as.vector( build.results$NCC.Import   )),
       rank( -as.vector( build.results$WCVI.Import  )),
       rank( -as.vector( build.results$QCS.Import   )),
       rank( -as.vector( build.results$SOG.Import   ))
))

#-- order according to the results ... 
olist <- c(1,11,2,10,3,9,6,4,5,7,8)
y <- y[, olist]
tab <- as.matrix(y)
tab[1,2] <- 0

row.names(y) <- nm
colnames(y)  <- names( build.results$HG.Import )[ olist ]
y

out.file <- 'Build_results_varImportance.csv'
write.csv( y, file = file.path(output.dir, out.file) )


#----------------- Heat maps of the tables?
#-- Install colour palette ... 
devtools::install_github("jakelawlor/PNWColors") 
library(PNWColors)

#-- Install superheat ... 
devtools::install_github("rlbarter/superheat")
library(superheat)

pal <- pnw_palette( 'Bay', 10 )

# 1) User accuracy  ... 
foo <- y.usr[,-1]
row.names( foo ) <- y.prd$Model
tab <- round( as.matrix(foo), 2)
superheat( foo, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
           X.text = tab, X.text.col = 'lightblue', X.text.size = 10,
           order.rows = rev(1:nrow(foo)),
           left.label.col = 'white', left.label.text.size = 10,
           bottom.label.col = 'White', bottom.label.size = 0.1, bottom.label.text.size = 10
)

# 2) Producer accuracy  ... 
foo <- y.prd[,-1]
row.names( foo ) <- y.prd$Model
tab <- round( as.matrix(foo), 2)
superheat( foo, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
           X.text = tab, X.text.col = 'lightblue', X.text.size = 10,
           order.rows = rev(1:nrow(foo)),
           left.label.col = 'white', left.label.text.size = 10,
           bottom.label.col = 'White', bottom.label.size = 0.1, bottom.label.text.size = 10
)

# 3) Predictor importance ... 
pal <- rev(pnw_palette( 'Bay', 11 ))

superheat( y, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
           X.text = tab, X.text.col = 'lightblue', X.text.size = 10,
           order.rows = rev(1:nrow(y)),
           left.label.col = 'white', left.label.text.size = 10,
           bottom.label.col = 'White', bottom.label.text.angle = 55, 
           bottom.label.size = 0.25, bottom.label.text.size = 10, bottom.label.text.alignment = 'right'
           )


#====================================================================================
#-- Independent Data Evaluation --
# Comparisons include:
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

#-- DON'T do this unless you do it for the regions also. Its easily confused 
#   with the build evaulation using the 30% test sample. 
# #-- Coast model vs. 100 m train ...
# compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = '100 m Train' )
# w <- Results.Row( rf.coast, x.test )
# x <- cbind( compare.what, w$Integrated )
# results.table <- rbind( results.table, x ) 


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
w <- Results.Row( rf.coast, x.test )
x <- cbind( compare.what, w$Integrated )
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
w <- Results.Row( rf.coast, x.test )
x <- cbind( compare.what, w$Integrated )
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
w <- Results.Row( rf.coast, x.test )
x <- cbind( compare.what, w$Integrated )
results.table <- rbind( results.table, x ) 


cat("\n\n---------------------------------\n")
cat("------ IDS evaluation Coastwide ---------\n\n")

row.names(results.table) <- NULL
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
        w <- Results.Row( m.model, x.test )
        x <- cbind( compare.what, w$Integrated )
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


out.file <- 'IDS_evaluation_by_region.csv'
write.csv( results.table, file = file.path(output.dir, out.file) )



#======================================================================================
# Examination of depth effect moved to its own script - 2020/03/19.




#-- FIN.
