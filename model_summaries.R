#*******************************************************************************
# Script:  model_summaries.R
# Created: 01 April 2020. EJG
#
# Purpose: Create a report on the model results.
#   Currently holds out-of-date code for writing to a text file. Update to markdown?

#*******************************************************************************



#====================================================================================
#-- Independent Data Evaluation --
# Comparisons include:
#   1. Coastwide model vs. all 3 ID sets.
#     ID merged from the regions to test coastwide.
#     summarizes the merged IDS
#   2. Regional models vs. all 3 ID sets. 
#     summarizes the regional IDS
# Build a table to hold the performance scores for the various comparisons.
# All models have been parameterized with the full set of training observations.
Do.Independent.Evaluation <- function( results=NULL ){
#-- Part 1: Coast model vs. each ID set. Regional IDE sets need to be assembled.  
  
  # Dive Data - pull the regional data together  ... 
  x.test <- Assemble.Coast.Test( dive.data )
  
  #-- Build the table piece ... 
  compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = 'Dive' )
  w <- Results.Row( rf.coast, x.test )
  x <- cbind( compare.what, w$Integrated )
  results <- rbind( results, x ) 
  
  # Repeat for the Cam IDS  ... 
  x.test <- Assemble.Coast.Test( cam.data )
  
  #-- Build the table piece ... 
  compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = 'Cam' )
  w <- Results.Row( rf.coast, x.test )
  x <- cbind( compare.what, w$Integrated )
  results <- rbind( results, x ) 
  
  # Repeat for the ROV IDS  ... 
  x.test <- Assemble.Coast.Test( ROV.data )
  
  #-- Build the table piece ... 
  compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = 'ROV' )
  w <- Results.Row( rf.coast, x.test )
  x <- cbind( compare.what, w$Integrated )
  results <- rbind( results, x ) 
  
  #------------------------------------------------------------
  #-- Part 2: Regional runs, each with all 3 ID sets.
  #   THIS IS 5 x 3 = 15 comparisions ... with ROV exception.
  #   Requres: Models to be loaded (rf.region.XX)
  
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
        results <- rbind( results, x ) 
      }
    }
  }
  return( results )
}



#-------------------------------------------------------------------------------
# Construct summary tables of build results. 
#   Writes to CSV files. Returns nothing. 
#     Build_results_Integrated.csv
#     Build_results_byClassStats.csv
#     Build_results_test_ClassPrevalence.csv
#     Build_results_varImportance.csv
#
Summarize.Build <- function( build.df ){
# Bunch of hacking here to pull the tables together ... 
# Requires: build.results data.frame as input.
#   Includes pairs of lists for each run. 
#   First is list of Integrated and perClass stats; Second is list of variable importance
# ASSUMES 6 regions done.
  
  # Build a list of names ... 
  a <- names( build.df )[c(1,3,5,7,9,11)]
  nm <- unlist( lapply( a, strsplit, split = '\\.' ))[c(1,3,5,7,9,11)]
  
  # Pull Integrated Stats ... number are the rows for each region-specific stat
  x <- do.call(rbind.data.frame, 
               build.df[ c(1,3,5,7,9,11) ] )
  
  # x has 2 components, the Integrated bit, .... 
  y <- do.call( rbind.data.frame, x$Integrated )
  
  z <- cbind( 'Model' = nm, y )
  row.names(z) <- NULL
  
  out.file <- 'Build_results_Integrated.csv'
  write.csv( z, file = file.path(output.dir, out.file) )
  
  # and the PerClass bit which includes:
  #   1) A table with both User and Producer. A df with 2 sets of rows.
  y <- do.call( rbind,  x$PerClass )
  z <- y[ (row.names(y) == 'User') | (row.names(y) == 'Prod') , ]
  
  # User first ... 
  y.usr <- data.frame( 'Region' = nm, z[ c(1,3,5,7,9,11), ])
  row.names(y.usr) <- NULL
  
  # Now producer  
  y.prd <- data.frame( 'Region' = nm, z[ c(2,4,6,8,10,12), ])
  row.names(y.prd) <- NULL
  
  zz <- rbind(
    cbind( 'Stat' = 'User', y.usr ),
    cbind( 'Stat' = 'Prod', y.prd ) )
  colnames(zz) <- c('Stat','Region','Hard','Mixed','Sand','Mud')
  out.file <- 'Build_results_byClassStats.csv'
  write.csv( zz, file = file.path(output.dir, out.file) )
  
 #   2) The prevalence of the testing obs vs. predicted.
  z <- y[ (row.names(y) == 'PrevObs') | (row.names(y) == 'PrevPred') , ]
  
  # Obs first ... 
  y.obs <- data.frame( 'Region' = nm, z[ c(1,3,5,7,9,11), ])
  row.names(y.obs) <- NULL
  
  # Now predicted  
  y.prd <- data.frame( 'Region' = nm, z[ c(2,4,6,8,10,12), ])
  row.names(y.prd) <- NULL
  
  zz <- rbind(
    cbind( 'Stat' = 'Obs', y.obs ),
    cbind( 'Stat' = 'Pred', y.prd ) )
  colnames(zz) <- c('Stat','Region','Hard','Mixed','Sand','Mud')
  out.file <- 'Build_results_test_ClassPrevalence.csv'
  write.csv( zz, file = file.path(output.dir, out.file) )
  
  #-- And finally variable importance ... 
  # Integrated stats first ... 
  # 2020/04/09: Moved the ranking to the plot function so values can go thru and be displayed 
  
  y <- data.frame(
    rbind( as.vector( build.df$Coast.Import ),
           as.vector( build.df$HG.Import    ),
           as.vector( build.df$NCC.Import   ),
           as.vector( build.df$WCVI.Import  ),
           as.vector( build.df$QCS.Import   ),
           as.vector( build.df$SOG.Import   )
    ))
  
  row.names(y) <- nm
  colnames(y)  <- names( build.df$HG.Import )
  
  out.file <- 'Build_results_varImportance.csv'
  write.csv( y, file = file.path(output.dir, out.file) )
}




###################################
### OLD CODE ####

#---Build Summary -----------------------------------------------------------------------
#-- Mar31: out-of-date code for writing Build results to a text file. Update to markdown?

# log.file <- 'build_ranger_models_log.txt'
# sink( file = file.path(output.dir, log.file))
# 
# cat('Ranger Model Building Results\n')
# cat('Run Date: ', date(), '\n')
# cat('-----------------------------\n')
# cat('Train partition: ', tpart, '\n')
# cat('Iterations:      ', iter,  '\n')
# 
# cat('\n------------------------------\n')
# cat('100 m Coastwide - Data Summary\n')
# x <- train.data.100m
# 
# x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
# cat( 'N =         ', nrow(x), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')
# 
# sink()
# foo <- Make.Ranger.Model( x, coast.formula, tpart, iter )
# imp <- foo$Model$variable.importance
# 
# #ASSIGN and REPORT on MODEL ... 
# rf.coast <- foo$Model
# build.results <- c( build.results, 'Coast' = list( 'stats' = foo$Stats, 'Import' = imp) )
# 
# sink( file = file.path(output.dir, log.file), append=TRUE)
# cat("-- Model fit --\n")
# foo$Stats[[1]]  #-- Integrated
# cat("\n--- by class ---\n")
# foo$Stats[[2]]
# cat("\n-- Variable importance --\n")
# sort( imp[ imp >= median( imp ) ], decreasing = TRUE)
# 
# 
# cat('\n-------------------- -\n')
# cat('20 m HG - Data Summary\n')
# x <- train.data.20m$HG
# x$rugosity <- x$rugosity / max( x$rugosity )
# 
# x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
# cat( 'N =         ', nrow(x), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')
# 
# sink()
# foo <- Make.Ranger.Model( x, shore.formula, tpart, iter)
# imp <- foo$Model$variable.importance
# 
# #ASSIGN MODEL ... 
# rf.region.HG <- foo$Model
# build.results <- c( build.results, 'HG' = list( 'stats' = foo$Stats, 'Import' = imp) )
# 
# 
# sink( file = file.path(output.dir, log.file), append=TRUE)
# cat("-- Model fit --\n")
# foo$Stats[[1]]  #-- Integrated
# cat("\n--- by class ---\n")
# foo$Stats[[2]]
# cat("\n-- Variable importance --\n")
# sort( imp[ imp >= median( imp ) ], decreasing = TRUE)
# 
# 
# cat('\n-----------------------\n')
# cat('20 m NCC - Data Summary\n')
# x <- train.data.20m$NCC
# x$rugosity <- x$rugosity / max( x$rugosity )
# 
# x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
# cat( 'N =         ', nrow(x), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')
# 
# sink()
# foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
# imp <- foo$Model$variable.importance
# 
# #ASSIGN MODEL ... 
# rf.region.NCC <- foo$Model
# build.results <- c( build.results, 'NCC' = list( 'stats' = foo$Stats, 'Import' = imp) )
# 
# sink( file = file.path(output.dir, log.file), append=TRUE)
# cat("-- Model fit --\n")
# foo$Stats[[1]]  #-- Integrated
# cat("\n--- by class ---\n")
# foo$Stats[[2]]
# cat("\n-- Variable importance --\n")
# sort( imp[ imp >= median( imp ) ], decreasing = TRUE)
# 
# 
# cat('\n------------------------\n')
# cat('20 m WCVI - Data Summary\n')
# x <- train.data.20m$WCVI
# x$rugosity <- x$rugosity / max( x$rugosity )
# 
# x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
# cat( 'N =         ', nrow(x), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')
# 
# sink()
# foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
# imp <- foo$Model$variable.importance
# 
# #ASSIGN MODEL ... 
# rf.region.WCVI <- foo$Model
# build.results <- c( build.results, 'WCVI' = list( 'stats' = foo$Stats, 'Import' = imp) )
# 
# sink( file = file.path(output.dir, log.file), append=TRUE)
# cat("-- Model fit --\n")
# foo$Stats[[1]]  #-- Integrated
# cat("\n--- by class ---\n")
# foo$Stats[[2]]
# cat("\n-- Variable importance --\n")
# sort( imp[ imp >= median( imp ) ], decreasing = TRUE)
# 
# 
# cat('\n-----------------------\n')
# cat('20 m QCS - Data Summary\n')
# x <- train.data.20m$QCS
# x$rugosity <- x$rugosity / max( x$rugosity )
# 
# x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
# cat( 'N =         ', nrow(x), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')
# 
# sink()
# foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
# imp <- foo$Model$variable.importance
# 
# #ASSIGN MODEL ... 
# rf.region.QCS <- foo$Model
# build.results <- c( build.results, 'QCS' = list( 'stats' = foo$Stats, 'Import' = imp) )
# 
# sink( file = file.path(output.dir, log.file), append=TRUE)
# cat("-- Model fit --\n")
# foo$Stats[[1]]  #-- Integrated
# cat("\n--- by class ---\n")
# foo$Stats[[2]]
# cat("\n-- Variable importance --\n")
# sort( imp[ imp >= median( imp ) ], decreasing = TRUE)
# 
# 
# cat('\n-----------------------\n')
# cat('20 m SOG - Data Summary\n')
# x <- train.data.20m$SOG
# x$rugosity <- x$rugosity / max( x$rugosity )
# 
# x.pvec <- summary(x$BType4) / sum( summary(x$BType4) )
# cat( 'N =         ', nrow(x), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n\n')
# 
# sink()
# foo <- Make.Ranger.Model( x, shore.formula, tpart, iter )
# imp <- foo$Model$variable.importance
# 
# #ASSIGN MODEL ... 
# rf.region.SOG <- foo$Model
# build.results <- c( build.results, 'SOG' = list( 'stats' = foo$Stats, 'Import' = imp) )
# 
# sink( file = file.path(output.dir, log.file), append=TRUE)
# cat("-- Model fit --\n")
# foo$Stats[[1]]  #-- Integrated
# cat("\n--- by class ---\n")
# foo$Stats[[2]]
# cat("\n-- Variable importance --\n")
# sort( imp[ imp >= median( imp ) ], decreasing = TRUE)
# 
# sink()
# 
# 
# #---IDE Results Summary -----------------------------------------------------------------------
# #-- 2020/03/31: As above, out-of-date code for writing Independent data evaluation results to a text file.
# 
# #-- DON'T use this test unless you do it for the regions also. Its easily confused 
# #   with the build evaulation using the 30% test sample. 
# # #-- Coast model vs. 100 m train ...
# # compare.what <- data.frame( 'Model' = '100m coast', 'Test Data' = '100 m Train' )
# # w <- Results.Row( rf.coast, x.test )
# # x <- cbind( compare.what, w$Integrated )
# # results.table <- rbind( results.table, x ) 
# 
# log.file <- 'independent_data_evaluation_log.txt'
# 
# sink( file = file.path(output.dir, log.file))
# cat("Independent Data Evaluation \n")
# cat("Run Date: ", date(), "\n")
# cat("---------------------------\n\n")
# 
# #-- Summarize the different observational data sets.
# cat("Summary of the observational data\n")
# cat("---------------------------------\n\n")
# 
# cat( '100 m train data - Summary\n')
# x.test <- train.data.100m
# x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
# cat( 'N =         ', nrow(x.test), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')
# 
# cat("\n---------------------------------\n\n")
# cat("20 m Dive IDS - Summary\n")
# x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
# cat( 'N =         ', nrow(x.test), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')
# 
# cat("\n---------------------------------\n\n")
# cat("20 m Camera IDS - Summary\n")
# x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
# cat( 'N =         ', nrow(x.test), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')
# 
# cat("\n---------------------------------\n\n")
# cat("20 m ROV IDS - Summary\n")
# x.pvec <- summary(x.test$BType4) / sum( summary(x.test$BType4) )
# cat( 'N =         ', nrow(x.test), '\n')
# cat( 'Prevalence =', round( x.pvec, 4), '\n')
# cat( 'Imbalance  =', round( Prev.Balance( x.pvec ), 4), '\n')
# 
# cat("\n\n---------------------------------\n")
# cat("------ IDS evaluation Coastwide ---------\n\n")
# 
# row.names(results.table) <- NULL
# results.table
# sink()
# 
# 
# sink( file = file.path(output.dir, log.file), append=TRUE)
# cat("\n\n----------------------------------------------\n")
# cat("---------- IDS evaluation by Region ----------\n")
# cat("    Regional models built w 20 m data \n\n")
# results.table
# sink()
# 
# #
# 
