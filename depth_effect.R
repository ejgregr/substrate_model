#*******************************************************************************
# Script:  depth_effect.R
# Created: 19 March 2020. EJG
#
# Purpose: Evaluate a suite of different RF models with Independent data (ID) across depth classes. 
#   Part of the substrate family of scripts. 
#   Part 1 evaluates 100 m coastwide model using training data, and aggregated IDS for entire coast
#   Part 2 evaluates regional and coastwide models by region and data set
#
# Revised: 2020/04/06: Focus on two tests, one with obs testing data, the other with focused IDE. 
#
# Notes:
#  - Requires loaded models and observational data from main script.
# Results: Makes 2 csv files:
#   Depth_compare_coast.csv : Statistic basket for depth zones, coastwide
#   Depth_compare_regions.csv : Statistic basket for depth zones, by Region
#********************************************************************************

# Define the depth classes ... 
#z.breaks <- c( -5000, -50, -20, -10, -5, 0, 1000)
z.breaks <- c( -1000, 0, 5, 10, 20, 50, 5000)
z.ribs <- c('ITD', '0-5', '5-10', '10-20', '20-50', '50+')

#======================= Testing with withheld Obs =========================
#-- Compare 100 m model and 20 m regional models, across depth zones
#   to see if there is a depth-resolution linkage.

#-----------------------------------
# Returns a row of test results ... 
# GLOBAL vars: z.breaks, z.ribs, rf models, obs.20mIV
Model.Fit.Obs.By.Depth <- function( regName, mant = 3 ){

  rf <- eval(parse( text=paste0( "rf.region.", regName ) ))
  
  if (regName == 'Coast') {
    tdat <- obs.100mIV
  } else
    tdat <- eval(parse( text=paste0( 'obs.20mIV$', regName ) ))

  tdat <- tdat[ tdat$TestDat == 1, ]
  tdat$zClass <- as.factor( findInterval( tdat$bathy, z.breaks) )

  results <- NULL  
  # loop thru each level, calculating performance of data subset 
  for (k in levels( tdat$zClass ) ){
    x.sub <- tdat[ tdat$zClass == k, ]

    # build the row label ... 
    compare.what <- data.frame( 'Ribbon' = z.ribs[ as.numeric(k) ], 
                                'meanZ' = round( mean( x.sub$bathy ), mant ) )
    
    # make the results ... 
    w <- Results.Row( rf, x.sub )
    x <- cbind( compare.what, w$Integrated )
    results <- rbind( results, x ) 
  }
  return (results)
}

#-----------------------------------
#-- For the 100 m model, how does it perform (using withheld obs) by Depth within Region?
# Model is fixed. 
# INPUT: Obs test data, partitioned by region ...  
Coast.Fit.By.Region.By.Depth <- function( tdat, mant = 3 ){
  
  rf <- rf.region.Coast
  
  tdat$zClass <- as.factor( findInterval( tdat$bathy, z.breaks) )
  
  results <- NULL  
  # loop thru each level, calculating performance of data subset 
  for (k in levels( tdat$zClass ) ){
    x.sub <- tdat[ tdat$zClass == k, ]
    
    # build the row label ... 
    compare.what <- data.frame( 'Ribbon' = z.ribs[ as.numeric(k) ], 
                                'meanZ' = round( mean( x.sub$bathy ), mant ) )
    
    # make the results ... 
    w <- Results.Row( rf, x.sub )
    x <- cbind( compare.what, w$Integrated )
    results <- rbind( results, x ) 
  }
  return (results)
}



#----- Duplicate pair of functions to extend depth analysis to a couple more depth classes.
# Define 2nd set of breaks to test depth effect more thoroughly ... 
z.breaks2 <- c( -1000, 0, 5, 10, 20, 50, 100, 200, 5000)
z.ribs2 <- c('ITD', '0-5', '5-10', '10-20', '20-50', '50-100', '100-200', '200+')

Model.Fit.Obs.By.Depth2 <- function( regName, mant = 3 ){
  
  rf <- eval(parse( text=paste0( "rf.region.", regName ) ))
  
  if (regName == 'Coast') {
    tdat <- obs.100mIV
  } else
    tdat <- eval(parse( text=paste0( 'obs.20mIV$', regName ) ))
  
  tdat <- tdat[ tdat$TestDat == 1, ]
  tdat$zClass <- as.factor( findInterval( tdat$bathy, z.breaks2) )
  
  results <- NULL  
  # loop thru each level, calculating performance of data subset 
  for (k in levels( tdat$zClass ) ){
    x.sub <- tdat[ tdat$zClass == k, ]
    
    # build the row label ... 
    compare.what <- data.frame( 'Ribbon' = z.ribs2[ as.numeric(k) ], 
                                'meanZ' = round( mean( x.sub$bathy ), mant ) )
    
    # make the results ... 
    w <- Results.Row( rf, x.sub )
    x <- cbind( compare.what, w$Integrated )
    results <- rbind( results, x ) 
  }
  return (results)
}

Coast.Fit.By.Region.By.Depth2 <- function( tdat, mant = 3 ){
  
  rf <- rf.region.Coast
  
  tdat$zClass <- as.factor( findInterval( tdat$bathy, z.breaks2) )
  
  results <- NULL  
  # loop thru each level, calculating performance of data subset 
  for (k in levels( tdat$zClass ) ){
    x.sub <- tdat[ tdat$zClass == k, ]
    
    # build the row label ... 
    compare.what <- data.frame( 'Ribbon' = z.ribs2[ as.numeric(k) ], 
                                'meanZ' = round( mean( x.sub$bathy ), mant ) )
    
    # make the results ... 
    w <- Results.Row( rf, x.sub )
    x <- cbind( compare.what, w$Integrated )
    results <- rbind( results, x ) 
  }
  return (results)
}


#========================== IDS Tests by depth ==============================
# 2020/05/10: Below needs updating based on objectives.

results.table <- NULL   # Re-use this name for the results here. 
z.table <- NULL         # A mean depth table, FWIW.

#---------------------------------------------------
# 1. DIVE DATA: Assemble the data across regions ... 
rm( 'x.sub', 'x.test')
x.test <- Assemble.NCoast.Test( dive.20mIV )
j <- "Dive"
z.row <- NULL

# Add the depthClass to the test data ... 
#-- what is this? - its [0,1,2,3]
# hist(x.test$DepthCat) 
x.test$zRibbon <- as.factor( 6 - findInterval( x.test$bathy, z.breaks) )

for (k in levels(x.test$zRibbon) ){
  
  x.sub <- x.test[ x.test$zRibbon == k, ]
  z.row <- cbind( z.row, round( mean(x.sub$bathy), 2) )
  
  # build the row label ... 
  compare.what <- data.frame( 'Test Data' = j, 'Ribbon' = z.ribs[ as.numeric(k) + 1 ] )
  
  # make the results ... 
  w <- Results.Row( m.model, x.sub )
  x <- cbind( compare.what, w$Integrated )
  results.table <- rbind( results.table, x ) 
}
z.row <- data.frame (z.row )
colnames( z.row ) <- z.ribs
z.row <- cbind('IDS' = j, z.row)
z.table <- rbind( z.table, z.row)

#---------------------------------------------------
# 2. CAM DATA: Assemble the data across regions ... 
rm( 'x.sub', 'x.test')
x.test <- Assemble.Coast.Test( cam.data )
i <- 1
j <- 'Cam'
z.row <- NULL

# Add the depthClass to the test data ... 
x.test$zRibbon <- as.factor( 6 - findInterval( x.test$bathy, z.breaks) )

for (k in levels(x.test$zRibbon) ){
  
  x.sub <- x.test[ x.test$zRibbon == k, ]
  z.row <- cbind( z.row, round( mean(x.sub$bathy), 2) )
  # build the row label ... 
  compare.what <- data.frame( 'Test Data' = j, 'Ribbon' = z.ribs[ as.numeric(k) + 1 ] )
  
  # make the results ... 
  w <- Results.Row( m.model, x.sub )
  x <- cbind( compare.what, w$Integrated )
  results.table <- rbind( results.table, x ) 
  i <- i+1
}
z.row <- data.frame (z.row )
colnames( z.row ) <- z.ribs
z.row <- cbind('IDS' = j, z.row)
z.table <- rbind( z.table, z.row)

#---------------------------------------------------
# 3. ROV DATA: Assemble the data across regions ... 
rm( 'x.sub', 'x.test')
x.test <- Assemble.Coast.Test( ROV.data )
j <- 'ROV'
z.row <- NULL

# Add the depthClass to the test data ... 
x.test$zRibbon <- as.factor( 6 - findInterval( x.test$bathy, z.breaks) )

for (k in levels(x.test$zRibbon) ){
  
  x.sub <- x.test[ x.test$zRibbon == k, ]
  z.row <- cbind( z.row, round( mean(x.sub$bathy), 2) )
  # build the row label ... 
  compare.what <- data.frame( 'Test Data' = j, 'Ribbon' = z.ribs[ as.numeric(k) + 1 ] )
  
  # make the results ... 
  w <- Results.Row( m.model, x.sub )
  x <- cbind( compare.what, w$Integrated )
  results.table <- rbind( results.table, x ) 
}

## MASSIVE HACK on the CBIND here to deal w missing ROV data... ###
z.row <- data.frame ( cbind( 0, 0, 0, z.row) )
colnames( z.row ) <- z.ribs
z.row <- cbind('IDS' = j, z.row)
z.table <- rbind( z.table, z.row)


#-------------------------------
# Report results ... 

cat("\n------------------------------------\n")
cat("Depth Ribbon Comparisons - Coastwide\n\n")

cat("Mean depths by ribbon\n")
z.table 

cat("\n\nStats by IDS and ribbon\n")
results.table

out.file <- 'Depth_compare_coast.csv'
write.csv( results.table, file = file.path(output.dir, out.file) )



#=====================================================================
#-- Part 2: Evaluation of the regional models ... 
#   zRibbons are on [0, 5], representing [ITD, 50+]

results.table <- NULL

IDS <- c("dive", "cam", "ROV")

for (i in bioregions) {
  # select the regional model ...
  m.model <- eval(parse( text=paste0( "rf.region.", i ) ))
  cat(i, "\n")
  
  for (j in IDS ){
  # Select the evaluation data set ... 
    cat("  ",j, "\n")  
    
    if (j == "ROV" && i %in% c('WCVI', 'QCS', 'SOG') ){
      # Bail ... 
      print( "skipping ... ")
    } else { # good to go ... 
      
      # grab the IDS data ... 
      x <- eval(parse( text=paste0( j, ".data" ) ))
      x.test <- x[[ i ]]
      
      # Add the depthClass to the test data ... 
      x.test$zRibbon <- as.factor( 6 - findInterval( x.test$bathy, z.breaks) )
      
      for (k in levels(x.test$zRibbon) ){
        
        x.sub <- x.test[ x.test$zRibbon == k, ]
        z.row <- cbind( z.row, round( mean(x.sub$bathy), 2) )
        # build the row label ... 
        compare.what <- data.frame( 'Model' = i, 'Test Data' = j, 'Ribbon' = z.ribs[ as.numeric(k) + 1 ] )
        
        # make the results ... 
        w <- Results.Row( m.model, x.sub )
        x <- cbind( compare.what, w$Integrated )
        results.table <- rbind( results.table, x ) 
      }
    }
  }
}

out.file <- 'Depth_compare_regions.csv'
write.csv( results.table, file = file.path(output.dir, out.file) )
#results.table <- read.csv(file = file.path(output.dir, 'Depth_compare_regions.csv') )




#-- FIN.