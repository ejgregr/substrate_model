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


#===========================================================================
#-- Part 1: Compare 100 m model to 20 m regional models, across depth zones
#   to see if there is a depth- resolution linkage. uses withheld obs data. 

results.table <- NULL
z.table <- NULL

# Define the depth classes ... 
z.breaks <- c( -5000, -50, -20, -10, -5, 0, 1000)
z.ribs <- c('ITD', "0-5", "5-10", "10-20", "20-50", "50+")

#---------------------------------------------------
# Grab 0. TRAINING DATA: 
rm( 'x.sub', 'x.test')
x.test <- train.data.100m
j <- "Train"
z.row <- NULL

#- confirm that depths have same sign ... 
#hist(x.test$bathy)
#hist(train.data.100m$bathy)

# Add the depthClass to the test data ... what are the levels? - they are [0,1,2,3]
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
# 1. DIVE DATA: Assemble the data across regions ... 
rm( 'x.sub', 'x.test')
x.test <- Assemble.Coast.Test( dive.data )
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