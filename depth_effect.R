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


#========================== IDS Tests by depth ==============================
# 2020/05/10: Below needs updating based on objectives.
# 2020/06/25: Revamped to do IDS by depth. IDs no longer pooled, and only NCC and HG 
#             regions considered. 


#=====================================================================
#-- Part 2: Evaluation of the regional models ... 
#   zRibbons are on [0, 5], representing [ITD, 50+]

# Not evaluating each region by each IDS for each ribbon anymore ... 
# code re-used elsewhere

xx <- NULL

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

#out.file <- 'Depth_compare_regions.csv'
#write.csv( results.table, file = file.path(output.dir, out.file) )
#results.table <- read.csv(file = file.path(output.dir, 'Depth_compare_regions.csv') )




#-- FIN.
