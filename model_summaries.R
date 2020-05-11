#*******************************************************************************
# Script:  model_summaries.R
# Created: 01 April 2020. EJG
#
# Purpose: Using the loaded observational data and the 5 RF models built, 
#   prepare a bunch of data summaries to then visualize/ plot.
# Methods: Hand-off to plots is via CSV files containing summary tables. 
#
# NOTES: 2020/05/10: Out-of-date code for writing to a text file deleted.
#
#*******************************************************************************







#-- Requires data structures and models to have been built ... 
# Different from Build predictions cuz uses full study area rasters. 
# DEFERRED - Hoping Cole can scrape the rasters. 
Summarize.Predictions <- function() {

  out <- NULL
  # coastwide ... 
  
    tdat <-  obs.100m.data[ obs.100m.data$TestDat == 1, ] #test
  x <- predict( rf.region.Coast, tdat )
  y <- table(x$predictions)
  out <- rbind( out, y)
              
  # HG: Load rasters, drop NAs, predict for region   ... 
  tdat <- Load.Predictors( paste0( predictor.dir, "/", 'HG' ))
  
  z <- tdat$bathy
  dim( z )
  str( z )
  z <- z[ !is.na( z  ), ]
  
  
  head(z)
  
  tdat <- tdat[ !is.na( tdat$bathy ), ]
  
  
  
  names(tdat)
  summary(tdat)
  
  
  
  x <- predict( rf.region.HG, tdat@data )
  y <- table(x$predictions)
  out <- rbind( out, y)
  


} 


