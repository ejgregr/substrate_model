#*******************************************************************************
# Script Name: 	  Plots
# Script Purpose:	Generate a series of plots summarizing Random Forest results
# Script Author: 	Cole Fields
# Script Date: 	  Jan. 28, 2019
# R Version:	    R version 3.5.2 (2018-12-20)
# R Packages:	    rgdal, ggplot2, raster, sp, magrittr, dplyr, caret, binr, tidyr, reshape2
#*******************************************************************************

# Adapted for substrate paper 
# Edward Gregr,
# 30 March 2020

# Much of this code originally written to address 1-step 2-step question. That's why 
# some of the proportion code scrapes rasters. 
# THe code has been adapted to look at proportions, statistics by class and depth, faceted across regions. 

# - For key metrics (Accuracy, TSS, TNR), bar plots of values by depth bin, faceted by region.
# - Also look at class prevalence bar plots for obs and predicted values, and class-based statistics. 
# NOTES:
#   - To ensure bars don't expand to fill gaps, need to fill the missing data with 'NA'.
#   - Consider whether bars are the best way to present these results.
#   - Contingency table as scatter plot - cool but hard to interpret without marginal values. 
#       Table itself seems better. 
#*******************************************************************************

#-------------------------------------------------------------------------------
# FUNCTIONS
# Revisions (EJG):
# Model building results ... 
# - data frame construction removed - they are built in the summary script
#   Build.Plots() : Build the heat maps for the build results.
#   Plot.Obs.Pred.Prevalence.Build() 
#   Plot.Stat.By.Depth.For.Regions() : User-specified statistic is plotted
#   Plot.Class.Stats.For.Regions() : User and Producer accuracies
# For the IDE ... 
#   Plot.Stat.By.IDS.For.Regions() : Plots selected integrated metrics by IDS and region
#-------------------------------------------------------------------------------

#===================== PLAYING WITH PALETTES =============================
# 2020/06/19: Standardized on some CB-friendly pallettes
library( RColorBrewer )
# display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE,
#                    colorblindFriendly=TRUE)


# A CB palette with 8 equally spaced colours. 
pal.cb8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#show.pal( pal.cb8 )

pal.RMSM <- pal.cb8[ c(7, 5, 1, 3) ]
#show.pal( pal.RMSM )

pal.cb2 <- pal.cb8[ c(5, 6) ]
#show.pal( pal.cb2 )

pal.cb3a <- pal.cb8[ c(3, 4, 7) ]
#show.pal( pal.cb2 )

pal.cb3b <- pal.cb8[ c(2, 5, 6) ]
#show.pal( pal.cb3 )

pal.cb4 <- brewer.pal(4, 'Dark2' )
#show.pal( pal.cb4 )

pal.heat.10 <- brewer.pal(10, "RdYlBu")
pal.heat.11 <- brewer.pal(11, "RdYlBu")

# show.pal( pal.cb8 )
# show.pal( pal.cb3b )

#-- Old palettes (Coles and PNW) dropped in favour of CB above.
# pal.cole <- c("#c6c6c6", "#c2f3f9", "#7fc9d8")

#-- devtools required above to properly install colours ... the first time?
#devtools::install_github("jakelawlor/PNWColors") 
#library(PNWColors)

#--- Could play with colours using the following code ... but No. 
# Use luminance=45 instead of default 65?
# Reduced luminance makes colours darker. Wld like that I think.

# To change luminance, convert from Hex RGB -->  HSL (Hue/saturation/luminance), 
# change Lum, and convert back. 
# Can also just change in ggplot() call:

#   ggplot(df, aes(x=cond, y=yval, fill=cond)) + geom_bar(stat="identity") +
#     scale_fill_hue(l=40)

#--------------------------------------------
# specify h as whole input degrees (e.g 0-360)
# s = 0.0 - 1 (0 - 100%)
# l = 0.0 - 1, (0 - 100%)
# returns output from R's rgb() functin
hsl_to_rgb <- function(h, s, l) {
  h <- h / 360
  r <- g <- b <- 0.0
  if (s == 0) {
    r <- g <- b <- l
  } else {
    hue_to_rgb <- function(p, q, t) {
      if (t < 0) { t <- t + 1.0 }
      if (t > 1) { t <- t - 1.0 }
      if (t < 1/6) { return(p + (q - p) * 6.0 * t) }
      if (t < 1/2) { return(q) }
      if (t < 2/3) { return(p + ((q - p) * ((2/3) - t) * 6)) }
      return(p)
    }
    q <- ifelse(l < 0.5, l * (1.0 + s), l + s - (l*s))
    p <- 2.0 * l - q
    r <- hue_to_rgb(p, q, h + 1/3)
    g <- hue_to_rgb(p, q, h)
    b <- hue_to_rgb(p, q, h - 1/3)
  }
  return(rgb(r,g,b))
}

#-----------------------------
# r, g, b = 0.0 - 1 (0 - 100%)
# returns h/s/l in a vector, h = 0-360 deg, s = 0.0 - 1 (0-100%), l = 0.0 - 1 (0-100%)
rgb_to_hsl <- function(r, g, b) {
  val_max <- max(c(r, g, b))
  val_min <- min(c(r, g, b))
  h <- s <- l <- (val_max + val_min) / 2
  if (val_max == val_min){
    h <- s <- 0
  } else {
    d <- val_max - val_min
    s <- ifelse(l > 0.5, d / (2 - val_max - val_min), d / (val_max + val_min))
    if (val_max == r) { h <- (g - b) / d + (ifelse(g < b, 6, 0)) }
    if (val_max == g) { h <- (b - r) / d/ + 2 }
    if (val_max == b) { h <- (r - g) / d + 4 }
    h <- (h / 6) * 360
  }
  return(c(h=h, s=s, l=l))
}

#-- To make the above work, Need to convert EACH HEX colour descrip into 3 numbers. 
#   Then back again. Defer since this is getting increasingly tangential. 
#   pal.RMSM
#   strtoi( substr( pal.RMSM[1], 2, 3 ), 16L )
#   rbb_to_hsl( pal.RMSM ) 

#-------------------------------
# Display the specified palette.
show.pal <- function( aa ){
  n <- length(aa)
  image( 1:n, 1, as.matrix(1:n), col=aa, axes=FALSE , xlab="", ylab="" )
}


#========================== PLOTTING FUNCTIONS ==========================

#-------------------------------
#-- Heat maps of the build tables.

#-Given a df with a column named 'Stat', pull the rows with a specified stat and create the heatmap.
# The heat map columns correspond to the substrate classes. These are columns in the df.
# Also, the Region column has attributes describing which model it is. 
Heat.Build.Class.Stats <- function( df, wtd = T, theStat, pal, w = 800, h = 600, txtCol ) {

  #--- Data prep ... 
  # Extract rows of interest 
  statDF <- df[ df$Stat == theStat, ] 
  
  # Remove the label colums, and name the rows after the regions 
  foo <- statDF[, -c(1:2) ]
  row.names( foo ) <- statDF$Region
  
  #The values to display on each cell 
  numtab <- round( data.matrix(foo), 2)
  
  today <- substr( Sys.time(), 1, 10)
  wtOrNot <- if(wtd==T) 'wtd' else 'noWts'
  fname <- paste0( 'heat_', theStat, '_Build_', wtOrNot, '_', today, '.png' )
  
  png( file.path( results.dir, fname), width = w, height = h )
       
  superheat( foo, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
             X.text = numtab, X.text.col = txtCol, X.text.size = 10,
             order.rows = rev(1:nrow(foo)),
             left.label.col = 'white', left.label.text.size = 10,
             bottom.label.col = 'White', bottom.label.size = 0.1, bottom.label.text.size = 10 )
  
  dev.off()
}

#---------------------------------------------------
#-- Heat map of predictor importance for all models. 
Heat.Build.Var.Import <- function( df, pal, w = 800, h = 600, txtCol ) {
  
  #-- create a rank table from the source df: 
  y <- data.frame(
    rbind( rank( -as.vector( df["Coast",] )),
           rank( -as.vector( df["HG",]    )),
           rank( -as.vector( df["NCC",]   )),
           rank( -as.vector( df["WCVI",]  )),
           rank( -as.vector( df["QCS",]   )),
           rank( -as.vector( df["SOG",]   ))
    ))
  row.names( y ) <- row.names(df) 
  
  #-- Define the text for the tigure using proportion of maximum ... 
  tab <- as.matrix(df)
  den <- apply( tab, 1, max)
  t2  <- round( apply( tab, 2, "/", den ), 2)
  
  # Zero out the NA in Coastal, text and value ... 
  t2[1,11] <- 'NA'
  y[1,11] <- 0

  #-- order according to the results ... someting funny about vectors require these 2 steps. 
  olist <- 1:ncol(y)
  olist  <- olist[ c(1,11,2,10,3,9,6,4,5,7,8) ]

  png( file.path( results.dir, 'heat_VariableImportance_Build.png' ), width = w, height = h )
  superheat( y, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
             X.text = t2, X.text.col = txtCol, X.text.size = 8, heat.lim = c(1,11),
             order.rows = rev(1:nrow(y)), order.cols = olist,
             left.label.col = 'white', left.label.text.size = 10,
             bottom.label.col = 'White', bottom.label.text.angle = 75, 
             bottom.label.size = 0.4, bottom.label.text.size = 10, bottom.label.text.alignment = 'right' )
  dev.off()
}


#-------------- FACETED Build Statistics --------------------

#-- Class prevalence of obs and pred during build, faceted by region.
# Source: csv built by build.summary(). 
# NOTE: the order of the Obs/Pred rows is irrelevant because melted here. 
Plot.Obs.Pred.Prevalence.Build <- function( dat.table, apal, sz=15, lx=0, ly=0 ){
  
  foo <- melt( dat.table, id.var = c('Region', 'Source'))
  
  a <- foo %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("Coast", "HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    
    ggplot(aes(x = variable, y = value, fill = Source)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = 'Prevalence' ) +
    facet_grid(. ~ Region) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(  text = element_text(size=sz ),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
            
            # legend stuff          
            legend.position = c(lx, ly),
            legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
    )
  
  return(a)
}

#2020/09/04: New plot to examine class-based stats across 3 RF models and all regions
# Needs to be applied oncve per statistic (accuracy, specificity, and reliability)
Plot.ClassStats.IDE <- function( dat.table, metr, apal, sz=15, lx=0, ly=0 ){
  
  foo <- melt( dat.table, id.var = c('Region', 'Model'))
  
  a <- foo %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("Coast", "HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    
    ggplot(aes(x = variable, y = value, fill = Model)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = metr ) +
    facet_grid(. ~ Region) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(  text = element_text(size=sz ),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
            
            # legend stuff          
            legend.position = c(lx, ly),
            legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
    )
  
  return(a)
}

#-- Class prevalence of obs vs. prediction across study area, faceted by region.
Plot.Obs.RegionPred.Prevalence <- function( dat.table, apal ){
  
  # change numbers to proportions
  y <- data.frame( dat.table[,c(3:6)]/rowSums( dat.table[,c(3:6)] ))
  # add results back to the ID columns
  y <- cbind( dat.table[,c(1,2)], y)
  names(y)[1] <- 'Source'
  
  foo <- melt( y, id.var = c('Stat', 'Region'))
  
  a <- foo %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("Coast", "HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    
    ggplot(aes(x = variable, y = value, fill = Source)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = 'Prevalence' ) +
    facet_grid(. ~ Region) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
  
  return(a)
}

#--- Producer and User stats (Build) by Class for Regions ---
Plot.Class.Stats.For.Regions <- function( dat.table, apal){
  
  foo <- melt( dat.table, id.var = c('Region', 'Stat') )
  
  a <- foo  %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("Coast", "HG", "NCC", "WCVI", "QCS", "SOG"))) %>%

        ggplot(aes(x = variable, y = value, fill = Stat)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = 'Score' ) +
    facet_grid(. ~ Region ) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
  
  return(a)
}

#--- TSS by Depth class within Regions --- includes 20 and 100 m models.
Plot.TSS.By.Depth.For.Regions <- function( dat.table, apal){

  a <- dat.table  %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    
    ggplot(aes(x = Ribbon, y = TSS, fill = Model)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = 'TSS' ) +
    facet_grid(. ~ Region ) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
  
  return(a)
}

#--- Non-random stat by Depth class withing Regions 
# Plots difference between the 20 and 100 m models for specified stat.
Plot.Stat.By.Depth.For.Regions <- function( dat.table, stat, apal){
  
  ylab <- paste0( stat, ' (100 m - 20 m )' )
  a <- dat.table  %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    
#    ggplot(aes(x = Ribbon, y = diffTSS, fill = Model)) +
    ggplot(aes(x = Ribbon, y = diff )) +
#    geom_bar(stat = "identity", width = .6, position = "dodge") +
    geom_bar(stat = "identity", width = .8 ) +
    labs( x = NULL, y = ylab) +
    facet_grid(. ~ Region ) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(  text = element_text(size=15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4), 
            
            # formatting the facet strip
            strip.text = element_text(size=rel(1.2))
        )
    return(a)
}

# Params: sz used to set text size; lx, ly set position of legend
Plot.Stats.By.IDS.For.Regions <- function( df, apal, sz=20, lx=0, ly=0 ){
  
  #  x <- df[, c( 'Region', 'IDS', 'TSS', 'Accuracy', 'TNRWtd' )]
  x <- df[, c( 'Region', 'IDS', 'TSS', 'Accuracy', 'TNR' )]
  
  # Adjust for different random baselines
  x$TSS <- x$TSS - 0.5
  x$Accuracy <- x$Accuracy - 0.25
  x$TNR <- x$TNR - 0.75
  
  #  names(x) <- c( 'Region', 'IDS', 'TSS', 'Accuracy', 'Specificity' )
  
  foo <- melt( x, id.vars = c('Region','IDS') )
  
  a <- foo %>%
    ggplot(aes(x = IDS, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = 'Score' ) +
    facet_grid(. ~ Region) +
    
    theme_bw() +
    theme(text         = element_text( size=sz ), 
          axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4),
          
          # legend stuff          
          legend.position = c(lx, ly),
          legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
          
    ) +
    #    scale_fill_manual(name = 'Independent\ntest data',
    scale_fill_manual(name = 'Metric',
                      values = apal)
  
  return( a )
}

Plot.Pontius.By.Depth.For.Regions <- function( df, apal, sz=20 ){
  
  y <- df[ , c('Region', 'Ribbon', 'Accuracy','Shift', 'Exchange', 'Quantity' )]
  y <- melt(y, id.vars = c('Ribbon','Region') )
  
  a <- ggplot(y, aes(x=Ribbon, y=value, fill=forcats::fct_rev(variable), order=desc(variable) )) +
  geom_bar(stat="identity") +
  labs( y = NULL, legend = 'Metric' ) +
  facet_grid(. ~ Region) +
    
  theme( text = element_text( size=sz ),
         axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4)
  ) +
#  scale_y_continuous( limits=c(0.5, 1.0) ) +
#  ylim( 0.5, 1.0 ) +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  scale_fill_manual( values = apal,
                     name = 'Metrics')
  
  return( a )
}


#-------------- FACETED IDE Statistics --------------------

#--------------------------------------------------------
#-- Simple facet of the sample size by region for context.
Plot.Obs.By.IDS.For.Regions <- function( df, apal, sz = 20, lx=0, ly=0 ){
  
  foo <- melt( df, id.vars = c('Region','IDS') )
  
  a <- foo %>%
    ggplot(aes(x = variable, y = value, fill = IDS)) +

        geom_bar(stat = "identity", width = .8, position = "dodge") +
    #geom_bar(stat = "identity", width = .8, position = position_dodge2(preserve = "single", padding = 0) ) +
    labs( x = NULL, y = 'N' ) +
    facet_grid(. ~ Region) +
    theme_bw() +
    theme( text         = element_text(size = sz), 
           axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4),   # vjust needed to centre rotated names on tick

           # formatting the facet strip
           strip.text = element_text(size=rel(1.2)),

           # legend position          
           legend.position = c(lx, ly), 
#           legend.box.margin = c(50, 50, 50, 50)
           legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
           ) + 
    
    scale_fill_manual(name = 'Independent\ntest data',
                      values = apal)
  return( a )
}




#------------------------------------------------
#-- Producer and User stats by Class for Regions
#   NEEDS UPDATED DATA TABLE ... 
IDS.Class.Stats.For.Regions  <- function( df, apal, sz=20, lx=0, ly=0 ){

  foo <- melt( df, id.vars = c('Region','IDS', 'Stat') )
  foo <- foo[ foo$Region == 'HG', ]
  
  a <- foo %>%
    ggplot(aes(x = variable, y = value, fill = Stat)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = 'Score' ) +
    facet_grid(. ~ IDS) +
    
    theme_bw() +
    theme(text         = element_text( size=sz ), 
          axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4),
          
          # legend stuff          
          legend.position = c(lx, ly),
          legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
          
    ) +
    #    scale_fill_manual(name = 'Independent\ntest data',
    scale_fill_manual(name = 'Metric',
                      values = apal)
  
  return( a )
}

#-------------------------------
#-- Map prevalence compared to model test prevalence
#rm('maprev', 'bs', 'a', 'b', 'c', 'y')
Plot.Pred.Map.Prevalence <- function( maprev, bs, apal, sz=20){
# TAKES: maprev: Map prevalence saved as part of study area predictions
#        bs: Summary of the build prevalences
  
  a <- rbind( maprev$Coast2, maprev$HG2, maprev$NCC2, maprev$WCVI2, maprev$QCS2, maprev$SOG2 )
  colnames( a ) <- c('Hard','Mixed','Sand','Mud')
  a <- a/100
  b <- factor( c('Coast','HG','NCC','WCVI','QCS','SOG'), 
               levels = c('Coast','HG','NCC','WCVI','QCS','SOG') )
  
  a <- cbind( Region = b, data.frame(a) )
  a <- cbind( Stat = 'Map', a )
  
  # Pull prediction results from build results
  b <- bs$build.results.ClassPrev[ bs$build.results.ClassPrev$Stat == 'Pred',]
  b <- cbind( b[ , c(1,2)], 
              b[ , c(-1,-2)] / rowSums(b[ , c(-1,-2)]) )

  c <- rbind(a, b)
  levels(c$Stat) <- c('Map','Pred', 'Points')
  c$Stat[ c$Stat == 'Pred'] <- 'Points'
  y <- melt(c, id.vars = c('Region','Stat') )
  y <- mutate(y, Scale = factor(Stat, levels=c("Points", "Map")))
  
  ggplot(y, aes(x=Scale, y=value, fill=variable, order=desc(variable) )) +
    
    geom_bar(stat="identity") +
    facet_grid(. ~ Region) +
    scale_fill_manual( values = apal,
                       name = 'Class') +
    theme_bw() +
    labs( y = NULL, legend = 'Prevalence' ) +
    theme( text = element_text( size=sz ),
           axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4)
    )
}

#-----------------------------------------------------------------
Plot.Pontius.By.IDS.For.Regions <- function( df, apal, sz=20 ){
  
  y <- melt(df, id.vars = c('IDS', 'Region') )
  
  a <- y %>%
  ggplot(aes(x=IDS, y=value, fill=forcats::fct_rev(variable), order=desc(variable) )) +
    geom_bar(stat="identity") +
    labs( y = NULL, legend = 'Metric' ) +
        facet_grid(. ~ Region) +
    
    theme( text = element_text( size=sz ),
           axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4)
           ) +
    
    scale_fill_manual( values = apal,
                       name = 'Metrics')
  
  return( a )
}

#------------------------------------------------------------------------
Plot.Pontius.By.IDS.Depth.For.Regions <- function( df, apal, sz=20 ){
  
  y <- melt(df, id.vars = c('Ribbon','IDS') )
  
  a <- y %>%
    ggplot(aes(x=Ribbon, y=value, fill=forcats::fct_rev(variable), order=desc(variable) )) +
    geom_bar(stat="identity") +
    labs( y = NULL, legend = 'Metric' ) +
    facet_grid(. ~ IDS) +
    
    theme( text = element_text( size=sz ),
           axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4)
    ) +
    
    scale_fill_manual( values = apal,
                       name = 'Metrics')
  
  return( a )
}

#--- Integrated Statistics by IDS faceted by Region V.2 ---
#   Compare TSS (diff from random) vs. Quantity and Allocation (is this diff from accuracy)
Plot.5Stats.By.IDS.For.Regions <- function( df, apal, sz=20, lx=0, ly=0 ){
  
  x <- df[, c( 'Region', 'IDS', 'TSS', 'Accuracy', 'TNRWtd', 'Quantity', 'Exchange' )]
  
  foo <- melt( x, id.vars = c('Region','IDS') )
  
  a <- foo %>%
    ggplot(aes(x = variable, y = value, fill = IDS)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = 'Score' ) +
    facet_grid(. ~ Region) +
    
    theme_bw() +
    theme(text         = element_text( size=sz ), 
          axis.text.x  = element_text(angle = 90, hjust = 1),
          
          # legend stuff          
          legend.position = c(lx, ly),
          legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
          
    ) +
    #    scale_fill_manual(name = 'Independent\ntest data',
    scale_fill_manual(name = 'Metric',
                      values = apal)
  
  return( a )
}


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
         fill = pal.RMSM, title= NA, bg = NA, box.col = NA)
  
  dev.off()
}


Error.Matrix <- function( rfModel, tpts ) {
  w <- predict(rfModel, tpts)
  a <- caret::confusionMatrix( w$predictions, tpts$BType4 )$table
  
  b <- sum(a)
  c <- cbind( a, 'Prev' = rowSums(a) / b )
  c <- cbind( c, 'User' = diag(a)/rowSums(a) )
  
  c <- rbind( c, 'Prev' = c( colSums(a) / b, 1, 1) ) 
  d <- rbind( c, 'Prod' = c( diag(a)/colSums(a), 1, 1) ) 
  
  return( d )
}


