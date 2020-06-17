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

#--- Build some palettes ... Try brightening Cole's colours a bit 
pal.cole <- c("#c6c6c6", "#c2f3f9", "#7fc9d8")
pal.3.win <- pnw_palette( 'Winter', 10 )[ c(9,6,3) ]
pal.3.bay <-pnw_palette( 'Bay', 6 )[ c(6,4,1) ]

#-- Want to re-order colours cuz colour ramp not required ... 
pal.6.shuk <- pnw_palette( 'Shuksan', 6 )
pal.6.bay  <- pnw_palette( 'Bay', 6 )
pal.5 <- pnw_palette( 'Bay', 5 )
pal.4 <- pnw_palette( 'Bay', 4 )
pal.10 <- pnw_palette( 'Bay', 10 )
pal.11 <- rev(pnw_palette( 'Bay', 11 ))

pal.map <- c( pal.10[10], pal.10[4], pal.3.win[1], pal.6.shuk[1] )



show.pal <- function( aa ){
  n <- length(aa)
  image( 1:n, 1, as.matrix(1:n), col=aa, axes=FALSE , xlab="", ylab="" )
}

#show.pal( pal.10 )
#show.pal( pal.map )




# Heat maps of the build tables
Plot.Build.Class.Stats <- function( df, pal, w = 800, h = 600 ) {
  
  # 1) User accuracy:
  usr <- df[ df$Stat == 'User', ] 
  foo <- usr[, -c(1:2) ]
  row.names( foo ) <- usr$Region
  
  #The cell values to display ... 
  tab <- round( data.matrix(foo), 2)
  
  png( file.path( results.dir, 'heat_UserAccuracy_Build.png'), width = w, height = h )
  superheat( foo, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
             X.text = tab, X.text.col = 'lightblue', X.text.size = 10,
             order.rows = rev(1:nrow(foo)),
             left.label.col = 'white', left.label.text.size = 10,
             bottom.label.col = 'White', bottom.label.size = 0.1, bottom.label.text.size = 10 )
  dev.off()
  
  # 2) Producer accuracy  ... 
  prd <- df[ df$Stat == 'Prod', ] 
  foo <- prd[, -c(1:2) ]
  row.names( foo ) <- prd$Region
  
  #The cell values to display ... 
  tab <- round( as.matrix(foo), 2)
  
  png( file.path( results.dir, 'heat_ProducerAccuracy_Build.png'), width = w, height = h )
  superheat( foo, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
             X.text = tab, X.text.col = 'lightblue', X.text.size = 10,
             order.rows = rev(1:nrow(foo)),
             left.label.col = 'white', left.label.text.size = 10,
             bottom.label.col = 'White', bottom.label.size = 0.1, bottom.label.text.size = 10 )
  dev.off()
}


#-- Heat map of predictor importance for all models ... 
Plot.Build.Var.Import <- function( df, pal, w = 800, h = 600 ) {

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
  
  #-- order according to the results ... 
  olist <- c(1,11,2,10,3,9,6,4,5,7,8)

  png( file.path( results.dir, 'heat_VariableImportance_Build.png' ), width = w, height = h )
  superheat( y, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
             X.text = t2, X.text.col = 'lightblue', X.text.size = 8, heat.lim = c(1,11),
             order.rows = rev(1:nrow(y)), order.cols = olist,
             left.label.col = 'white', left.label.text.size = 10,
             bottom.label.col = 'White', bottom.label.text.angle = 75, 
             bottom.label.size = 0.4, bottom.label.text.size = 10, bottom.label.text.alignment = 'right' )
  dev.off()
}


#-- Heat map of model performance by region and depth class. 
Plot.Region.Class.Stats <- function( csv, res, pal, w = 800, h = 600 ) {
  
  # grab the Variable imporatance  table ... 
  # csv <- 'Build_results_varImportance.csv'
  a <- read.csv( file = file.path( results.dir, csv ))
  a <- a[ a$Model == res, -c(1,2)]

  b <- rbind( a[ a$Region =="HG", "TSS"],
         a[ a$Region =="NCC", "TSS"],
         a[ a$Region =="WCVI", "TSS"],
         a[ a$Region =="QCS", "TSS"],
         a[ a$Region =="SOG", "TSS"]
    )
  rownames( b ) <- bioregions
  colnames( b ) <- c( 'ITD','0-5','5-10','10-20','20-50','50+' )
  
  #Reverse the rows so it plots in correct order.  
  b <- b[seq(dim(b)[1],1),]

  png( file.path( results.dir, paste0( 'heat_TSS_byClassForResgions_', res, '.png' )), width = w, height = h )
  superheat( b, heat.pal = pal, legend = F, grid.hline = F, grid.vline = F, scale = F,
             X.text = round(b,2), X.text.col = 'grey50', X.text.size = 8,
             pretty.order.cols = F, pretty.order.rows = F,
             left.label.col = 'white', left.label.text.size = 10,
             bottom.label.col = 'White', bottom.label.text.angle = 45, 
             bottom.label.size = 0.4, bottom.label.text.size = 10, bottom.label.text.alignment = 'right' )
  dev.off()
}


#-------------- FACETED Build Statistics --------------------

#-- Class prevalence of obs and pred during build, faceted by region.
# Source: csv built by build.summary(). 
# NOTE: the order of the Obs/Pred rows is irrelevant because melted here. 
Plot.Obs.Pred.Prevalence.Build <- function( dat.table, apal ){
  
  # change numbers to proportions
  y <- data.frame( dat.table[,c(3:6)]/rowSums( dat.table[,c(3:6)] ))
  # add results back to the ID columns
  y <- cbind( dat.table[,c(1,2)], y)
  names(y)[1] <- 'Source'
  
  foo <- melt( y, id.var = c('Region', 'Source'))
  
  a <- foo %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("Coast", "HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    
    ggplot(aes(x = variable, y = value, fill = Source)) +
    geom_bar(stat = "identity", width = .6, position = "dodge") +
    labs( x = NULL, y = 'Prevalence' ) +
    facet_grid(. ~ Region) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
  
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
    geom_bar(stat = "identity", width = .6, position = "dodge") +
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
    geom_bar(stat = "identity", width = .6, position = "dodge") +
    labs( x = NULL, y = 'Score' ) +
    facet_grid(. ~ Region ) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
  
  return(a)
}


#--- TSS by Depth class withing Regions --- includes 20 and 100 m models.
Plot.TSS.By.Depth.For.Regions <- function( dat.table, apal){

    #foo <- melt( dat.table )
  
  a <- dat.table  %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    
    ggplot(aes(x = Ribbon, y = TSS, fill = Model)) +
    geom_bar(stat = "identity", width = .6, position = "dodge") +
    labs( x = NULL, y = 'TSS' ) +
    facet_grid(. ~ Region ) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(text = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
  
  return(a)
}


#-------------- FACETED IDE Statistics --------------------

#--- Simple facet of the sample size by region for context. ---
Plot.Obs.By.IDS.For.Regions <- function( df, apal, sz = 20, lx=0, ly=0 ){
  
  foo <- melt( df, id.vars = c('Region','IDS') )
  
  a <- foo %>%
    ggplot(aes(x = variable, y = value, fill = IDS)) +
    # dodge2 allows the bar size to be preserved for missing values, but they still get recentred if one is missing
    # that's why resorted to padding the input data. Retained for posterity. position=dodge wld currently be enough.
    geom_bar(stat = "identity", width = .8, position = position_dodge2(preserve = "single", padding = 0) ) +
    labs( x = NULL, y = 'N' ) +
    facet_grid(. ~ Region) +
    
    theme_bw() +
    theme( text         = element_text(size = sz), 
           axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.4),   # vjust needed to centre rotated names on tick

           # formatting the facet strip
           strip.text = element_text(size=rel(1.2)),
           strip.background = element_rect(fill="lightblue", colour="black", size=1),
  
           # legend position          
           legend.position = c(lx, ly), 
#           legend.box.margin = c(50, 50, 50, 50)
           legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
           ) + 
    
    scale_fill_manual(name = 'Independent\ntest data',
                      values = apal)
  return( a )
}


#--- Integrated Statistics by IDS faceted by Region ---
#   Statistics are hard-coded as list of 3.
#   Params: sz used to set text size; lx, ly set position of legent
Plot.Stats.By.IDS.For.Regions <- function( df, apal, sz=20, lx=0, ly=0 ){
  
  x <- df[, c( 'Region', 'IDS', 'TSS', 'Accuracy', 'TNRWtd' )]
  
  foo <- melt( x, id.vars = c('Region','IDS') )
  
  a <- foo %>%
    ggplot(aes(x = IDS, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = .6, position = "dodge") +
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


#--- Producer and User stats by Class for Regions
#   NEEDS UPDATED DATA TABLE ... 
IDS.Class.Stats.For.Regions  <- function( df, apal, sz=20, lx=0, ly=0 ){

  foo <- melt( df, id.vars = c('Region','IDS', 'Stat') )
  foo <- foo[ foo$Region == 'HG', ]
  
  a <- foo %>%
    ggplot(aes(x = variable, y = value, fill = Stat)) +
    geom_bar(stat = "identity", width = .6, position = "dodge") +
    labs( x = NULL, y = 'Score' ) +
    facet_grid(. ~ IDS) +
    
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


  # a <- df %>%
  #   ggplot(aes(x = Test.Data, y = Accuracy, fill = Test.Data)) +
  #   geom_bar(stat = "identity", width = .6, position = "dodge") +
  #   facet_grid(. ~ Model) +
  #   scale_fill_manual(values = my.colours) +
  #   theme_bw() +
  #   theme(text = element_text( size=sz )) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))



#-- Map prevalence compared to model test prevalence
Plot.Pred.Map.Prevalence <- function( maprev, bs){
# TAKES: maprev: Map prevalence saved as part of study area predictions
#        bs: Summary of the build prevalences
  
  a <- rbind( maprev$Coast2, maprev$HG2, maprev$NCC2, maprev$WCVI2, maprev$QCS2, maprev$SOG2 )
  colnames( a ) <- c('Hard','Mixed','Sand','Mud')
  a <- round( a/100, 2)
  a <- cbind( Region = c('Coast','HG','NCC','WCVI','QCS','SOG'), data.frame(a) )
  a <- cbind( Stat = 'Map', a )
  
  b <- bs$build.results.ClassPrev[ bs$build.results.ClassPrev$Stat == 'Pred',]
  b <- cbind( b[ , c(1,2)],
              round( b[ , c(-1,-2)] / rowSums(b[ , c(-1,-2)]), 2) )
  
  c <- rbind(a, b)
  
  y <- melt(c, id.vars = c('Region','Stat') )
  y <- mutate(y, Scale = factor(Stat, levels=c("Pred", "Map")))
  
  ggplot(y, aes(x=Scale, y=value, fill=variable, order=desc(variable) )) +
    geom_bar(stat="identity") +
    facet_grid(. ~ Region) +
    scale_fill_manual( values = pal.map,
                       name = 'Class') +
    labs( y = NULL, legend = 'Prevalence' ) +
    theme( text = element_text( size=20 ),
    )
}



Plot.Pontius.By.IDS.For.Regions <- function( df, apal, sz=20 ){
  
  y <- melt(df, id.vars = c('Region','IDS') )
  
  a <- y %>%
  ggplot(aes(x=IDS, y=value, fill=forcats::fct_rev(variable), order=desc(variable) )) +
    geom_bar(stat="identity") +
    labs( y = NULL, legend = 'Metric' ) +
        facet_grid(. ~ Region) +
    
    theme( text = element_text( size=sz ),
           axis.text.x  = element_text(angle = 90, hjust = 1)
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
    geom_bar(stat = "identity", width = .6, position = "dodge") +
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


#-- Integrated Statistics by class, for each IDS by region.
#-- Nice idea, but lots  of small sample sizes here. Are things really comparable this decomposed?
#   DEFERRED for NOW.

    # # Build the data ...
    # a.stat <- 'Quantity'
    # dat.table <- read.csv(file = file.path(output.dir, 'Depth_compare_regions.csv') )
    # x <- dat.table[[ 'Quantity' ]]
    # x <- cbind( dat.table[, c( 'Model', 'Test.Data', 'Ribbon', 'N')], 'Statistic' = x)
    # 
    # #-- Structure the data a bit ... 
    # levels( results.table$Ribbon ) <- c( 'ITD','0-5','5-10','10-20','20-50','50+' )
    # 
    # # Plot and save ...
    # a <- Plot.Stat.By.Depth.For.Regions( x, a.stat, pal.pnw )
    # # save plot to wd
    # ggsave( paste0( 'Stat ', a.stat, ' by IDE by Depth for Regions.png'), a, dpi = 300, width = 16, height = 10, path = output.dir)
    # 
    # #---Statistic by depth faceted by Region ---
    # # Table is custom built for the specified Statistic ... 
    # Plot.Stat.By.Depth.For.Regions <- function( a.table, theStat, apal ){
    #   
    #   a <- a.table %>%
    #     ggplot(aes(x = Ribbon, y = Statistic, fill = Test.Data)) +
    #     geom_bar(stat = "identity", width = .6, position = "dodge") +
    #     ylab( theStat ) +
    #     facet_grid(. ~ Model) +
    #     
    #     theme_bw() +
    #     theme(text         = element_text(size=15), 
    #           axis.text.x  = element_text(angle = 90, hjust = 1),
    #     ) +
    #     scale_fill_manual(name = 'Independent\ntest data',
    #                       values = apal)
    #   
    #   return( a )
    # }



# facetted plots for proportions
facet.props <- function(master.df, working.dir){
  # hex colours for plots
  my.colours <- c("#c6c6c6", "#c2f3f9", "#7fc9d8")

  # plot proportions of observations and predicted substrate by class
  sub.bars <- master.df %>%
    ggplot(aes(x = Substrate, y = Percent, fill = Source)) +
    geom_bar(stat = "identity", width = .6, position = "dodge") + 
    scale_fill_manual(values = my.colours) +
    ylab("Percent") + 
    ylim(0, 100) + 
    ggtitle("Proportions of Modelled Substrates (One-Step and Two-Step) and Observations") +
    theme_bw() +
    theme(text = element_text(size=15)) +
    facet_grid(. ~ Region) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # save plot to wd
  ggsave("Proportions of Modelled Substrates (One-Step and Two-Step) and Observations.png", sub.bars, dpi = 300, width = 16, height = 10, path = working.dir)
  
  return(sub.bars)
}



#-------------------------------------------------------------------
#------ not yet implemented ------

# # plots from single dataframes of prop obs vs single-step and two-step
# plots <- function(dataframe, working.dir, region){
#   bars <- dataframe %>%
#     ggplot(aes(x = Substrate, y = Percent, fill = Source)) +
#     geom_bar(stat = "identity", width = .4, position = "dodge") + 
#     scale_fill_brewer(palette = "PuBuGn") +
#     ylab("Percent") + 
#     xlab(paste(region, "Substrate")) + 
#     ylim(0, 100) + 
#     ggtitle("Proportions of Substrate Observations and Modelled Substrate") +
#     theme_bw() +
#     theme(text = element_text(size=15))
#   
#   # save plot to wd
#   ggsave(paste(region, "Substrate.png"), bars, dpi = 300, width = 12, height = 10, path = working.dir)
#   # return plot
#   return(bars)
# }
# 
# # aggregates performance by model resolution and then creates plots for kappa and accuracy - GROUP DID NOT WANT
# overall.plots.means <- function(df, working.dir){
#   # get mean kappa and accuracy with standard deviations
#   df100 <- df[df$Resolution == "100", ]
#   df20 <- df[df$Resolution == "20", ]
#   means.100 <- aggregate(df100[, c("Accuracy", "Kappa")], list(df100$depth.bin, df100$Method), FUN = function(x) c(mean = mean(x), sd = sd(x)))
#   means.100$Accuracy.mean <- means.100$Accuracy[, 1]
#   means.100$Accuracy.sd   <- means.100$Accuracy[, 2]
#   means.100$Kappa.mean    <- means.100$Kappa[, 1]
#   means.100$Kappa.sd      <- means.100$Kappa[, 2]
#   means.100$Resolution <- "100"
#   means.20 <- aggregate(df20[, c("Accuracy", "Kappa")], list(df20$depth.bin, df20$Method), FUN = function(x) c(mean = mean(x), sd = sd(x)))
#   means.20$Accuracy.mean  <- means.20$Accuracy[, 1]
#   means.20$Accuracy.sd    <- means.20$Accuracy[, 2]
#   means.20$Kappa.mean     <- means.20$Kappa[, 1]
#   means.20$Kappa.sd       <- means.20$Kappa[, 2]
#   means.20$Resolution <- "20"
#   
#   # bind means together
#   means <- rbind(means.100, means.20)
#   means$Accuracy <- NULL
#   means$Kappa <- NULL
#   names(means)[names(means) == 'Group.1'] <- 'depth.bin'
#   names(means)[names(means) == 'Group.2'] <- 'Method'
#   
#   # add means to df
#   df <- merge(df, means, by = c("Resolution", "depth.bin", "Method"), sort = FALSE)
#   
#   # accuracy plot
#   accuracy.p <- df %>%
#     ggplot(aes(x = depth.bin, y = Accuracy.mean, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     geom_errorbar(aes(ymin = Accuracy.mean - Accuracy.sd, ymax = Accuracy.mean + Accuracy.sd), width=.2,
#                   position = position_dodge(.9)) +
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Accuracy") + 
#     ylim(0, 1) + 
#     ggtitle("Accuracy by Model Resolution") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_wrap(. ~ as.factor(Resolution), scales="free_x", nrow=1) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Accuracy by Model Resolution.png", accuracy.p, dpi = 300, width = 16, height = 10, path = working.dir)
#   
#   # kappa plot
#   kappa.p <- df %>%
#     ggplot(aes(x = depth.bin, y = Kappa.mean, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     geom_errorbar(aes(ymin = Kappa.mean - Kappa.sd, ymax = Kappa.mean + Kappa.sd), width=.2,
#                   position = position_dodge(.9)) +
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Kappa") + 
#     ylim(0, 1) + 
#     ggtitle("Kappa by Model Resolution") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_wrap(. ~ as.factor(Resolution), scales="free_x", nrow=1) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   
#   # save plot to wd
#   ggsave("Kappa by Model Resolution.png", kappa.p, dpi = 300, width = 16, height = 10, path = working.dir)
# 
#   # return accuracy and kappa plots
#   return(list(accuracy.p, kappa.p))
# }
# 
# # creates plots for kappa and accuracy by region and depth bin
# overall.plots <- function(df, working.dir){
#   # hex colours for plots
#   my.colours <- c("#c2f3f9", "#7fc9d8")
#   
#   # accuracy plot
#   accuracy.p <- df[df$Metric == "Accuracy", ] %>%
#     ggplot(aes(x = depth.bin, y = Value, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_manual(values = my.colours) +
#     ylab("Accuracy") + 
#     ylim(0, 1) + 
#     ggtitle("Accuracy by Region and Depth Bin") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_wrap(. ~ as.factor(Region), scales="free_x", nrow=1) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Accuracy by Region and Depth Bin.png", accuracy.p, dpi = 300, width = 16, height = 10, path = working.dir)
#   
#   # kappa plot
#   kappa.p <- df[df$Metric == "Kappa", ] %>%
#     ggplot(aes(x = depth.bin, y = Value, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_manual(values = my.colours) +
#     ylab("Kappa") + 
#     ylim(0, 1) + 
#     ggtitle("Kappa by Region and Depth Bin") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_wrap(. ~ as.factor(Region), scales="free_x", nrow=1) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Kappa by Region and Depth Bin.png", kappa.p, dpi = 300, width = 16, height = 10, path = working.dir)
#   
#   # return accuracy and kappa plots
#   return(list(accuracy.p, kappa.p))
# }
# 
# class.plots <- function(df, working.dir){
#   
#   # get 
#   class.df %>% group_by(Resolution, depth.bin, Method, Metric) %>% mutate(Mean.Value = mean(Value), Std.Dev = sd(Value))
#   
#   
#   # sensitivity plot
#   sensSpec.p <- df[df$Metric %in% c("Sensitivity", "Specificity"), ] %>%
#     ggplot(aes(x = depth.bin, y = Value, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Sensitivity") + 
#     ylim(0, 1) + 
#     ggtitle("Sensitivity by Region and Substrate") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_wrap(vars(Region, Class, Metric)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Sensitivity by Region and Substrate.png", sensitivity.p, dpi = 300, width = 12, height = 12, path = working.dir)
# 
#   # specificity plot
#   specificity.p <- df %>%
#     ggplot(aes(x = depth.bin, y = Specificity, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Specificity") + 
#     ylim(0, 1) + 
#     ggtitle("Specificity by Region and Substrate") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_grid(Region ~ Class) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Specificity by Region and Substrate.png", specificity.p, dpi = 300, width = 12, height = 12, path = working.dir)
#   
#   # positive pred value plot
#   ppv.p <- df %>%
#     ggplot(aes(x = depth.bin, y = PosPredValue, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Positive Pred Value") + 
#     ylim(0, 1) + 
#     ggtitle("Positive Pred Value by Region and Substrate") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_grid(Region ~ Class) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Positive Pred Value by Region and Substrate.png", ppv.p, dpi = 300, width = 12, height = 12, path = working.dir)
#   
#   # negative pred value plot
#   npv.p <- df %>%
#     ggplot(aes(x = depth.bin, y = NegPredVal, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Negative Pred Value") + 
#     ylim(0, 1) + 
#     ggtitle("Negative Pred Value by Region and Substrate") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_grid(Region ~ Class) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Negative Pred Value by Region and Substrate.png", npv.p, dpi = 300, width = 12, height = 12, path = working.dir)
#   
#   # prevalence plot
#   prev.p <- df %>%
#     ggplot(aes(x = depth.bin, y = Prevalence, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Prevalence") + 
#     ylim(0, 1) + 
#     ggtitle("Prevalence by Region and Substrate") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_grid(Region ~ Class) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Prevalence by Region and Substrate.png", prev.p, dpi = 300, width = 12, height = 12, path = working.dir)
#   
#   # prevalence plot
#   balAcc.p <- df %>%
#     ggplot(aes(x = depth.bin, y = BalancedAccuracy, fill = Method)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_brewer(palette = "Paired") +
#     ylab("Balanced Accuracy") + 
#     ylim(0, 1) + 
#     ggtitle("Balanced Accuracy by Region and Substrate") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_grid(Region ~ Class) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Balanced Accuracy by Region and Substrate.png", balAcc.p, dpi = 300, width = 12, height = 12, path = working.dir)
#   
#   # return accuracy and kappa plots
#   return(list(accuracy.p, kappa.p))
# }
# 
# # creates plots for kappa and accuracy by region and depth bin
# compare.plots <- function(df, working.dir){
#   # hex colours for plots
#   my.colours <- c("#c2f3f9", "#7fc9d8")
#   
#   # accuracy plot
#   accuracy.p <- df[df$Metric == "Accuracy", ] %>%
#     ggplot(aes(x = depth.bin, y = Value, fill = Resolution)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_manual(values = my.colours) +
#     ylab("Accuracy") + 
#     ylim(0, .75) + 
#     ggtitle("Accuracy by Region and Depth Bin") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_wrap(as.factor(Method) ~ as.factor(Region), scales="free_x", nrow=2) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Accuracy by Region and Depth Bin.png", accuracy.p, dpi = 300, width = 16, height = 10, path = working.dir)
#   
#   # kappa plot
#   kappa.p <- df[df$Metric == "Kappa", ] %>%
#     ggplot(aes(x = depth.bin, y = Value, fill = Resolution)) +
#     geom_bar(stat = "identity", width = .6, position = "dodge") + 
#     scale_fill_manual(values = my.colours) +
#     ylab("Kappa") + 
#     ylim(0, .75) + 
#     ggtitle("Accuracy by Region and Depth Bin") +
#     theme_bw() +
#     theme(text = element_text(size=15)) +
#     facet_wrap(as.factor(Method) ~ as.factor(Region), scales="free_x", nrow=1) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   
#   # save plot to wd
#   ggsave("Kappa by Region and Depth Bin.png", kappa.p, dpi = 300, width = 16, height = 10, path = working.dir)
#   
#   # return accuracy and kappa plots
#   return(list(accuracy.p, kappa.p))
# }
# 
# # table of obs w/source (input should be master.dataframe)
# data.summary <- function(obs.list, working.dir){
#   obs.list[[1]]$Region <- 'Coast-Wide'
#   obs.list[[2]]$Region <- 'North Central Coast'
#   obs.list[[3]]$Region <- 'Strait of Georgia'
#   obs.list[[4]]$Region <- 'Queen Charlotte Strait'
#   obs.list[[5]]$Region <- 'Haida Gwaii'
#   obs.list[[6]]$Region <- 'West Coast Vancouver Island'
#   
#   # list of columns to keep
#   keeps <- c('BT_Sorc', 'BType4', 'TestDat', 'Rock', 'Region') 
#   
#   # drop all columns except for the keeps
#   obs.list <- lapply(obs.list, FUN = function(obs.df){
#     obs.df[, (names(obs.df) %in% keeps)]
#   })
#   
#   obs.tables <- lapply(obs.list, FUN = function(obs.df){
#     addmargins(table(obs.df$BT_Sorc, obs.df$BType4))
#   })
#   
#   temp.tables <- lapply(obs.list, FUN = function(obs.df){
#     table(obs.df$BT_Sorc, obs.df$BType4)
#   })
#   
# 
#   names(obs.tables) <- c('Coast-Wide', 'North Central Coast', 'Strait of Georgia', 'Queen Charlotte Strait', 'Haida Gwaii', 'West Coast Vancouver Island')
#   names(temp.tables) <- c('Coast-Wide', 'North Central Coast', 'Strait of Georgia', 'Queen Charlotte Strait', 'Haida Gwaii', 'West Coast Vancouver Island')
#   
#   prop.tables <- lapply(temp.tables, FUN = function(obs.tbl){
#     round(100 * addmargins(prop.table(obs.tbl)), 2)
#   })
#   
#   return(list(obs.tables, prop.tables))
#   
# }
# 
# 
# # coast-wide performance metrics by region (20m substrate regions)
# regional.divisions <- function(one.step, two.step, reg.dir, reg.name){
#   # make raster of coast-wide
#   coast.r.one <- raster(one.step)
#   coast.r.two <- raster(two.step)
#   
#   # rename substrate rasters
#   names(coast.r.one) <- "pred_oneStep"
#   names(coast.r.two) <- "pred_twoStep"
#   
#   # make spatial polygons layer from shapefile
#   region <- readOGR(dsn = reg.dir,layer = reg.name)
#   
#   # mask the shelf model with regional polygon
#   mask.one <- mask(coast.r.one, region)
#   mask.two <- mask(coast.r.two, region)
#   
#   # make raster stack of substrate models
#   r.stack <- stack(mask.one, mask.two)
#   
#   # # subset test data
#   # test.sp <- obs.sp[obs.sp$TestDat == 1, ]
#   # 
#   # # attach predicted raster values to test points 
#   # preds <- raster::extract(r.stack, test.sp, df = TRUE)
#   # 
#   # # convert to data frame; combine observations w/preds
#   # test.sp <- cbind(as.data.frame(test.sp), preds)
#   # 
#   # # add depth.bin column and populate with NAs
#   # test.sp$depth.bin <- NA
#   # 
#   # # sort df based on depth value
#   # test.sp <- test.sp[order(test.sp$bathy),]
#   # 
#   # # create depth bins and add to dataframe
#   # test.sp[test.sp$bathy <= 5, ]$depth.bin                       <- "Min - 5m"
#   # test.sp[test.sp$bathy > 5 & test.sp$bathy <= 10, ]$depth.bin  <- "5m - 10m"
#   # test.sp[test.sp$bathy > 10 & test.sp$bathy <= 20, ]$depth.bin <- "10m - 20m"
#   # test.sp[test.sp$bathy > 20 & test.sp$bathy <= 50, ]$depth.bin <- "20m - 50m"
#   # test.sp[test.sp$bathy > 50, ]$depth.bin                       <- "50m - Max"
#   # 
#   # # split test dataframe into list of dataframes by depth.bin values (5 dataframes)
#   # test.list <- split(test.sp, f = test.sp$depth.bin)
#   # 
#   # # make confusion matrices for each of the df splits
#   # eval.1step <- lapply(test.list, FUN = function(test.df){
#   #   confusionMatrix(as.factor(test.df[, "pred_oneStep"]), as.factor(test.df$BType4))
#   # })
#   # 
#   # # make confusion matrices for each of the df splits
#   # eval.2step <- lapply(test.list, FUN = function(test.df){
#   #   confusionMatrix(as.factor(test.df[, "pred_twoStep"]), as.factor(test.df$BType4))
#   # })
#   # 
#   # # get accuracies 1 step
#   # accuracies.1 <- lapply(eval.1step, FUN = function(eval.bin){
#   #   as.data.frame(eval.bin$overall)[1, ]
#   # })
#   # 
#   # # get kappas 1 step
#   # kappas.1 <- lapply(eval.1step, FUN = function(eval.bin){
#   #   as.data.frame(eval.bin$overall)[2, ]
#   # })
#   # 
#   # # get accuracies 2 step
#   # accuracies.2 <- lapply(eval.2step, FUN = function(eval.bin){
#   #   as.data.frame(eval.bin$overall)[1, ]
#   # })
#   # 
#   # 
#   # # get kappas 2 step
#   # kappas.2 <- lapply(eval.2step, FUN = function(eval.bin){
#   #   as.data.frame(eval.bin$overall)[2, ]
#   # })
#   # 
#   # # make overall dataframe that has kappas and accuracies
#   # overall.df.1 <- data.frame(depth.bin = names(kappas.1), Kappa = unlist(kappas.1), row.names = NULL)
#   # overall.df.1$Accuracy <- data.frame(Accuracy = unlist(accuracies.1), depth.bin = names(accuracies.1), row.names = NULL)[, 1]
#   # overall.df.1$Method <- "One-Step"
#   # 
#   # # make overall dataframe that has kappas and accuracies
#   # overall.df.2 <- data.frame(depth.bin = names(kappas.2), Kappa = unlist(kappas.2), row.names = NULL)
#   # overall.df.2$Accuracy <- data.frame(Accuracy = unlist(accuracies.2), depth.bin = names(accuracies.2), row.names = NULL)[, 1]
#   # overall.df.2$Method <- "Two-Step"
#   # 
#   # # bind dataframe together
#   # overall.df <- rbind(overall.df.1, overall.df.2)
#   # overall.df$Region <- as.factor(region)
#   # 
#   # # add resolution as factor to df
#   # overall.df$Resolution <- "100 m"
#   # overall.df$Resolution <- as.factor(overall.df$Resolution)
#   # overall.df$Region <- as.factor(region_name)
#   # overall.df <- gather(overall.df, Metric, Value, Kappa:Accuracy)
#   
#   return(r.stack)
#   
# }
# 
# # create evaluation metrics with test data for both 1 and 2 step substrates
# performance.metrics <- function(one.step, two.step, obs.sp, region){
#   # make rasters of one and two step
#   one.step.r <- raster(one.step)
#   two.step.r <- raster(two.step)
#   
#   # rename substrate rasters
#   names(one.step.r) <- "pred_oneStep"
#   names(two.step.r) <- "pred_twoStep"
# 
#   # make raster stack of substrate models
#   r.stack <- stack(one.step.r, two.step.r)
#   
#   # subset test data
#   test.sp <- obs.sp[obs.sp$TestDat == 1, ]
#   
#   # attach predicted raster values to test points 
#   preds <- raster::extract(r.stack, test.sp, df = TRUE)
#   
#   # convert to data frame; combine observations w/preds
#   test.sp <- cbind(as.data.frame(test.sp), preds)
#   
#   # drop extra columns --> causing issue w/gather()
#   test.sp <- test.sp[, c("BType4", "pred_oneStep", "pred_twoStep")]
#   
#   # create evaluation stats from confusion matrix
#   matx.one <- confusionMatrix(as.factor(test.sp$pred_oneStep), as.factor(test.sp$BType4))
#   matx.two <- confusionMatrix(as.factor(test.sp$pred_twoStep), as.factor(test.sp$BType4))
#   
#   # use gather() to merge 2 columns into one to be used for facet in plots
#   test.sp <- gather(test.sp, key = "Pred", value = "Pred.Value", c(pred_oneStep, pred_twoStep))
#   
#   # add region name
#   test.sp$Region <- region
#   
#   return(list(test.sp, matx.one, matx.two))
#     
# }
# 
# # main plot function for regional performance metrics using plotly
# plotly.plots <- function(df){
#   p <- plot_ly()
#   p <- add_trace(p,
#                  x = df$depth.bin, 
#                  y = df$Value, 
#                  color = df$facet, 
#                  showlegend = TRUE, 
#                  text = paste("Percent Change by Resolution:", df$per_change),
#                  colors = "RdBu", alpha = .6,
#                  name = paste0(df$Region, ": ", df$facet, " ", df$Metric)) %>%
#     # layout(yaxis = list(title = df$Metric),
#     #        xaxis = list(title = df$Region))
#     layout(yaxis = list(range = c(0, 1.0),
#                         title = df$Metric),
#            xaxis = list(title = df$Region))
#   return(p)
# }
# 
# # subplot - split df by metric and subplot
# metric.subplot <- function(df){
#   df %>%
#     mutate(facet = paste0(Resolution, " m ", Method)) %>%
#     split(.$Metric) %>%
#     lapply(plotly.plots) %>%
#     subplot(nrows = NROW(.), titleX = TRUE, shareY = TRUE, shareX = TRUE)
# }
# 
# jitter_box <- function(df){
#   # sample 20% of the data to plot
#   df <- df[sample(nrow(df), (nrow(df) * .2 )), ]
#   # jitter plot sample test data
#   p <- qplot(as.factor(BType4), 
#              as.factor(Pred.Value), 
#              data = df,
#              colour = as.factor(BType4), 
#              facets = . ~ Pred,
#              geom = c("boxplot", "jitter"), 
#              main = paste0(df$Region, ": predicted vs. observed in 20% of validation data"), 
#              xlab = "Observed Class", 
#              ylab = "Predicted Class",
#              alpha = .1) + theme_bw()
#   
#   ggsave(filename = file.path("D:/Projects/DeepWaterBackgroundSubstrate/Documents/_SubstratePackage", paste0(df$Region, "_boxplot.png")), plot = p, device = "png")
# }
