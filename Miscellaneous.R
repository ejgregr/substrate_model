#*******************************************************************************
# Script:  miscellaneous.R
# Created: July 2020. EJG

# A collection of code to test for various things that we don't want
# cluttering up the main scripts.


#----------------------------------
# 2020/11/23: Removing more unneeded functions including:
#   Build.IDE.BoP.Results, Partition.By.Density, Plot.ClassStats.IDE


#-----------------------------------------------------
# Test how regional BoP models perform against the IDS
# models Input: BoP geodatabase and layer, the corresponding region name in the region shape file.
# Requires: IDS point data to exist ... 
# Returns: table of how well specified BoP fc predicted all IDS
Build.IDE.BoP.Results <- function( bop, lyr, nm ){
  
  # load the regional bottom patches ... 
  bp <- file.path("d:/spacedata2019/BoPs/Delivered/", bop)
  bops <- readOGR(dsn = bp,layer = lyr )
  
  # load the region file to pull the IDS points ... 
  pgons <- pgon.data$Regions
  
  # Get name of region pgon ... this call took 30 min to put together! :\
  # 2020/07/22: factor order messed up. just pass the name in
  #idx <- which( bioregions %in% strsplit(bop, '_')[[1]][1] )
  #nm <- pgons$Region[order( pgons$Region )][ idx ]
  
  bop.IDE <- NULL
  for (i in c('Dive', 'Cam')) {
    print( i )
    
    # Select just the IDS points in the region.
    #[ Not working!! ]
    #length( point.data$Dive )
    
    #x <- pts[ pgons[ pgons$Region == "Haida Gwaii",], ]@data
    
    rpts <- point.data[[ i ]][ pgons[ pgons$Region == nm,], ]
    print( length( rpts) )
    
    # Ensure projection of BoPs agree with points (can be subtly different) ...
    crs( bops ) <- crs( rpts )
    
    # Now pull BTypes to the points 
    y <- over( rpts, bops[ c('BType1', 'BType2') ] )
    
    # drop the spatial deets and combine obs with BoP pred 
    z <- cbind(rpts@data, y)
    
    z <- z[ !is.na(z$BType1), ]
    # Transform BTypes into new BType comparable w RMSM
    z$BType4 <- with(z, ifelse( BType1 != 3, BType1, 
                                ifelse( BType2 == 'a', 3, 4))
    )
    
    z$RMSM  <- factor( z$RMSM, levels = c('1','2','3','4') )
    z$BType4 <- factor( z$BType4, levels = c('1','2','3','4') )
    
    #print( caret::confusionMatrix( z$RMSM, z$BType )$table )
    bop.IDE <- rbind( bop.IDE, 
                      cbind( 'IDS' = i, 'Region' = strsplit(bop, '_')[[1]][1],
                             Results.Row( z$RMSM, z, paired = T )$Integrated ))
  }
  return( bop.IDE )
}


#---------------------------------------------
# Partition input points into low and high density piles
# Returns: List of 2 dataframes, one for each of low/hi density regions
# Requires: names.100m as a global (defined during data load)
#           Obs points, and the density shape file. 
Partition.By.Density <- function( pts ){
  
  out <- list()
  
  #-- Load the regions shape file containing the High Density pgons ... 
  pgons <- readOGR( file.path(source.dir, "/regions/hi_density_area.shp") )
  
  # Split the points into those that fall in the High density area ... 
  dens.pts <- x[ pgons[ pgons$name == "High_Density_Area",], ]@data
  # and those that don't ... (confirmed IDs are unique)
  sparse.pts <- x[ !(x$ID %in% dens.pts$ID), ]@data
  
  # Now fix the attribute lists 
  dens.pts$BType4   <- as.factor( dens.pts$BType4 )
  sparse.pts$BType4 <- as.factor( sparse.pts$BType4 )
  
  x <- dens.pts[ , !colnames(dens.pts) %in% drop.list ]
  names(x)[3:12] <- names.100m 
  out <- c( out, list( 'Dense' = x) )
  
  x <- sparse.pts[ , !colnames(sparse.pts) %in% drop.list ]
  names(x)[3:12] <- names.100m 
  out <- c( out, list( 'Sparse' = x) )
  
  
  return( out )
}


#----------------------------------
# Nov 2020: Seeking the best approach to showing class-based metrics for IDE. 
# This is class, by IDE, by region (Coast, HG, NCC), by metric. 
# Tried to match both Fig 2 heatmaps - SD says it looks like a quilt.
# Also tried the accuracy variables but Reliability doesnt have a baseline and TSS
#   can't be calculated (I don't think) for classes. 

# Right now, the temporary figure uses only baseline differenced TPR/TNR for illustration. 
# 2020/11/23: Dropped this idea as the figure doesn't seem to tell us anything useful. 
#   Dropped RMD code is below ... 

a <- IDE.results.wtd$PerClass
b <- a[ a$Stat   %in% c('TPR', 'TNR', 'User'), ]
c <- b[ b$Region %in% c('Coast', 'HG', 'NCC'), ]

# Take apart, and put back together after adjusting for different random baselines

d <- rbind( 
  cbind( 
    c[ c$Stat == 'TPR', c('Region', 'IDS')], 'Stat' = 'Accuracy',  
    c[ c$Stat == 'TPR', c('1', '2', '3', '4') ] - 0.25 ),
  cbind( 
    c[ c$Stat == 'TNR', c('Region', 'IDS', 'Stat')],  
    c[ c$Stat == 'TNR', c('1', '2', '3', '4') ] - 0.75 ),
  cbind( 
    c[ c$Stat == 'User', c('Region', 'IDS')], 'Stat' = 'Reliability',
    c[ c$Stat == 'User', c('1', '2', '3', '4') ])
)

d
str(d)
Plot.ClassStats.IDE( d[ d$IDS == 'Dive', -grep('IDS', colnames(c)) ], 'Dive', pal.cb3b )
Plot.ClassStats.IDE( d[ d$IDS == 'Cam', -grep('IDS', colnames(c)) ], 'Cam', pal.cb3b )
Plot.ClassStats.IDE( d[ d$IDS == 'ROV', -grep('IDS', colnames(c)) ], 'ROV', pal.cb3b )
row.names(d) <- NULL

d <- IDE.results.wtd$PerClass[ , c('IDS', 'Stat', 'Region', '1', '2', '3', '4')]
d <- d[ d$Stat   %in% c( 'TPR', 'TNR', 'User'),]
d <- d[ d$Region %in% c( 'Coast', 'HG', 'NCC'),]
row.names( d ) <- NULL
d[ 22, 5] <- 0

d[ d$IDS == 'Dive',]

str(d)
str(build.sum$build.results.ByClass)

Heat.Build.Class.Stats( d[ d$IDS == 'Dive', -1], T, 'TPR', rev( pal.heat.10 ), 800, 600, 'black', 'Dive' )
Heat.Build.Class.Stats( d[ d$IDS == 'Dive', -1], T, 'TNR', rev( pal.heat.10 ), 800, 600, 'black', 'Dive' )
Heat.Build.Class.Stats( d[ d$IDS == 'Dive', -1], T, 'User', rev( pal.heat.10 ), 800, 600, 'black', 'Dive' )

Heat.Build.Class.Stats( d[ d$IDS == 'Cam', -1], T, 'TPR', rev( pal.heat.10 ), 800, 600, 'black', 'Cam' )
Heat.Build.Class.Stats( d[ d$IDS == 'Cam', -1], T, 'TNR', rev( pal.heat.10 ), 800, 600, 'black', 'Cam' )
Heat.Build.Class.Stats( d[ d$IDS == 'Cam', -1], T, 'User', rev( pal.heat.10 ), 800, 600, 'black', 'Cam' )

Heat.Build.Class.Stats( d[ d$IDS == 'ROV', -1], T, 'TPR', rev( pal.heat.10 ), 800, 600, 'black', 'ROV' )
Heat.Build.Class.Stats( d[ d$IDS == 'ROV', -1], T, 'TNR', rev( pal.heat.10 ), 800, 600, 'black', 'ROV' )
Heat.Build.Class.Stats( d[ d$IDS == 'ROV', -1], T, 'User', rev( pal.heat.10 ), 800, 600, 'black', 'ROV' )


# RMD Code to do the simple version of the above ...

```{r FigXXa_IDEbyClass, echo=FALSE, out.width='120%', fig.fullwidth=TRUE}
# A simple region-based facet of class-based accuracy and TNR. 
# Build the data, then 3 plots.

a <- IDE.results.wtd$PerClass
b <- a[ a$Stat   %in% c('TPR', 'TNR'), ]
c <- b[ b$Region %in% c('Coast', 'HG', 'NCC'), ]
row.names(c) <- NULL
colnames(c)[ colnames(c) %in% c('1', '2', '3', '4') ] <- c('Hard', 'Mixed', 'Sand', 'Mud')

# Take apart and put back together after adjusting for different random baselines

d <- rbind( 
  cbind( 
    c[ c$Stat == 'TPR', c('Region', 'IDS')], 'Stat' = 'Accuracy',  
    c[ c$Stat == 'TPR', c('Hard', 'Mixed', 'Sand', 'Mud') ] - 0.25 ),
  cbind( 
    c[ c$Stat == 'TNR', c('Region', 'IDS', 'Stat')],  
    c[ c$Stat == 'TNR', c('Hard', 'Mixed', 'Sand', 'Mud') ] - 0.75 )
)

Plot.ClassStats.IDE( d[ d$IDS == 'Dive', -grep('IDS', colnames(c)) ], 'Dive', pal.cb2, sz = 25 )
Plot.ClassStats.IDE( d[ d$IDS == 'Cam', -grep('IDS', colnames(c)) ], 'Cam', pal.cb2, sz = 25 )
Plot.ClassStats.IDE( d[ d$IDS == 'ROV', -grep('IDS', colnames(c)) ], 'ROV', pal.cb2, sz = 25 )


#2020/09/04: New plot to examine class-based stats across 3 RF models and all regions
# Needs to be applied once per statistic (accuracy, specificity, and reliability)
#2020/11/05: Adapted to do TPR, TNR, Reliability by class faceted by region. Done for each IDS. 
Plot.ClassStats.IDE <- function( dat.table, ylab, apal, sz=20, lx=0, ly=0 ){
  
  foo <- melt( dat.table, id.var = c('Region', 'Stat'))
  
  a <- foo %>%
    # Adjust levels for correct faceting ... 
    mutate(Region = factor(Region, levels=c("Coast", "HG", "NCC", "WCVI", "QCS", "SOG"))) %>%
    mutate(Stat = factor(Stat, levels=c("Accuracy", "TNR", "Reliability"))) %>%
    
    ggplot(aes(x = variable, y = value, fill = Stat)) +
    geom_bar(stat = "identity", width = .8, position = "dodge") +
    labs( x = NULL, y = ylab ) +
    facet_grid(. ~ Region) +
    scale_fill_manual(values = apal) +
    theme_bw() +
    theme(  text = element_text(size=sz ),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
            
            # legend stuff          
            #            legend.position = c(lx, ly),
            #            legend.background = element_rect(fill="gray90", size=1, linetype="dotted")
    )
  return(a)
}




#------------------------------------------
#-- Testing TSS Calculation ...
test <- matrix( c(10,0,0,0, 0,10,0,0, 0,0,10, 0, 0, 0, 0, 10), 4, 4 )
test <- matrix( c(.25,0,0,0, 0,.25,0,0, 0,0,.25, 0, 0, 0, 0, .25), 4, 4 )

#-- A simplee bad model ...
test <- matrix( c(0,10,0,0, 0,0,0,10, 10,0,0,0, 0,0,10,0), 4, 4 )

#-- Shouldn't these have greater badness?
test <- matrix( c(0,100,0,0, 0,0,0,10, 10,0,0,0, 0,0,10,0), 4, 4 )
test <- matrix( c(0,33,33,34, 10,0,10,10, 10,10,0,10, 10,10,10,0), 4, 4 )

test
TSS.Calc(test)

#-- Random matrix has values around 0 which is good ...
test <- matrix( sample(1:100, 16), 4, 4)
TSS.Calc(test)


#--------------------------------------
#-- Testing diffeR stats ... 
z <- caret::confusionMatrix( rf.region.HG$predictions, train.data.20m$HG$BType4)
diffr.Stats( z$table )

x <- matrix( z$table, 4,4 )
rownames( x ) <- c(1,2,3,4); colnames( x ) <- c(1,2,3,4)
y <- diffTablej( x, digits = 0, analysis = "error" )


#-------------------------------------------
a <- build.sum$build.results.ClassPrev
b <- nw.build.sum$build.results.ClassPrev[7:12,]

# change numbers to proportions
x <- a[, names(a) %in% c('Hard', 'Mixed', 'Sand', 'Mud')]
y <- data.frame(  x / rowSums( x ))

# add results back to the Region column
z <- cbind( 'Region' = a$Region, 'Source' = c(rep('Observed',6), rep('Weighted',6)), y)

#--- do it again with the non-wtd model, this time keep only the predictions ... 

x <- b[, names(b) %in% c('Hard', 'Mixed', 'Sand', 'Mud')]
y <- data.frame(  x / rowSums( x ))
zz <- cbind( 'Region' = b$Region, 'Source' = rep('No Weights',6), y)

foo <- rbind( z, zz)


Plot.Obs.Pred.Prevalence.Build( foo, pal.cb3b, sz=30, lx=0.95, ly=0.875 ) 



#---------------------------
# Build some summary tables.

table( obs.100mIV$BType4 )

x <- point.data$Obs@data
table( x$BT_Sorc )


w <- NULL
for (y in c('CHS_Grabs', 'NRCan_Observations', 'Dive', 'ROV', 'Marsh') ) {
  
  v <- table( x[ x$BT_Sorc == y,]$BType4 )
  w <- rbind( w, v )
  
}
w <- cbind( 'Source' = c('CHS_Grabs', 'NRCan_Observations', 'Dive', 'ROV', 'Marsh'), 
            data.frame(w) )
w
colSums(w[, -1])
sum( colSums(w[, -1]) )

#-------


x <- rbind( table( point.data$Dive$RMSM ),
            table( point.data$Cam$RMSM ),
            table( point.data$ROV$RMSM) )

x <- cbind( 'Source' = c('Dive', 'Cam', 'ROV'), 
            data.frame( x) )
x
colSums(x[, -1])
sum( colSums(x[, -1]) )


#----- Pull some results for display ... 

# get the obs test points
a <- obs.100mIV
a <- a[ a$TestDat == 1, ]

b <- predict( rf.region.Coast, a )

# made a df with ID and prediction, and error
c <- cbind('ID' = a$ID, 'PredCst' = b$predictions, 'Err_Cst' = b$predictions == a$BType4 )

# pull relevant point data ... 
d <- point.data$Obs
d <- d[ d$TestDat == 1, ]

# attach the coast results ... 
e <- cbind( d, c)

#---- Get NCC results
a <- obs.20mIV$NCC
a <- a[ a$TestDat == 1, ]

# predict and build result 
b <- predict( rf.region.NCC, a )
c <- cbind('ID' = a$ID, 'PredNCC' = b$predictions, 'Err_NCC' = b$predictions == a$BType4 )

# further subset the points ... 
f <- e[ e$ID %in% c[,'ID'], ]


f <- cbind( f, c )

dim(f)[[1]]
sum( f$Err_Cst == 1 ) / dim(f)[[1]]

sum( f$Err_NCC == 1 )/ dim(f)[[1]]

sum( f$Err_Cst != f$Err_NCC )


writeOGR(obj=f, dsn=results.dir, layer="obspreds", driver="ESRI Shapefile", overwrite_layer = TRUE) # this is in geographical projection


rm('a', 'b', 'c', 'd', 'e', 'f')


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

y <- Partition.By.Density( point.data$Obs[ point.data$Obs$TestDat == 1, ] )

#-- Build the table piece ... 
compare.what <- data.frame( 'Region' = 'Coast', 'IDS' = 'Dive' )
w <- Results.Row( rf.region.Coast, x.test )
x <- cbind( compare.what, w$Integrated )
results.int <- rbind( results.int, x ) 


#----------------------------------------
z <- IDE.depths[ IDE.depths$Region == 'HG', c( 'IDS', 'Ribbon', 'Accuracy', 'Shift','Exchange','Quantity' )]


#----------------------------------------------------------------
# This code confirms the base rates for TPR and TNR, solidifying
# the difference between Accuracy and TPR in a multiclass context.

w <- obs.20mIV$HG[ obs.20mIV$HG$TestDat ==1, ]
x <- predict( rf.region.HG, w )
y <- caret::confusionMatrix( x$predictions, w$BType4 )
z <- diffr.Stats( y$table )

z <- data.frame(y$byClass)
z$Balanced.Accuracy * z$Prevalence

# Accuracy ... 
y$table / sum( y$table )

.248+.055+.336+.09
# 
# x <- predict( rf.region.HG, w )
# y <- caret::confusionMatrix( x$predictions, w$BType4 )
str( y$byClass )
mean( y$byClass[, 11 ] )
weighted.mean( y$byClass[, 11 ], y$byClass[, 8 ])
y$overall


# generate a random sample wtd by prevalence
length( w$BType4 )
x <- sample( as.factor( c(1,2,3,4)), length( w$BType4 ), replace=T, prob = y$byClass[, 8 ])
y <- caret::confusionMatrix( x, w$BType4 )
y


set.seed(42)
x1 <- sample( as.factor( c(1,2,3,4)), 10000, replace=T )
set.seed(24)
x2 <- sample( as.factor( c(1,2,3,4)), 10000, replace=T )
y <- caret::confusionMatrix( x1, x2 )
y$table / sum( y$table )
y$overall

z <- as.data.frame( y$byClass ) 
sum( z$Sensitivity * z$Prevalence ) + 0.25 
sum( z$Specificity * z$Prevalence ) - 0.25

sum( z$Sensitivity * z$Prevalence)



#-- NOTE: Balanced accuracy is corrected for prevalence. 

sum( z$`Balanced Accuracy` * z$Prevalence )
z$Specificity


#---------------------------------------------------------------------------------
# RMD code I'm not ready to throw out. Flensed after review of Sept 2020 draft ... 

```{r IDE_byClass_allRegions_table, echo=FALSE}
library( knitr )

T.cap <- "Table S2: Class-based assessment of independent data evaluation showing true positive (TPR/Accuracy) and true negative (TNR/Specificity) rates, and reliability (User Accuracy) for all regions, for all independent data sets."

x <- rbind(
  cbind( 'Model' = 'Wtd',   IDE.results.wtd$PerClass[ IDE.results.wtd$PerClass$IDS == 'Dive', ]),
  cbind( 'Model' = 'NoWts', IDE.results.nowt$PerClass[ IDE.results.nowt$PerClass$IDS == 'Dive', ]),
  cbind( 'Model' = 'Trim',  IDE.results.trm$PerClass[ IDE.results.trm$PerClass$IDS == 'Dive', ]),
  cbind( 'Model' = 'Wtd',   IDE.results.wtd$PerClass[ IDE.results.wtd$PerClass$IDS == 'Cam', ]),
  cbind( 'Model' = 'NoWts', IDE.results.nowt$PerClass[ IDE.results.nowt$PerClass$IDS == 'Cam', ]),
  cbind( 'Model' = 'Trim',  IDE.results.trm$PerClass[ IDE.results.trm$PerClass$IDS == 'Cam', ])
)

x <- x[ x$Stat %in% c( 'TPR', 'TNR', 'User'), ]
colnames(x) <- c('Model', 'Region', 'IDS', 'Stat', 'Hard', 'Mixed', 'Sand', 'Mud')

y <- x[ , c('IDS', 'Region', 'Stat', 'Model', 'Hard', 'Mixed', 'Sand', 'Mud') ]

#y <- mutate(y, Scale = factor(Stat, levels=c("Points", "Map")))
y <- mutate( y, Stat = forcats::fct_relevel( y$Stat, c('TPR', 'TNR', 'User') ))

z <- y[ with(y, order(IDS, Region, Stat, Model)), ]
row.names(z) <- NULL

kable( z, digits=2, caption = T.cap )

```   


```{r Model_comparison_integrated_table, echo=FALSE}
library( knitr )

T.cap <- "Table S3: Integrated performance metrics for three models (the full random forest (RF) model, a RF model trimmed to use only the 6 most influential predictors, and the bottom patch (BoP) object-based model) against the Dive and Camera (Cam) independent data sets, for all 20 m regional models." 

# USES model.compare built above for Fig 5
kable( model.compare, digits=2, caption = T.cap )

```   
```{r Fig6_Model_Compare, echo=FALSE, out.width='120%', fig.fullwidth=TRUE}
#-- Does using just the top 6 predictors improve RF performance against the IDS?
#   Rolled up with BoPs data as figure

#-- COMBINE the alternate models into a table of IDE results ... 
a <- IDE.results.wtd$Integrated
a <- a[, c('IDS', 'Region', 'N', 'Imbalance', 'TSS', 'Accuracy', 'TNR', 'Quantity', 'Exchange', 'Shift') ]
a <- a[ a$IDS %in% c('Dive', 'Cam'), ]

a <- data.frame( a[ a$Region %in% c('HG', 'NCC', 'WCVI', 'QCS', 'SOG'), ])
a <- cbind( a[, c(1,2)], 'Model' = 'Full', a[, c(3:10)] )

b <- IDE.results.trm$Integrated
b <- b[, c('IDS', 'Region', 'N', 'Imbalance', 'TSS', 'Accuracy', 'TNR', 'Quantity', 'Exchange', 'Shift') ]
b <- b[ b$IDS %in% c('Dive', 'Cam'), ]

b <- data.frame( b[ b$Region %in% c('HG', 'NCC', 'WCVI', 'QCS', 'SOG'), ])
b <- cbind( b[, c(1,2)], 'Model' = 'Trimmed', b[, c(3:10)] )

c <- IDE.BoP[, c('IDS', 'Region', 'N', 'Imbalance', 'TSS', 'Accuracy', 'TNR', 'Quantity', 'Exchange', 'Shift') ]
c <- cbind( c[, c(1,2)], 'Model' = 'BoP', c[, c(3:10)] )

d <- rbind( a, b, c )

model.compare <- d[ with(d, order(IDS, Region)), ]
row.names(model.compare) <- NULL


Plot.TSS.By.IDS.For.Regions( model.compare, pal.cb3b, sz = 25, lx=0.85, ly=0.87 )

```

**Figure 5: True skill statistic (scaled to a random baseline of 0.5) assessing the forecast skill of three models (the full random forest (RF) model, a RF model trimmed to use only the 6 most influential predictors, and the bottom patch (BoP) object-based model) against the Dive and Camera (Cam) independent data sets, for all 20 m regional models.  **  
  


  ```{r Fig7_IDS_Two, echo=FALSE, out.width='120%', fig.fullwidth=TRUE}

# What about TSS vs. Quantity and Allocation?
# Remember: Accuracy = 1 – (Quantity + Exchange + Shift).

y <- IDE.results.wtd$Integrated
z <- y[ y$Region %in% c('Coast', 'HG', 'NCC'), ]

z$Allocation <- z$Shift + z$Exchange
z <- z[, c( 'Region', 'IDS', 'Accuracy', 'Quantity', 'Exchange','Shift' )]

Plot.Pontius.By.IDS.For.Regions( z, rev(pal.cb4), sz = 25 )
```

**Figure 7: Error assessment statistics for independent data sets by region.**  
  


#-------------------------
# Playing with colours ...  

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

  
  
