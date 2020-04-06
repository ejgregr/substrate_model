#*******************************************************************************
# Script:  sample_effect.R
# Created: 26 March 2020. EJG
#
# Purpose: Compare data across scales and possibly across regions.
#   Evaluate effect of sample size and prevalence on performance metrics. 

# Notes:
#  - Requires loaded models and observationl data from main script.
#********************************************************************************


#-----------------------------------------------
#-- Contrast predictors at 100m and 20 m depths.
#   A few extra contortions required because the data are in different structures. 
#   Also, the IDs in the 20m are not unique across regions. 

x <- rbind(
    data.frame( Region = 'HG',   train.data.20m$HG ),
    data.frame( Region = 'NCC',  train.data.20m$NCC ),
    data.frame( Region = 'WCVI', train.data.20m$WCVI ),
    data.frame( Region = 'QCS',  train.data.20m$QCS ),
    data.frame( Region = 'SOG',  train.data.20m$SOG )
)
x <- x[ unique( x$ID ), ]
y <- train.data.100m[ match( x$ID, train.data.100m$ID ), ]


#-- Bathymetry 
plot( y$bathy * -1 ~ x$bathy, xlim=c(-10, 1000), ylim=c(-10, 1000), xlab="20 m", ylab = "100 m"  )

#-- Rugosity
plot( y$rugosity ~ x$rugosity, xlab="20 m", ylab = "100 m" )

#-- Curvature
plot( y$curvature ~ x$curvature, xlim=c(-2, 2), ylim=c(-0.4, 0.4), xlab="20 m", ylab = "100 m" )

#-- Slope 
plot( y$slope ~ x$slope )



#------------------------------------------------
#-- Compare 20 m bathy across regions and ID sets.
#   Previously done for partitioned 100 m data, which showed much less variability.
#   NOTE: Not reviewed since done originally in Feb 2020.

x <- rbind( 
  data.frame( region = 'HG', type = 'Train', depth = train.data.20m$HG$bathy ),
  data.frame( region = 'HG', type = 'Dive', depth = dive.data$HG$bathy ),
  data.frame( region = 'HG', type = 'Cam', depth = cam.data$HG$bathy ),
  data.frame( region = 'HG', type = 'ROV', depth = ROV.data$HG$bathy ),
  
  data.frame( region = 'NCC', type = 'Train', depth = train.data.20m$NCC$bathy ),
  data.frame( region = 'NCC', type = 'Dive', depth = dive.data$NCC$bathy ),
  data.frame( region = 'NCC', type = 'Cam', depth = cam.data$NCC$bathy ),
  data.frame( region = 'NCC', type = 'ROV', depth = ROV.data$NCC$bathy ),
  
  data.frame( region = 'WCVI', type = 'Train', depth = train.data.20m$WCVI$bathy ),
  data.frame( region = 'WCVI', type = 'Dive', depth = dive.data$WCVI$bathy ),
  data.frame( region = 'WCVI', type = 'Cam', depth = cam.data$WCVI$bathy ),
  # No ROV on WCVI
  
  data.frame( region = 'QCS', type = 'Train', depth = train.data.20m$QCS$bathy ),
  data.frame( region = 'QCS', type = 'Dive', depth = dive.data$QCS$bathy ),
  data.frame( region = 'QCS', type = 'Cam', depth = cam.data$QCS$bathy ),
  # No ROV on QCS
  
  data.frame( region = 'SOG', type = 'Train', depth = train.data.20m$SOG$bathy ),
  data.frame( region = 'SOG', type = 'Dive', depth = dive.data$SOG$bathy ),
  data.frame( region = 'SOG', type = 'Cam', depth = cam.data$SOG$bathy )
  # No ROV on SOG
)

#-- Rotate x axis labels ... 
par(las=2)

#-- For the 20 m data ... 
boxplot( depth ~ type*region, data = x, ylim=c(0,200), ann = FALSE)
# Add the grey vertical lines
for(i in seq(0.5 , 20 , 4)){ 
  abline(v=i,lty=1, col="grey")
}

boxplot( depth ~ region, data = x, ylim=c(0,500))
boxplot( depth ~ type, data = x, ylim=c(0,400))


#---------------------------------------------------------------------
# Compare how performance stats vary with  sample size and prevalence. 
# Uses sub-sampling. 

#-- Use only a portion of the training data otherwise it takes too long ....
#   ALSO, sample size effects plateau well below 30k points. 
#x <- train.data.100m
#x <- train.data.100m[ sample( 1:nrow(train.data.100m), 1000), ]
x <- train.data.100m[ sample( 1:nrow(train.data.100m), 35000), ]


#-- TEST Sample Size ---
rm( foo )
t1 <- Sys.time()
foo <- Test.Sample.Size( x, 20, 0.7, coast.formula )
t2 <- Sys.time()
cat("Run time:", round( difftime(t2, t1, units="mins"), 2), "minutes." )
foo

#-- Extract the results. 
#   NOTE: unpacking the results looks different here than in IDE_MAin script ... ?? 
bar <- do.call( rbind.data.frame, foo[ , 'Integrated' ] )

#-- Joininng multiple runs; Some data edits ...  
bar1 <- bar
bar <- rbind( bar1[-38,], bar[ 18:20, ] )
bar$Allocation <- bar$Exchange + bar$Shift

sample.table  <- bar

bar2 <- reshape2::melt( bar[, -c(2, 3, 7, 10, 11) ], id.vars = c('N') )
sort( colMeans(bar), decreasing = T)

  
#-- Plot the results  ... 
txtsize <- 25
a <- ggplot(bar2, aes(x=N, y=value, colour=variable, shape=variable)) +
  geom_line( size = 1.5 ) +
  geom_point( size = 3 ) +
  
  theme(axis.text = element_text(size = txtsize ),
        axis.title.x = element_text(size = rel(2), face='bold'),
        axis.title.y = element_blank(),
        legend.key = element_blank(),                               # removes the shading behind legend keys
        legend.title = element_text(size = rel(1.75), face='bold'),
        legend.text = element_text(size = txtsize ),
        legend.title.align = 0.5,                                   
        legend.key.size = unit(1, 'cm')                             # legend line spacing
  ) + 
  
  scale_shape_discrete(name  = "Metric") +
  scale_colour_discrete(name  = "Metric")

#-- Fix the background ... 
a <- a + theme(panel.background = element_rect(fill = "white"))
a <- a + theme(panel.grid.major = element_line(colour = "darkgrey"))

a


#--------------------------------
#-- And with prevalence ... 

#-- Let the function pull the sample you want.  
set.seed( 42 )
# there is some sort of maximum sample effect here ... 20k fine; 30k breaks at prev=0.5 
#   23519 limit for some reason in function? 25k works. 
foo <- Test.Prevalence( train.data.100m, 25000, 0.7, coast.formula )

#-- Extract and prepare the results. 
bar <- do.call( rbind.data.frame, foo[ , 'Integrated' ] )
bar$Allocation <- bar$Exchange + bar$Shift

prev.table <- bar

bar2 <- reshape2::melt( bar[, -c(1, 3, 7, 10, 11) ], id.vars = c('Imbalance') )
#bar2 <- reshape2::melt( bar[, -c(1, 3) ], id.vars = c('Imbalance') )

txtsize <- 25
a <- ggplot(bar2, aes(x=Imbalance, y=value, colour=variable, shape=variable)) +
  geom_line( size = 1.5 ) +
  geom_point( size = 3 ) +

  theme(axis.text = element_text(size = txtsize ),
        axis.title.x = element_text(size = rel(2), face='bold'),
        axis.title.y = element_blank(),
        legend.key = element_blank(),                               # removes the shading behind legend keys
        legend.title = element_text(size = rel(1.75), face='bold'),
        legend.text = element_text(size = txtsize ),
        legend.title.align = 0.5,                                   
        legend.key.size = unit(1, 'cm')                             # legend line spacing
  ) + 
  
  scale_shape_discrete(name  = "Metric") +
  scale_colour_discrete(name  = "Metric")


#-- Fix the background ... 
a <- a + theme(panel.background = element_rect(fill = "white"))
a <- a + theme(panel.grid.major = element_line(colour = "darkgrey"))

a


save( sample.table, prev.table, file = file.path(model.dir, 'sample_test_results_29Mar2020.RData' ))




#-- Took 30 min to run 10 samples from the full train dataset. 
save( foo, file = file.path(model.dir, 'sample_size_test_largeN_Feb03.RData' ))
# Above has 19758 as first sample size. To catenate smaller values, max N should be
# 19758/0.7 = 28225. 



### Fin


