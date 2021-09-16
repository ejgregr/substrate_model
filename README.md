# Random Forest Substrate Model at Two Resolutions

__Created:__      2020/02/20 Forked from earlier version of the project created by Cole Fields
__Main author:__  Edward Gregr  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        SciTech Environmental Consulting   
__Location:__     Vancouver, BC   
__Contact:__      e-mail: ed@scitechconsulting.com | tel: 604-612-8324  
__Last Update:__  2021/09/16   
__Version:__      R version 3.6.3 x64

- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Contents and methods](#contents-and-methods)
- [Data management](#data-management)
- [Requirements](#requirements)
- [Caveats](#caveats)
- [Uncertainty](#uncertainty)
- [Acknowledgements](#acknowledgements)
- [References](#references)

## Objective
This project was created to develop coast-wide predictions of bottom type for Pacific Canada at two resolutions (20 and 100 m) for all Pacific Canadian waters. The 100 m model covers the entire shelf area, the 20 m models extend from the high water line to 50 m depth, and are constructed for 5 different subregions. The predictive performance of the models is tested using independent data sets. The work is presented in Gregr et al. 2021. 

Gregr, EJ, Haggarty D, Davies, S, Fields, C, and Lessard, J. 2021. Comprehensive marine substrate classification applied to Canada’s Pacific shelf. PLOS One. 

## Summary   

The analysis combines point observations of substrate with a suite of bathymetric and energy-based predictors to predict 4 classes of substrate: Rocky, Mixed, Sand, and Mud using a random forest modelling approach. This project includes code to load and prepare the data, develop and summarize a number of random forest models, and test the models using indepenendent substrate observations. The R Markdown script generates most of the results that appear in Gregr et al. 

This archived version of the project does not include the original, but the analysis can be re-produced using the data contained in the associated RData files. The RData files include data structures of the observations used to build test the models with the 20 m and 100 m predictor variables attached. This allows the analysis to be reproduced, but does not support reproduction of the predictions, because several of the predictor layers require data licensing agreements. Please contact the authors to arrange use of these data if required. 

Substrate observations were assembled from field data collected by Fisheries and Oceans Canada (DFO) and Natural Resources Canada. The observations were pushed through 20m and 100 m raster stacks of gridded environmental predictors to assign predictor values to each observation, thereby creating the analytic data sets included herein. For the 20 m models, polygons outlining the regional boundaries were used to partition the observations. To build the models, a large set of almost 200,000 historic observations were included in a 'Build' data set, which was divided into training and testing partitions. Three more recent independently collected sets of substrate observations were used to test predictive power. Details on the observations and predictors are provided in Gregr et al. 

Running the main script (IDE_Main.R) will variously load or re-build the data andthe analysis depending on the flags set at the top of the script. The archived version of the code is configured to load prepared source data and re-build the analysis. This will take some time. File path and program constants are contained in the substrate_functions.R script. Ensure the rData files are available and pathed correctly.

Running the RMarkdown script will create an MS Word document containing most of the results (the heatmaps are produced in the IDE_Main.R script). The project runs in the global environment so the most straightforward way to run the RMD script is from the console (>) as:

x <- substr(  Sys.time(), 1, 10)   
rmarkdown::render( "substrate_figures.Rmd",   
  output_file = paste0('SubstrateFigures_', x ),  
  output_dir = results.dir )  

## Status
Complete and finalized for distribution. 

## Contents and methods
The repository contains the R code needed to reproduce the analysis, the archived data sets needed for the analysis, and the substrate maps resulting from the analysis. The technical details of running the scripts are describede herein. See Gregr et al. for a broader  description of methods. 

The high level analysis is controlled by the *IDE_Main.R* script. It sources all necessary functions and libraries, and controls the loading of the data and building the analysis via user settings (flags) by sourcing the main supporting script *build_substrate.R*. Supporting functions for data loading, data analysis, and result building are in *substrate_functions.R*. Results are summarized by *substrate_figures.RMD*, which calls plot functions contained in *Plot_Functions.R*. 

### Data loading
Source data include the build data (points) with the 100 m predictors attached, rasters of regional 20 m predictors, and the independent data (Dive, Cam, ROV) for evaluation as points.

Data structures resulting from the data load include: *obs.100mIV, obs.20mIV, Dive.20mIV, Cam.20mIV, ROV.20mIV*. These processed substrate observations are used for the analysis and are included in the distributed code package.  

### Data analysis
Analytical steps include: creating the weighted and unweighted RF models; summarizing the model fit to the build data; testing of model predictive power using the independent substrate observations; and summarizing model performance results across regions and depths.   

Weighted random forest data structures: *rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI,	rf.region.QCS, rf.region.SOG*    

Non-weighted random forest data structures: *nwrf.region.Coast, nwrf.region.HG, nwrf.region.NCC, nwrf.region.WCVI,	nwrf.region.QCS, nwrf.region.SOG*   

*build.results* is a multi-layered list including agregated and class-based results. *build.sum* contains tablulated summaries of the build results, and is used for several tables and figures. *IDE.results.wtd* and *IDE.results.nowt* are multi-layerd lists including agregated and class-based results for both the weighted, and non-weighted models. These are also used by the RMD file to produce outputs. 

## Data management  

The results of the above steps are saved as .RData files. These include:  

*loaded_data_2020-11-19.RData* contains: point.data, pgon.data, dive.20mIV, cam.20mIV, ROV.20mIV, dive.100mIV, cam.100mIV, ROV.100mIV, obs.100mIV, obs.20mIV.  

*rf_allModels_2020-11-19.RData* contains: rf.region.Coast, rf.region.HG, rf.region.NCC, rf.region.WCVI, rf.region.QCS, rf.region.SOG.  

*buildResults_2020-11-19.RData* contains: build.results.  

*nwrf_allModels_2020-11-19.RData* contains: nwrf.region.Coast, nwrf.region.HG, nwrf.region.NCC, nwrf.region.WCVI, nwrf.region.QCS, nwrf.region.SOG.  

*nwbuildResults_2020-11-19.RData* contains: build.results.nw.  

*rasterMapObjects_2020-11-20.RData* contains: map.prev (despite the name, the map objects are too large to save, but the prevalence is needed for a figure so is retreivable).  

## Requirements
Working directories are required to re-build the data. At a minimum, a source and output directory are needed. These are located near the top of the substrate_function.R script.

## Caveats
Versioning of tidyverse packages has been an issue during development, as functions continue to be depreciated. 

## Acknowledgements
See Gregr et al.

## References
See Supporting Information for Gregr et al.

