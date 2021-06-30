# Random Forest Substrate Model at Two Resolutions

__Created:__      2020/02/20 from earlier version of project created by Cole Field
__Main author:__  Edward Gregr  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        SciTech Environmental Consulting   
__Location:__     Vancouver, BC   
__Contact:__      e-mail: ed@scitechconsulting.com | tel: 604-612-8324
__Last Update:__  2020/11/09   
__Version:__      R version 3.6.3 x64

- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Contents](#contents)
  + [Subsections within contents](#subsections-within-contents)
- [Methods](#methods)
  + [Subsections within methods](#subsections-within-methods)
- [Requirements](#requirements)
- [Caveats](#caveats)
- [Uncertainty](#uncertainty)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
This project uses a coast-wide data set of observations to build random forest models of bottom substrate using predictors at two resolutions (20 and 100 m) for all Pacific Canadian waters and 5 subregions. The predictive performance of the models is tested using independent data sets.

## Summary
Substrate observations were assembled in ArcGIS from field data collected by Fisheries and Oceans Canada. The observations were pushed through raster stacks of predictor variables. We used 20m and 100 m raster stacks (collections of gridded environmental values) comprised of bathymetric and energy-based possible predctors to predict 4 classes of substrate: Rocky, Mixed, Sand, and Mud. Observations include the 'Build' data set, divided into training and testing partitions, and three independent data sets collected for validation. 

We imported the ArcGIS raster stacks and used them to 1) estimate the model, and 2) interpolate it across the study area. To build the random forst model, we assigned raster values to each observations. This analytic data set is disributed with this script. The raster stacks contain proprietrary data, and are available on request from Fisheries and Oceans Canada. Details on the observations and predictors are provided in Gregr et al. (2021). The data to replicate the analysis are included as rData files.  

Running the main script (IDE_Main.R) will produce an RMD file containing most of the results. Some functions seem fussy within RMD (the heatmaps) so these need to be pre-generated and pointed to by the RMD file. IDE_Main.R controls the loading or re-building of the analysis using flags set at the top of the script. The script is archived so that sourcing this script will load the source data, and re-build the analysis. This will take some time. File path and program constants are contained in ... Ensure the rData files are available and pathed correctly.

Running the RMarkdown script will create an MS Word document file of tables and figures, including supplemental materials on the analysis. To run the Markdown script, Knitr needs to see the necessary data objects. As this entire project runs in the global environment, the most straightforward way to do this is to run the script from the console (>) as:   


## Status
Stablized. Last revision prior to re-submission focused on updating documentation, and clarifying the partial data pathway. 

## Contents
The repository contains the R code needed to reproduce the analysis, the archived data sets needed for the analysis, and the substrate maps resulting from the analysis.

### R code
The R code includes the following files:  
*IDE_Main.R*: The control script that sources necessary functions and libraries, and controls the loading of the data and building the analysis. 
*build_substrate.R*: A high level script that does most of the work (function calling and result building) for IDE_Main.R.   
*substrate_functions.R*: Code for loading packages, defining constants (including file directories), and defining all necessary functions for data loading, data analysis, and result building.   
*plot_functions.R*: Plotting functions for all the figures in the manuscript.   
*substrate_figures.RMD*: This Markdown script generates tables and figures from existing data objects in the environment. Some plots have supporting code chunks to prepare the summary data for plotting (e.g., scaling of the TSS) of the summarized data.   

### Data sets
Subsections within contents

### Substrate maps
Use subsections to describe the purpose of each script if warranted.


## Methods
Please see Gregr et al. (2021) for a detailed description of the source data, analytical methods, and results. 

## Requirements
List the input data requirements or software requirements to successfully execute the code.

Working directories are required to re-build the data. At a minimum, and output directory should be provided. These are located near the top of the substrate_function.R script.


## Caveats
Anything other users should be aware of including gaps in the input or output data and warnings about appropriate use.

Versioning of tidyverse packages has been an issue during development, as functions continue to be depreciated. 


## Uncertainty
*Optional section.* Is there some uncertainty associated with the output? Assumptions that were made?


## Acknowledgements
*Optional section.* List any contributors and acknowledge relevant people or institutions


## References
*Optional section.*


--------------
UPDATE HISTORY

Nov 2020
Mopping up time. Cleaning up RMD file, plotting script, and stabilizing the data packaging. 

Oct 2020
Figures and tables for manuscript finalized.
Bottom patch analysis, and the testing of predictor variable trimming (i.e., the overfitting story) dropped.
Lots of extra text removed. Spacing standardized. Testing code moved to Miscellaneous.R script.
RMD file heavily revised including: table content, and table and figure numbers updated; contingency tables added for 6 main models; Draft of supplementary material (code documentation) added; earlier development notes deleted.

Sept 2020
Rmd now being used to organise all R-generated tables and figs. And also to link in figs (pngs) made elsewhere.
Loading or re-building of analytic data frames now automated. Just run IDE_Main.R and go ... hopefully. 
Added analysis of non-weighted model, and IDS evaluation of BoPs.
Writing to CSV files removed. 
Data density analysis removed because sampling strongly correlated and hence conflated with depth
Hence, regional figures showing build data density also deemed unnecessary. 
		
June 2020 (b)
Added rMarkdown file; moved all plotting here to organize figures and tables. 
Split IDE_Main.R into a control script and a data and summary building script (build_substrate.R)
Subsumed depth_effects.R script into the main function script (substrate_functions.R) and the markdown file. 
More plots added; Palettes turned CB friendly and standardized.

June 2020 (a)
Getting close to finalizing figures. Integrated code from Cole to generate and plot predictions. Having trouble plotting these from the saved results. Other than that, just need to adjust/standardize colour palettes, and tweak some plots.
- 2020/06/15: re-structured to place all plotting in RMD file. 

May 2020 
Deep cleaning of data loading and model preparation after code review with Sarah Davies (@sare-ah). This included loading data straight from the shape files provided; centralizing most of the attribute manipulation; purged unncessary functions, data structures, and code; updated some function documentation. Any older stuff that breaks as a result would likely need to be re-written anyway. 
Model_summaries.R script dropped and the few functions absorbed in substrate.functions.R

April 2020
Update predictor influence plots to show proportion of best variable explained
Add some of Cole's graphs to inform results. 
Re-write Make.Ranger.Model() to produce obs/pred pairs of build testing partition.
Re-structure scripts to better separate model building, summarizing, and plotting.

March 2020
Replace model building section with the building/fitting/reporting code. Save testing results to file along with models. 
Review statistics used to evaluate models.
Finalized model fitting report for the principal models (100 m Coastwide; 20 m regional (N=6)).
Examine model performance by depth class. New script = depth_effect.R.
Formalized the testing of whether/how statistics change with prevalence/sample size.

February 2020
Moved and renamed fully operational version to GitHub.
Intial repo commit includes a control file (IDE_main.R) and a constants and functions file (Sediment_Functions.R). The control file runs to completion. It imports and structures data, creating models, and tabulates some intial independent data evaluation results. 
All data inspection and comparison now in IDE_data_compare_V1.R script.

NOTES: 
- Data load and model build can be time consuming. Separate local folders (/DATA, /MODELS, /RESULTS) are used to store the most recent versions of these.
   Always order parameters as (predicted, observed)
- ranger() has inbag.counts - do we need to calculate our own OOB error?
- 2020/04/09: sample_effect.R DEFERRED UNTIL SUBSTRATE PAPER COMPLETE


