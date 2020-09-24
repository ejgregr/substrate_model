# substrate_model
Repository created: 2020/02/20
Updates: Continuous thru 2020/09/03
Version: R version 3.6.2 x64

CONTACT
Edward Gregr
SciTech Environmental Consulting
ed@SciTechConsulting.com
604-612-8324

DESCRIPTION
This project uses a coast-wide data set of observations to build random forest models of bottom substrate using predictors at two resolutions (20 and 100 m) for all Pacific Canadian waters and 5 subregions. The predictive performance of the models is then tested using independent data sets.

To load/run the model:

1) From IDE_Main.R, source substrate_functions.R and Plot_Functions.R
2) Load the 4 RData files. 
3) Re-build the summary data structures
4) Run the RMarkdown script (from the console so the environment is visible)


UPDATES

Sept 2020
Been awhile. Rmd now being used to organise all R-generated tables and figs. And also to link in figs (pngs) made elsewhere.
Loading or re-building of analytic data frames now automated. Just run IDE_Main.R and go ... hopefully. 
Added analysis of non-weighted model, and IDS evaluation of BoPs.
Writing to CSV files removed. 


June 2020 (b)
Added rMarkdown file; moved all plotting here to organize figures and tables. 
Split IDE_Main.R into a control script and a data and summary building script (build_substrate.R)
Subsumed depth_effects.R script into the main function script (substrate_functions.R) and the markdown file. 
More plots added; Palettes turned CB friendly and standardized.

June 2020 (a)
Getting close to finalizing figures. Integrated code from Cole to generate and plot predictions. Having trouble plotting these from the saved results. Other than that, just need to adjust/standardize colour palettes, and tweak some plots.

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
- Use ConfusionMatrix() for consistency of statistics. 
   Always order parameters as (predicted, observed)
- ranger() has inbag.counts - do we need to calculate our own OOB error?
- Standardize SIGNS on the 100 m and 20 m bathymetries. Currently flipping sign on 20 m during data load.

-------------------------------------------------------------------------------
Description of data processing/analysis

IDE_Main.R

Part 1: Build Data Sets
obs.data		All 4 observational data sets, loaded from the geodatabases, with attributes attached.
train.data.100m 	The obs data with the predictors already assigned. Read from gdb file. (n=197587)
train.data.20m		The obs points partitioned by region, with the 20 m predictors attached. (list of 5)
dive.data		One of thee IDS with 20 m predictors attached.
cam.data		One of thee IDS with 20 m predictors attached.
ROV.data		One of thee IDS with 20 m predictors attached.

Names on the train.data.100m attributes updated to match the more complete names from the 20m raster stack. 
This allows the built models to predict with different data sets. 

Part 2: Build and evaluate 6 RF models
rf.coast	Built using train.data.100m
rf.region.HG	Built using train.data.20m$HG
rf.region.NCC	Built using train.data.20m$NCC
rf.region.WCVI	Built using train.data.20m$WCVI
rf.region.QCS	Built using train.data.20m$QCS
rf.region.SOG	Built using train.data.20m$SOG

Make.Range.Model() returns: the ranger() model structure; Build statistics; and the prevalence of observed and predicted values of the testing partition.
This function also supports re-sampling of the training data to create confidence intervals but doesn't look like things are going that way.

