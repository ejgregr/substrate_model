# substrate_model
Repository created: 2020/02/20
Updates: Continuous thru 2020/10/06
Version: R version 3.6.3 x64

CONTACT
Edward Gregr
SciTech Environmental Consulting
ed@SciTechConsulting.com
604-612-8324

DESCRIPTION
This project uses a coast-wide data set of observations to build random forest models of bottom substrate using predictors at two resolutions (20 and 100 m) for all Pacific Canadian waters and 5 subregions. The predictive performance of the models is tested using independent data sets.

DATA SOURCES

The code is provided with data sources in rData files.  

THIS CODE

IDE_Main.R controls the loading or re-building of the analysis based on flags at the top of the script. This code is archived so that sourcing this script will load the source data, and re-build the analysis. This will take some time.

Ensure the rData files are available and pathed correctly.

Running the RMarkdown script using the following code snippet will run it so the environment is visible. An MS Word document file of tables and figures, including supplemental materials on the analysis, will be produced. 


--------------
UPDATE HISTORY

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
- Use ConfusionMatrix() for consistency of statistics. 
   Always order parameters as (predicted, observed)
- ranger() has inbag.counts - do we need to calculate our own OOB error?
- Standardize SIGNS on the 100 m and 20 m bathymetries. Currently flipping sign on 20 m during data load.
- 2020/04/09: sample_effect.R DEFERRED UNTIL SUBSTRATE PAPER COMPLETE




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

