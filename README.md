# substrate_model
Repository created: 2020/02/20
Last Update: 2020/03/31

CONTACT
Edward Gregr
SciTech Environmental Consulting
ed@SciTechConsulting.com
604-612-8324

DESCRIPTION
This project uses a coast-wide data set of observations to build random forest models of bottom substrate using predictors at two resolutions (20 and 100 m) for all Pacific Canadian waters and 5 subregions. The predictive performance of the models is then tested using independent data sets.

R version 3.6.2 x64

UPDATES

April 2020
Add some of Cole's graphs to inform results. 
  Needs a re-write of Make.Ranger.Model() to produce obs/pred pairs during testing
Re-structure scripts to better separate data build, summarize/plot, and report (reflecting Cole's structure)

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
