# substrate_model
Repository created 2020/02/20.

CONTACT
Edward Gregr
SciTech Environmental Consulting
ed@SciTechConsulting.com
604-612-8324

DESCRIPTION
This project uses a coast-wide data set of observations to build random forest models of bottom substrate using predictors at two resolutions (20 and 100 m) for all Pacific Canadian waters and 5 subregions. The predictive performance of the models is then tested using independent data sets.

R version 3.6.2 x64

UPDATES

March 2020
Replace model building section with the building/fitting/reporting code. Save testing results to file along with models. 
Review statistics used to evaluate models.
Finalized model fitting report for the principal models (100 m Coastwide; 20 m regional (N=6)).
Examine model performance by depth class. New script = depth_effect.R.
Formalized the testing of whether/how statistics change with prevalence/sample size.

February 2020
Intial repo commit includes a control file (IDE_main.R) and a constants and functions file (Sediment_Functions.R). The control file runs to completion. It imports and structures data, creating models, and tabulates some intial independent data evaluation results. 
NOTE: Data load and model build can be time consuming. Separate local folders (/DATA, /MODELS, /RESULTS) are used to store the most recent versions of these.


