# The following repository is comprised of the following code sections, intended to use in this order:
## 1. Preprocess_data based on matlab 2022 with eeglab 2022
   ### Run the main.m
   a. loads the eeg data using the index tables and seizure onset/offset times in the format provided in example data
   b. Applies the automated preprocessing pipeline
   c. saves a file containing the seizure signals and a control segment from the same file with similar length.

   This project is based on the data from the EPILEPSIAE surface database, described in this [publication](https://pubmed.ncbi.nlm.nih.gov/20863589/). 

## 2. Visual_analysis
   ### Run the label_all.m
   This file sequentially loads the files saved in the previous step and provides a visual graphical user interface to label the data
   Prior to using this step on new data it is recommended that the rater practices and understands the ICA representations provided
   A good tutorial and practice option can be found in: https://labeling.ucsd.edu/tutorial/practice

   if the user is interested in labeling bifurcation morphology it is recommended to practice on the simulation data
   this is provided with the paper: ["A taxonomy of seizure dynamotypes"](https://doi.org/10.7554%2FeLife.55632) and can be downloaded [here](https://doi.org/10.7302/ejhy-5h41)

   The information presented includes:
   ----TODO: add description



## 3. Label analysis
   ### Run the label_validation_analysis.py (tested up to python 3.12)
   #### This file will:
   * Analyze detectability rates based on the lobe containing the epileptogenic zone and seizure classification.
   * Run the inter rater agreement analysis
   * Produce the label distribution plots and tables
   * Compare the labeler selection to automated measures extracted from the signal
   * Run the selection of a single component per seizure based on: TODO: add flow chart


   ### Bayesian modeling of clinical factors, use: (tested up to R 4.4 + [brms 2.21.0](https://paul-buerkner.github.io/brms/) + cmdstanr 2.35.0) 
   1. "brms_clinical_factor_model_onset.R", to fit the model to the onset labels
   2. "model_analysis_onset_summary.R", to run onset model inference analysis
   3. brms_clinical_factor_model_offset.R, to fit the model to the offset labels
   4. model_analysis_offset_summary.R, to run offset model inference analysis

These scripts will save the model, inference result tables and plots in the "example_data" folder


## 4. Automated classification (tested up to python 3.12)
   ### Includes a file to run for each classification task
   1. Identify components with a clear transition
   2. Identify slowing of the inter spike intervals towards seizure edges
   
   ### Each file will:
   1. Train and evaluate the model in a leave-one-out cross validation approach
   2. Save performance metrics including balanced accuracy, ROC-AUC ans Cohen's Kappa
   3. Train a full model for feature contribution analysis and run [SHAP analysis](https://shap.readthedocs.io/en/latest/) on the model to examine th feature contribution to the classification decision.
'