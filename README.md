# The following repository is comprised of the following code sections, intended to be used in this order:

## The analysis overview:
The following figure presents the analysis flow in this project for more detail see ["link to publication"]
![image](https://github.com/miriamguen/miriamguen-surface_eeg_seizure_analysis/blob/main/figures/paper_figures/Figure%201%20-%20The%20project%20work%20flow.png)



## 1. preprocess_data based on Matlab 2022a with EEGLAB 2022
   ### Run the main.m
   * loads the EEG data using the index tables and seizure onset/offset times in the format provided in the example data
   * Applies the automated preprocessing pipeline
   * saves a file containing the seizure signals and a control segment from the same file with a similar length.

   This project is based on data from the EPILEPSIAE surface database described in this [publication](https://pubmed.ncbi.nlm.nih.gov/20863589/). 

## 2. visual_analysis
   ### Run the label_all.m
   This file sequentially loads the files saved in the previous step and provides a visual graphical user interface to label the data.
   Before using this step on new data, it is recommended that the rater practices and understands the ICA representations provided
   A good tutorial and practice option can be found at https://labeling.ucsd.edu/tutorial/practice

   if the user is interested in labeling bifurcation morphology, it is recommended to practice on the simulation data
   this is provided with the paper ["A taxonomy of seizure dynamotypes"](https://doi.org/10.7554%2FeLife.55632) and can be downloaded [here](https://doi.org/10.7302/ejhy-5h41)
   To run this code, install the 'viewprops' update for EEGLAB from this ["link"](https://github.com/sccn/viewprops/tree/9db7a1119a1d3da1ac0847f3ce3026842843e8fa).
   You can adjust the opening screen sizes by modifying lines 56 and 62.

   ### The information presented includes:
   
   ![image](https://github.com/miriamguen/miriamguen-surface_eeg_seizure_analysis/blob/main/figures/paper_figures/Figure%202%20-%20Manual%20labeling%20interface.png)



## 3. label_analysis
   ### Run the label_validation_analysis.py (tested up to Python 3.12)
   * Analyze detectability rates based on the lobe containing the epileptogenic zone and seizure classification.
   * Run the inter-rater agreement analysis
   * Produce the label distribution plots and tables
   * Compare the labeler selection to automated measures extracted from the signal
   * Run the selection of a single component per seizure based on: TODO: add flow chart

![image](https://github.com/miriamguen/miriamguen-surface_eeg_seizure_analysis/blob/main/figures/paper_figures/Figure%203%20-%20label%20proportions.png)

   ### Run the patient_metadata.py
   * get a printout of the clinical factors of the patients in the data

   ### Bayesian modeling of clinical factors, use: (tested up to R 4.4 + [brms 2.21.0](https://paul-buerkner.github.io/brms/) + cmdstanr 2.35.0) 
   * "onset_modeling_and_analysis.R" to fit the model and run inference analysis with the onset labels
   * "offset_modeling_and_analysis.R to fit the model and run inference analysis with the offset labels

These scripts will save the model, inference result tables, and plots in the "example_data" folder.


## 4. Automated classification (tested up to Python 3.12)
   ### Includes a file to run for each classification task
   * Identify components with a clear transition
   * Identify slowing of the inter-spike intervals towards seizure edges
   
   ### Each file will:
   * Train and evaluate the model in a leave-one-out cross-validation approach
   * Save performance metrics, including balanced accuracy, ROC-AUC, and Cohen's Kappa
   * Train a full model for feature contribution analysis and run [SHAP analysis](https://shap.readthedocs.io/en/latest/) on the model to examine the feature contribution to the classification decision.
'
![image](https://github.com/miriamguen/miriamguen-surface_eeg_seizure_analysis/blob/main/figures/paper_figures/Figure%204%20-%20classification%20analysis.png)