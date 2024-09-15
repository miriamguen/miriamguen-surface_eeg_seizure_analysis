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


(A)	Displays a time series of the independent components (ICs) of an example seizure. The ICs highlighted in dark blue contain a clear transition into a seizure that includes data mainly from brain sources. The PCA-ICA procedure affects the components' scale; hence, the units are presented as procedure-determined units (PDU).
 
(B)	For the selection of the ICs, we used custom code and EEGLAB visualization functions to present each reviewer with the component time series, the topographic distribution of the electrode weights, the spectral profile, and the IClabel classifier output. The topographic maps represent the normalized weight (NW) distribution used to transform the EEG time series into the IC time series. These maps help evaluate the scalp distribution of each component, helping with the recognition of patterns typical for brain-related sources (IC1, 3, 5), eye movement (IC2, 6), or other non-brain sources (IC4, 7, 8) (Pion-Tonachini et al. 2019b). Using the power spectrum, we could better distinguish ICs with typical brain frequency profiles (IC1, 3, 5) from components containing mainly muscular activity (IC4, 6, 8). These and additional signal features are input to the IClabel classifier (Pion-Tonachini et al. 2019b). This classifier is built into EEGLAB and is designed to identify and remove noise components in surface EEG preprocessing (Pion-Tonachini et al. 2019a). To support our decision and estimate the consistency of this classifier on seizure data, we presented the classification probabilities to the user. Extended Figure 2A displays the overall brain score distribution given the different labels and shows that clear transitions are generally associated with a higher brain score. However, exceptions exist, and this does not provide a clear-cut role to identify seizure data.
 
(C)	Examples of the three onset bifurcation morphologies (BMs) are as follows: (a) SupH, which shows an amplitude increase (from zero); (b) SNIC, characterized by increased spike frequency (from zero) or decreasing inter-spike intervals (ISIs) from the seizure onset; (c) SN (expected DC shift) or SubH, in which both amplitude and ISI remain constant or random; due to the lack of DC recording, these are indistinguishable. Examples of the three offset BMs include (d) SNIC or SH (expected DC shift), characterized by increasing ISIs or reduction in frequency (to zero) and are indistinguishable like the previous example; (e) SupH, with an amplitude decrease (to zero); and (f) FLC, with amplitude and ISI that are constant or random towards seizure offset. See Extended Figure 1 for additional examples.
 
(D)	An increasing amplitude at seizure onset characterizing the SupH BM may be attributed to propagation in non-invasive EEG. Therefore, it is important to analyze the component time series in the context of additional information beyond the time series alone. For example, if we observe only IC2, it shows a SupH BM starting at the 20-second mark. However, when examining IC3, we see a lower constant amplitude beginning at the 16-second mark and undergoing a second change at the 20-second mark. The topographic maps of these components reveal a more localized weight distribution in IC3 than in IC2. This suggests that the pattern in IC2 may result from propagation, emphasizing the importance of selecting the earliest and most localized spatial pattern (in focal seizures) as the seizure onset BM.  


## 3. label_analysis
   ### Run the label_validation_analysis.py (tested up to Python 3.12)
   * Analyze detectability rates based on the lobe containing the epileptogenic zone and seizure classification.
   * Run the inter-rater agreement analysis
   * Produce the label distribution plots and tables
   * Compare the labeler selection to automated measures extracted from the signal
   * Run the selection of a single component per seizure based on: TODO: add flow chart

![image](https://github.com/miriamguen/miriamguen-surface_eeg_seizure_analysis/blob/main/figures/paper_figures/Figure%203%20-%20label%20proportions.png)
 
The distribution of BM in our data, encompassing all seizures marked by onset or offset bifurcations, aligns broadly with previously reported intracranial data. However, the proportions of SN/SubH and FLC morphologies are higher in surface data. These morphologies do not necessitate specific patterns, such as ISI or amplitude changes. Additionally, the higher noise levels in surface data may obscure other morphologies, potentially leading to overestimating these patterns.

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
### General Classifier Analytics
![image](https://github.com/miriamguen/miriamguen-surface_eeg_seizure_analysis/blob/main/figures/paper_figures/Figure%204%20-%20classification%20analysis.png)


This figure presents the analysis of two classification tasks. The first classifier identifies components with a clear transition among all given components, while the second classifier identifies components exhibiting inter-spike interval (ISI) slowing based on bifurcation morphology (BM). Both classifiers were trained using the same feature set and algorithm hyperparameters, with class weights adjusted differently to account for class imbalance and achieve balanced performance.
Each classifier was evaluated using a leave-one-subject-out cross-validation approach. The confusion matrices (panels A and C) display cumulative classification results on the test sets. To understand the classification process, feature contributions were analyzed using SHAP analysis (Rodríguez-Pérez and Bajorath 2020), highlighting the most influential features of each classifier and their directional effects.
 
(A)	The confusion matrix for the transition detection task.
 
(B)	Dominant features for the transition detection task. High brain probability scores were associated with clear components, while high probability scores for undefined noise sources ("other") indicated lower likelihoods of being seizure components. High spectral entropy, a feature previously used in seizure detection tasks (Boonyakitanont et al. 2020), supported positive decisions. High mean peak amplitude in the first 5 seconds of detected onset/offset was related to clear seizure transitions.

(C)	The confusion matrix for the ISI slowing classification task.

(D)	Top features for the ISI slowing classification task. A value indicating an offset increased the probability of ISI slowing due to the higher likelihood of this pattern occurring at offsets. Seizure length was positively correlated with ISI probability, although no significant trend was observed in multilevel analysis. Lower dominant frequency and higher standard deviation in ISI were associated with ISI slowing. The ISI linear intercept estimate (higher CI limit) and lower fit quality (RMSE) to the ISI trend curves were positively related to slowing. This may be due to low error in cases where the fitted line is flat and the observation that slowing usually occurs in components with higher amplitudes, as suggested by the peak trend estimation relation.

