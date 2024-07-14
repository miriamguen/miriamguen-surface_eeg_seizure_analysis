# use cmdstanr backend
library(brms)
library(dplyr)
library(data.table)
library(forcats)
library(cmdstanr)
set_cmdstan_path("C:/Users/USER/Documents/.cmdstan/cmdstan-2.35.0")

# Load / Clean data -------------------------------------------------------
setwd("D:/Ben Gurion University of Negev Dropbox/Miriam Guen/Miriam Guendelman/Projects/surface_eeg_seizure_analysis")
DATA_PATH = "example_data/labels/analysis_results/"  
MODEL_DATA_PATH = "example_data/labels/analysis_results/models"

onset_labels <- read.csv(paste(DATA_PATH, "agreed_onset.csv", sep="/"), header=TRUE, stringsAsFactors=FALSE)
onset_labels <- onset_labels[c("age", "gender" , "etiology", "lobe",  
                               "patient_id",
                               "seizure_length", "vigilance","classification", 
                               "bifurcation_start")]


onset_labels = onset_labels[onset_labels["vigilance"] != "unclear",] #  seizures without clear vigilance [n=6]

# make sure that all categorical features are defined as factors
onset_labels$bifurcation_start <- factor(onset_labels$bifurcation_start)
onset_labels$patient_id <- factor(onset_labels$patient_id)

onset_labels$classification <- factor(onset_labels$classification)
onset_labels$vigilance <- factor(onset_labels$vigilance)

onset_labels$lobe <- fct_other(onset_labels$lobe,
                                   keep = c('temporal', 'frontal'),
                                   other_level = 'Other')

onset_labels$lobe <- factor(onset_labels$lobe)

onset_labels$gender <- factor(onset_labels$gender)

onset_labels$etiology <- fct_other(onset_labels$etiology,
                                   keep = c('cryptogenic', 'hippocampal sclerosis', 'cortical dysplasia'),
                                   other_level = 'Other')

onset_labels$etiology <- factor(onset_labels$etiology)



# Priors and assumptions -----------------------------------------------------------------
# flat priors 

# We are ignoring hospital because there are too few levels to be able to say
# something meaningful

b_formula <- 
  bifurcation_start ~ 
  # Subject level
  age + gender + lobe + etiology + 
  # Within Subject
  vigilance + seizure_length + classification + 
  (vigilance + seizure_length + classification | patient_id)

b_fam <- categorical("logit")


mod <- brm(
  b_formula,
  family = b_fam,
  prior = get_prior(b_formula, family = b_fam, data = onset_labels),
  data = onset_labels,
  backend = getOption("brms.backend", "cmdstanr"),
  cores = 12,
  # threads = 2,
  warmup = 2000, iter = 6000, chains = 6,
)

saveRDS(mod, paste(MODEL_DATA_PATH, "onset_model.RDS", sep="/"))

# assert the model loads and test fit 
mod = readRDS( paste(MODEL_DATA_PATH, "onset_model.RDS", sep="/"))

# diagnostics: Rhat ESS
summary(mod, waic = TRUE)

# pp check
pp_check(mod, type = "bars") # YAY

