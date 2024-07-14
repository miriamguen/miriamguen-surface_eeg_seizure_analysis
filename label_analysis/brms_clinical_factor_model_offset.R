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

offset_labels <- read.csv(paste(DATA_PATH, 'agreed_offset.csv', sep="/"), header=TRUE, stringsAsFactors=FALSE)

offset_labels <- offset_labels[c("age", "gender" , "lobe", "etiology",  
                                 "patient_id",
                                 "seizure_length", "vigilance","classification", 
                                 "bifurcation_end")]

offset_labels = offset_labels[offset_labels["vigilance"] != "unclear",] # drop segments without clear vigilance [n=2]
offset_labels = offset_labels[offset_labels$bifurcation_end != 'SupH', ]# drop the SupH seizures [n=2]
# 467 obs remaining

#make sure the relevant features are set as factors 
offset_labels$patient_id <- factor(offset_labels$patient_id)
offset_labels$bifurcation_end <- factor(offset_labels$bifurcation_end)

offset_labels$classification <- factor(offset_labels$classification)
offset_labels$vigilance <- factor(offset_labels$vigilance)

offset_labels$lobe <- fct_other(offset_labels$lobe,
                                keep = c('temporal', 'frontal'),
                                other_level = 'Other')

offset_labels$lobe <- factor(offset_labels$lobe)
offset_labels$gender <- factor(offset_labels$gender)

offset_labels$etiology <- fct_other(offset_labels$etiology,
                                   keep = c('cryptogenic', 'hippocampal sclerosis', 'cortical dysplasia'),
                                   other_level = 'Other')

offset_labels$etiology <- factor(offset_labels$etiology)



# two class model -----------------------------------------------------------------

b_formula <- 
  bifurcation_end ~ 
  # Subject level
  age + gender + lobe + etiology + 
  # Within Subject
  vigilance + seizure_length + classification + 
  (vigilance + seizure_length + classification | patient_id)

b_fam_bi <-  bernoulli() # for 2 levels


mod <- brm(
  b_formula,
  prior = get_prior(b_formula, family = b_fam_bi, data = offset_labels),
  family = b_fam_bi,
  data = offset_labels,
  backend = getOption("brms.backend", "cmdstanr"),
  cores = 12,
  warmup = 2000, iter = 6000, chains = 6,
)

saveRDS(mod, paste(MODEL_DATA_PATH, "offset_model.RDS", sep="/"))


# test the saved model loads 
mod = readRDS( paste(MODEL_DATA_PATH ,  "offset_model.RDS", sep="/"))


# diagnostics: Rhat ESS
summary(mod, waic = TRUE)


# pp check
pp_check(mod, type = "bars") # YAY
