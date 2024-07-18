# Load necessary libraries
library(brms)
library(dplyr)
library(data.table)
library(forcats)
library(cmdstanr)
library(posterior)
library(bayestestR)
library(ggplot2)
library(ggdist)
library(emmeans)

# Set paths and cmdstan path
set_cmdstan_path("C:/Users/USER/Documents/.cmdstan/cmdstan-2.35.0")

# Define paths
setwd("D:/Ben Gurion University of Negev Dropbox/Miriam Guen/Miriam Guendelman/Projects/surface_eeg_seizure_analysis")
DATA_PATH <- "example_data/labels/analysis_results/"
MODEL_DATA_PATH <- "example_data/labels/analysis_results/models"
FIGURE_PATH <- "figures/clinical_analysis"

# Load and clean data
offset_labels <- read.csv(file.path(DATA_PATH, 'agreed_offset.csv'), header = TRUE, stringsAsFactors = FALSE)
offset_labels <- offset_labels %>%
  select(age, gender, lobe, etiology, patient_id, seizure_length, vigilance, classification, bifurcation_end) %>%
  filter(vigilance != "unclear", bifurcation_end != 'SupH') %>%
  mutate(
    patient_id = factor(patient_id),
    bifurcation_end = factor(bifurcation_end),
    classification = factor(classification),
    vigilance = factor(vigilance),
    lobe = fct_other(lobe, keep = c('temporal', 'frontal'), other_level = 'Other'),
    lobe = factor(lobe),
    gender = factor(gender),
    etiology = fct_other(etiology, keep = c('cryptogenic', 'hippocampal sclerosis', 'cortical dysplasia'), other_level = 'Other'),
    etiology = factor(etiology)
  )

# Define model formula and family
b_formula <- bifurcation_end ~ age + gender + lobe + etiology + vigilance + seizure_length + classification + (vigilance + seizure_length + classification | patient_id)
b_fam_bi <- bernoulli() 

# Calculate log-odds for priors
bifurcation_proportions <- prop.table(table(offset_labels$bifurcation_end))
log_odds_FLC <- log(bifurcation_proportions["FLC"] / (1 - bifurcation_proportions["FLC"]))

# Set prior
prior_intercept <- set_prior(paste0("normal(", log_odds_FLC, ", 1)"), class = "Intercept")


validate_prior(b_formula,
               prior = prior_intercept,
               family = b_fam_bi,
               data = offset_labels)

# Fit the model
mod <- brm(
  b_formula,
  prior = prior_intercept,
  family = b_fam_bi,
  data = offset_labels,
  backend = getOption("brms.backend", "cmdstanr"),
  cores = 12,
  warmup = 2000, iter = 6000, chains = 6
)

# Save the model
saveRDS(mod, file.path(MODEL_DATA_PATH, "offset_model.RDS"))

# Load the model for diagnostics
mod <- readRDS(file.path(MODEL_DATA_PATH, "offset_model.RDS"))

# Summarize the model
summary(mod, waic = TRUE)

# Posterior predictive checks
pp_check(mod, type = "bars")

# Define significant margin and statistical analysis parameters
significant_margin <- 0.1
significant_abs_change <- min(bifurcation_proportions) * significant_margin
rope_rang_diff <- c(-significant_abs_change, significant_abs_change)
rope_rang_ratio <- c(1 / (1 + significant_margin), 1 + significant_margin)
metrics <- c("ci", "pd", "rope", "p_map", "p_significant")
rope_ci <- 1

# Create posterior draws from the data
nd <- expand.grid(
  seizure_length = c(0, 60),
  age = c(10, 25, 40),
  gender = unique(offset_labels$gender),
  lobe = unique(offset_labels$lobe),
  etiology = unique(offset_labels$etiology),
  vigilance = unique(offset_labels$vigilance),
  classification = unique(offset_labels$classification),
  patient = "NA"
)

nd$posterior <- rvar(posterior_epred(mod, newdata = nd, re_formula = NA))

# Helper functions for Bayesian analysis
bayes_factor <- function(rvar_diff, null) {
  set.seed(1)
  diff_samples <- draws_of(rvar_diff)
  density_diff <- density(diff_samples)
  posterior_density_at_null <- approx(density_diff$x, density_diff$y, xout = null)$y
  
  null_density <- density(rnorm(length(diff_samples), mean = null, sd = sd(diff_samples)))
  prior_density_at_null <- approx(null_density$x, null_density$y, xout = null)$y
  
  bf <- prior_density_at_null / posterior_density_at_null
  return(bf)
}

get_diff_row <- function(data, param, value.1, value.2) {
  ind1 <- which(data[[param]] == value.1)
  ind2 <- which(data[[param]] == value.2)
  row <- data %>%
    summarise(
      parameter = paste(c(param, value.1, value.2), collapse = "_"),
      prob_change = prob[ind1] - prob[ind2],
      describe_posterior(prob_change, centrality = c('map', 'median'), test = metrics, null = 0, rope_range = rope_rang_diff, rope_ci = rope_ci),
      bayes_factor(prob_change, 0)
    )
  return(row)
}

get_ratio_row <- function(data, param, value.1, value.2) {
  ind1 <- which(data[[param]] == value.1)
  ind2 <- which(data[[param]] == value.2)
  row <- data %>%
    summarise(
      parameter = paste(c(param, value.1, value.2), collapse = "_"),
      ratio_change = prob[ind1] / prob[ind2],
      describe_posterior(ratio_change, centrality = c('map', 'median'), test = metrics, null = 1, rope_range = rope_rang_ratio, rope_ci = rope_ci),
      bayes_factor(ratio_change, 1)
    )
  return(row)
}


# Analysis by seizure length
nd_sl <- nd %>%
  mutate(as.data.frame(posterior)) %>%
  group_by(seizure_length) %>%
  summarise(posterior = rvar_mean(posterior)) %>%
  tidyr::pivot_longer(posterior, values_to = "prob")

p <- ggplot(nd_sl, aes(x = seizure_length, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(10))
ggsave(file = file.path(FIGURE_PATH, "seizure_length_offset.svg"), plot = p)

results <- get_diff_row(nd_sl, 'seizure_length', 0, 60)
results_ratio <- get_ratio_row(nd_sl, 'seizure_length', 0, 60)

# Analysis by vigilance
nd_vig <- nd %>%
  mutate(as.data.frame(posterior)) %>%
  group_by(vigilance) %>%
  summarise(posterior = rvar_mean(posterior)) %>%
  tidyr::pivot_longer(posterior, values_to = "prob")

vig_order <- factor(nd_vig$vigilance, levels = c('REM', 'awake', 'sleep stage I', 'sleep stage II', 'sleep stage III/IV'))

p <- ggplot(nd_vig, aes(x = vig_order, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "vigilance_offset.svg"), plot = p)

v1 <- c('sleep stage I', 'sleep stage II', 'sleep stage III/IV', 'sleep stage III/IV', 'sleep stage I', 'sleep stage III/IV')
v2 <- c('awake', 'awake', 'awake', 'sleep stage I', 'sleep stage II', 'sleep stage II')

for (i in seq_along(v1)) {
  results <- rbind(results, get_diff_row(nd_vig, 'vigilance', v1[i], v2[i]))
  results_ratio <- rbind(results_ratio, get_ratio_row(nd_vig, 'vigilance', v1[i], v2[i]))
}

# Analysis by classification
nd_class <- nd %>%
  mutate(as.data.frame(posterior)) %>%
  group_by(classification) %>%
  summarise(posterior = rvar_mean(posterior)) %>%
  tidyr::pivot_longer(posterior, values_to = "prob")

p <- ggplot(nd_class, aes(x = classification, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "classification_offset.svg"), plot = p)

v1 <- c('SP', 'SG', 'UC', 'SG', 'SP', 'SG')
v2 <- c('CP', 'CP', 'CP', 'UC', 'UC', 'SP')

for (i in seq_along(v1)) {
  results <- rbind(results, get_diff_row(nd_class, 'classification', v1[i], v2[i]))
  results_ratio <- rbind(results_ratio, get_ratio_row(nd_class, 'classification', v1[i], v2[i]))
}

# Analysis by etiology
nd_etiology <- nd %>%
  mutate(as.data.frame(posterior)) %>%
  group_by(etiology) %>%
  summarise(posterior = rvar_mean(posterior)) %>%
  tidyr::pivot_longer(posterior, values_to = "prob")

p <- ggplot(nd_etiology, aes(x = etiology, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "etiology_offset.svg"), plot = p)

results <- rbind(results, get_diff_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis'))
results_ratio <- rbind(results_ratio, get_ratio_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis'))

# Analysis by age
nd_age <- nd %>%
  mutate(as.data.frame(posterior)) %>%
  group_by(age) %>%
  summarise(posterior = rvar_mean(posterior)) %>%
  tidyr::pivot_longer(posterior, values_to = "prob")

p <- ggplot(nd_age, aes(x = age, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(10))
ggsave(file = file.path(FIGURE_PATH, "age_offset.svg"), plot = p)

ages <- list(c(10, 25), c(10, 40), c(25, 40))

for (age_pair in ages) {
  results <- rbind(results, get_diff_row(nd_age, 'age', age_pair[1], age_pair[2]))
  results_ratio <- rbind(results_ratio, get_ratio_row(nd_age, 'age', age_pair[1], age_pair[2]))
}

# Analysis by lobe
nd_lobe <- nd %>%
  mutate(as.data.frame(posterior)) %>%
  group_by(lobe) %>%
  summarise(posterior = rvar_mean(posterior)) %>%
  tidyr::pivot_longer(posterior, values_to = "prob")

p <- ggplot(nd_lobe, aes(x = lobe, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "lobe_offset.svg"), plot = p)

results <- rbind(results, get_diff_row(nd_lobe, 'lobe', 'frontal', 'temporal'))
results_ratio <- rbind(results_ratio, get_ratio_row(nd_lobe, 'lobe', 'frontal', 'temporal'))

# Analysis by gender
nd_gender <- nd %>%
  mutate(as.data.frame(posterior)) %>%
  group_by(gender) %>%
  summarise(posterior = rvar_mean(posterior)) %>%
  tidyr::pivot_longer(posterior, values_to = "prob")

p <- ggplot(nd_gender, aes(x = gender, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "gender_offset.svg"), plot = p)

results <- rbind(results, get_diff_row(nd_gender, 'gender', 'female', 'male'))
results_ratio <- rbind(results_ratio, get_ratio_row(nd_gender, 'gender', 'female', 'male'))

# Save results
results <- results %>% select(-prob_change) %>% select(-Parameter)
results_ratio <- results_ratio %>% select(-ratio_change) %>% select(-Parameter)

fwrite(results, file.path(DATA_PATH, "clinical_measures", "offset_results_summary.csv"))
fwrite(results_ratio, file.path(DATA_PATH, "clinical_measures", "offset_results_ratio_summary.csv"))
# save.image(file.path(DATA_PATH, "clinical_measures", "offset_data_summary.RData"))
