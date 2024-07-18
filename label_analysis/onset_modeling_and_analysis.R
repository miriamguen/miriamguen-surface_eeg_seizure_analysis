# Load necessary libraries
library(brms)
library(dplyr)
library(data.table)
library(tidyr)
library(posterior)
library(bayestestR)
library(ggplot2)
library(ggdist)
library(emmeans)
library(forcats)
library(cmdstanr)

# Set paths and seed
set.seed(1)
set_cmdstan_path("C:/Users/USER/Documents/.cmdstan/cmdstan-2.35.0")

# Define paths
setwd("D:/Ben Gurion University of Negev Dropbox/Miriam Guen/Miriam Guendelman/Projects/surface_eeg_seizure_analysis")
DATA_PATH <- "example_data/labels/analysis_results/"
MODEL_DATA_PATH <- "example_data/labels/analysis_results/models"
FIGURE_PATH <- "figures/clinical_analysis"

# Load and clean data
onset_labels <- read.csv(file.path(DATA_PATH, "agreed_onset.csv"), header = TRUE, stringsAsFactors = FALSE)
onset_labels <- onset_labels %>%
  select(age, gender, etiology, lobe, patient_id, seizure_length, vigilance, classification, bifurcation_start) %>%
  filter(vigilance != "unclear") %>%
  mutate(
    bifurcation_start = factor(bifurcation_start),
    patient_id = factor(patient_id),
    classification = factor(classification),
    vigilance = factor(vigilance),
    lobe = fct_other(lobe, keep = c('temporal', 'frontal'), other_level = 'Other'),
    gender = factor(gender),
    etiology = fct_other(etiology, keep = c('cryptogenic', 'hippocampal sclerosis', 'cortical dysplasia'), other_level = 'Other')
  )

# Define model formula and family
b_formula <- bifurcation_start ~ age + gender + lobe + etiology + vigilance + seizure_length + classification + 
  (vigilance + seizure_length + classification | patient_id)
b_fam <- categorical("logit")

# Calculate log-odds for priors
bifurcation_proportions <- prop.table(table(onset_labels$bifurcation_start))
log_odds_SN_SubH <- log(bifurcation_proportions["SN/SubH"] / (1 - bifurcation_proportions["SN/SubH"]))
log_odds_SNIC <- log(bifurcation_proportions["SNIC"] / (1 - bifurcation_proportions["SNIC"]))
log_odds_SupH <- log(bifurcation_proportions["SupH"] / (1 - bifurcation_proportions["SupH"]))

# Set bifurcation intercept priors
priors <- c(
  set_prior(paste0("normal(", log_odds_SNIC, ", 1)"), class = "Intercept", dpar = "muSNIC"),
  set_prior(paste0("normal(", log_odds_SupH, ", 1)"), class = "Intercept", dpar = "muSupH")
)

validate_prior(b_formula,
               family = b_fam,
               prior = priors,
               data = onset_labels,)

# Fit the model
mod <- brm(
  b_formula,
  family = b_fam,
  prior = priors,
  data = onset_labels,
  backend = getOption("brms.backend", "cmdstanr"),
  cores = 12,
  warmup = 2000, iter = 6000, chains = 6
)

# Save the model
saveRDS(mod, file.path(MODEL_DATA_PATH, "onset_model.RDS"))

# Load the model for diagnostics
mod <- readRDS(file.path(MODEL_DATA_PATH, "onset_model.RDS"))

# Summarize the model
summary(mod, waic = TRUE)

# Posterior predictive checks
pp_check(mod, type = "bars")

# Define significant margin
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
  gender = unique(onset_labels$gender),
  lobe = unique(onset_labels$lobe),
  etiology = unique(onset_labels$etiology),
  vigilance = unique(onset_labels$vigilance),
  classification = unique(onset_labels$classification),
  patient = "NA"
)

nd$posterior <- rvar(posterior_epred(mod, newdata = nd, re_formula = NA))

# Helper functions for the analysis

#NOTICE - this is not the formal bayes factor - it just serves as a sanity check
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

get_diff_row <- function(data, name, value.1, value.2) {
  data['param'] <- data[name]
  row <- data %>%
    group_by(Class) %>%
    summarise(
      parameter = paste(c(name, value.1, value.2), collapse = "_"),
      prob_change = prob[param == value.1] - prob[param == value.2],
      describe_posterior(prob_change, centrality = c('map', 'median'), test = metrics, null = 0, rope_range = rope_rang_diff, rope_ci = rope_ci),
      bayes_factor(prob_change, 0)
    )
  return(row)
}

get_ratio_row <- function(data, name, value.1, value.2) {
  data['param'] <- data[name]
  row <- data %>%
    group_by(Class) %>%
    summarise(
      parameter = paste(c(name, value.1, value.2), collapse = "_"),
      ratio_change = prob[param == value.1] / prob[param == value.2],
      describe_posterior(ratio_change, centrality = c('map', 'median'), test = metrics, null = 1, rope_range = rope_rang_ratio, rope_ci = rope_ci),
      bayes_factor(ratio_change, 1)
    )
  return(row)
}

# Seizure length analysis
nd_sl <- nd %>%
  mutate(as.data.frame(posterior), posterior = NULL) %>%
  group_by(seizure_length) %>%
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`), SNIC = rvar_mean(SNIC), SupH = rvar_mean(SupH)) %>%
  pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

p <- ggplot(nd_sl, aes(x = seizure_length, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(10))
ggsave(file = file.path(FIGURE_PATH, "seizure_length_onset.svg"), width = 4, height = 4, plot = p)

results <- get_diff_row(nd_sl, 'seizure_length', 0, 60)
results_ratio <- get_ratio_row(nd_sl, 'seizure_length', 0, 60)

# Vigilance analysis
nd_vig <- nd %>%
  mutate(as.data.frame(posterior), posterior = NULL) %>%
  group_by(vigilance) %>%
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`), SNIC = rvar_mean(SNIC), SupH = rvar_mean(SupH)) %>%
  pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

vig_order <- factor(nd_vig$vigilance, levels = c('REM', 'awake', 'sleep stage I', 'sleep stage II', 'sleep stage III/IV'))

p <- ggplot(nd_vig, aes(x = vig_order, ydist = prob, color = Class)) + 
  stat_slabinterval(side = "both", position = position_dodge(0.4)) + 
  labs(x = "Vigilance state", y = "Estimated distribution") + 
  scale_x_discrete(labels = c('awake' = 'Awake', 'sleep stage I' = 'NREM I', 'sleep stage II' = 'NREM II', 'sleep stage III/IV' = 'NREM III')) +            
  scale_color_manual(values = c("SN/SubH" = "blue", "SNIC" = "gold", "SupH" = "darkgreen")) +
  theme_minimal() + 
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(color = "#dcdcdc"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(margin = margin(r = 0, b = 12), size = 11),
    axis.text.y = element_text(margin = margin(t = 0, l = 12), size = 11),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

ggsave(file = file.path(FIGURE_PATH, "vigilance_onset.svg"), width = 6, height = 4, plot = p)

v1 <- c('sleep stage I', 'sleep stage II', 'sleep stage III/IV', 'sleep stage III/IV', 'sleep stage I', 'sleep stage III/IV')
v2 <- c('awake', 'awake', 'awake', 'sleep stage I', 'sleep stage II', 'sleep stage II')

for (i in seq_along(v1)) {
  results <- rbind(results, get_diff_row(nd_vig, 'vigilance', v1[i], v2[i]))
  results_ratio <- rbind(results_ratio, get_ratio_row(nd_vig, 'vigilance', v1[i], v2[i]))
}

# Classification analysis
nd_class <- nd %>%
  mutate(as.data.frame(posterior), posterior = NULL) %>%
  group_by(classification) %>%
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`), SNIC = rvar_mean(SNIC), SupH = rvar_mean(SupH)) %>%
  pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

class_order <- factor(nd_class$classification, levels = c('SP', 'UC', 'CP', 'SG'))

p <- ggplot(nd_class, aes(x = class_order, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4)) +
  labs(x = "Seizure classification", y = "Estimated distribution") + 
  scale_x_discrete(labels = c('SP' = 'FAS', 'UC' = 'UC', 'CP' = 'FIAS', 'SG' = 'FTBCS')) +            
  scale_color_manual(values = c("SN/SubH" = "blue", "SNIC" = "gold", "SupH" = "darkgreen")) +
  theme_minimal() + 
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(color = "#dcdcdc"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(margin = margin(r = 0, b = 12), size = 11),
    axis.text.y = element_text(margin = margin(t = 0, l = 12), size = 11),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

ggsave(file = file.path(FIGURE_PATH, "classification_onset.svg"), width = 6, height = 4, plot = p)

v1 <- c('SP', 'SG', 'UC', 'SG', 'SP', 'SG')
v2 <- c('CP', 'CP', 'CP', 'UC', 'UC', 'SP')

for (i in seq_along(v1)) {
  results <- rbind(results, get_diff_row(nd_class, 'classification', v1[i], v2[i]))
  results_ratio <- rbind(results_ratio, get_ratio_row(nd_class, 'classification', v1[i], v2[i]))
}

# Etiology analysis
nd_etiology <- nd %>%
  mutate(as.data.frame(posterior), posterior = NULL) %>%
  group_by(etiology) %>%
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`), SNIC = rvar_mean(SNIC), SupH = rvar_mean(SupH)) %>%
  pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

p <- ggplot(nd_etiology, aes(x = etiology, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "etiology_onset.svg"), width = 6, height = 4, plot = p)

results <- rbind(results, get_diff_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis'))
results_ratio <- rbind(results_ratio, get_ratio_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis'))

# Age analysis
nd_age <- nd %>%
  mutate(as.data.frame(posterior), posterior = NULL) %>%
  group_by(age) %>%
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`), SNIC = rvar_mean(SNIC), SupH = rvar_mean(SupH)) %>%
  pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

p <- ggplot(nd_age, aes(x = age, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(10)) + 
  labs(x = "Age group", y = "Estimated distribution") + 
  scale_color_manual(values = c("SN/SubH" = "blue", "SNIC" = "gold", "SupH" = "darkgreen")) +
  theme_minimal() + 
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(color = "#dcdcdc"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(margin = margin(r = 0, b = 12), size = 11),
    axis.text.y = element_text(margin = margin(t = 0, l = 12), size = 11),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

ggsave(file = file.path(FIGURE_PATH, "age_onset.svg"), width = 6, height = 4, plot = p)

ages <- list(c(10, 25), c(10, 40), c(25, 40))

for (age_pair in ages) {
  results <- rbind(results, get_diff_row(nd_age, 'age', age_pair[1], age_pair[2]))
  results_ratio <- rbind(results_ratio, get_ratio_row(nd_age, 'age', age_pair[1], age_pair[2]))
}

# Lobe analysis
nd_lobe <- nd %>%
  mutate(as.data.frame(posterior), posterior = NULL) %>%
  group_by(lobe) %>%
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`), SNIC = rvar_mean(SNIC), SupH = rvar_mean(SupH)) %>%
  pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

p <- ggplot(nd_lobe, aes(x = lobe, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "lobe_onset.svg"), width = 4, height = 4, plot = p)

results <- rbind(results, get_diff_row(nd_lobe, 'lobe', 'frontal', 'temporal'))
results_ratio <- rbind(results_ratio, get_ratio_row(nd_lobe, 'lobe', 'frontal', 'temporal'))

# Gender analysis
nd_gender <- nd %>%
  mutate(as.data.frame(posterior), posterior = NULL) %>%
  group_by(gender) %>%
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`), SNIC = rvar_mean(SNIC), SupH = rvar_mean(SupH)) %>%
  pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

p <- ggplot(nd_gender, aes(x = gender, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file = file.path(FIGURE_PATH, "gender_onset.svg"), width = 4, height = 4, plot = p)

results <- rbind(results, get_diff_row(nd_gender, 'gender', 'female', 'male'))
results_ratio <- rbind(results_ratio, get_ratio_row(nd_gender, 'gender', 'female', 'male'))

# Save results
results <- results %>% select(-prob_change) %>% select(-Parameter)
results_ratio <- results_ratio %>% select(-ratio_change) %>% select(-Parameter)

fwrite(results, file.path(DATA_PATH, "clinical_measures", "onset_results_summary.csv"))
fwrite(results_ratio, file.path(DATA_PATH, "clinical_measures", "onset_results_ratio_summary.csv"))

# save.image(file.path(DATA_PATH, "clinical_measures", "onset_data_summary_new.RData"))
