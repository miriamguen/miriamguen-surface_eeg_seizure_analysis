library(brms)
library(dplyr)
library(data.table)

library(posterior)
library(bayestestR)
library(ggplot2)
library(ggdist)
library(emmeans)
library(forcats)

# Load model -------------------------------------------------------

# Load / Clean data -------------------------------------------------------

setwd("D:/Ben Gurion University of Negev Dropbox/Miriam Guen/Miriam Guendelman/Projects/surface_eeg_seizure_analysis")
DATA_PATH = "example_data/labels/analysis_results/"  
FIGURE_PATH = "figures/clinical_analysis"
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

# _________________________________________________________________



summary(offset_labels$bif_end)

# load the binary model for the offset 
mod = readRDS(paste(MODEL_DATA_PATH,"offset_model.RDS", sep="/"))
# diagnostics: Rhat ESS
summary(mod, waic = TRUE)

# pp check
pp_check(mod, type = "bars") # YAY

# inference + plots

# create a posterior draw from the data 
nd <- expand.grid(
  seizure_length = c(0 ,60),
  age = c(10, 25, 40),
  gender =  unique(offset_labels$gender),
  lobe = unique(offset_labels$lobe),
  etiology = unique(offset_labels$etiology),
  vigilance = unique(offset_labels$vigilance),
  classification =  unique(offset_labels$classification),
  patient = "NA"
)

nd$posterior <- rvar(posterior_epred(mod, newdata = nd, re_formula = NA))
### define help functions -----
nd$posterior <- rvar(posterior_epred(mod, newdata = nd, re_formula = NA))
rope_rang_diff = c(-0.05, 0.05) 
rope_rang_ratio = c((1/1.2), 0.2)
metrics_diff = c("p_map", "pd")
metrics_ratio = c("rope","pd", "p_map")
rope_ci = 1


get_diff_row <- function(data, param, value.1, value.2){
 ind1 = which(data[param]==value.1)
 ind2 = which(data[param]==value.2)
 row <- data |>
    summarise(
      parameter = paste(c(param, value.1, value.2), collapse="_"),
      prob_change = prob[ind1] - prob[ind2],
      describe_posterior(centrality='map', prob_change, test = metrics_diff),
    )
  
  return(row)
}
  
p_map_ratio <- function(x, precision = 4^10, method = "kernel") {
  # Density at MAP
  map_x = map_estimate(x, precision = precision,method = method)$MAP_Estimate
  map <- density_at(x, map_x, precision = precision, method = method)
  print(c('map_estimate: ', map))
  # Density at 1
  d_1 <- density_at(x, 1, precision = precision, method = method)
  print(c('d_1 estimate: ',d_1 ))
  if (is.na(d_1)) d_1 <- 0
  
  # Odds
  p <- d_1 / map
  (print(p))
  #class(p) <- c("p_map", class(p))
  return(p) 
}

get_ratio_row <- function(data, param, value.1, value.2){
  ind1 = which(data[param]==value.1)
  ind2 = which(data[param]==value.2)
  row <- data |>
     summarise(
       parameter = paste(c(param, value.1, value.2), collapse="_"),
       ratio_change = prob[ind1] / prob[ind2],
       describe_posterior(centrality='median', ratio_change, test = metrics_ratio, 
                          null = 1, rope_range = rope_rang_ratio, rope_ci=rope_ci),
       p_map_sanity = p_map_ratio(ratio_change))
      
  return(row)
}



## seizure_length  -----

nd_sl <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = posterior) |> 
  group_by(seizure_length) |> 
  summarise(posterior = rvar_mean(posterior)) |> 
  tidyr::pivot_longer(posterior,  values_to = "prob")

nd_sl

p <- ggplot(nd_sl, aes(x = seizure_length, ydist = prob) )+ 
  stat_slabinterval(position = position_dodge(10))

ggsave(file=paste(FIGURE_PATH, "seizure_length_offset.svg", sep="/"), plot=p)


results <- get_diff_row(nd_sl, 'seizure_length', 0, 60)
results_ratio <- get_ratio_row(nd_sl, 'seizure_length', 0, 60)



## vigilance  -----
summary(factor(offset_labels$vigilance))

nd_vig <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = posterior) |> 
  group_by(vigilance) |> 
  summarise(posterior = rvar_mean(posterior)) |> 
  tidyr::pivot_longer(posterior,  values_to = "prob")


nd_vig
vig_order = factor(nd_vig$vigilance,
                   level = c('REM','awake', 'sleep stage I','sleep stage II','sleep stage III/IV' ))


p <- ggplot(nd_vig, aes(x = vig_order, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file=paste(FIGURE_PATH, "vigilance_offset.svg", sep="/"), plot=p)


v1 = c('sleep stage I', 'sleep stage II', 'sleep stage III/IV', 'sleep stage III/IV', 'sleep stage I',  'sleep stage III/IV')
v2 = c('awake',         'awake',          'awake',              'sleep stage I',      'sleep stage II', 'sleep stage II')
  
for (i in 1:length(v1)){
  n = nrow(results)+1
  results[n,] <- get_diff_row(nd_vig, 'vigilance', v1[i], v2[i])
  results_ratio[n,] <- get_ratio_row(nd_vig, 'vigilance', v1[i], v2[i])
}

## classification -----
summary(factor(offset_labels$classification))

nd_class <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = posterior) |> 
  group_by(classification) |> 
  summarise(posterior = rvar_mean(posterior)) |> 
  tidyr::pivot_longer(posterior,  values_to = "prob")


p <-  ggplot(nd_class, aes(x = classification, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4)) #+ theme(panel.background = element_rect(fill = "white", colour = "grey50"))

ggsave(file=paste(FIGURE_PATH,"classification_offset.svg", sep="/"), plot=p)


v1 = c('SP', 'SG', 'UC', 'SG', 'SP', 'SG')
v2 = c('CP', 'CP', 'CP', 'UC', 'UC', 'SP' )

for (i in 1:length(v1)){
  n = nrow(results)+1
  results[n,] <- get_diff_row(nd_class, 'classification', v1[i], v2[i])
  results_ratio[n,] <- get_ratio_row(nd_class, 'classification', v1[i], v2[i])
}



## etiology ----- 


nd_etiology <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = posterior) |> 
  group_by(etiology) |> 
  summarise(posterior = rvar_mean(posterior)) |> 
  tidyr::pivot_longer(posterior, values_to = "prob")

nd_etiology

p <- ggplot(nd_etiology, aes(x = etiology, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file=paste(FIGURE_PATH,  "etiology_offset.svg", sep="/"), plot=p)

summary(factor(offset_labels$etiology))
# only cortical dysplasia and Ts are clearly defined 

n = nrow(results)+1
results[n,] <- get_diff_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis')
results_ratio[n,] <- get_ratio_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis')


## age -----

nd_age <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = posterior) |> 
  group_by(age) |> 
  summarise(posterior = rvar_mean(posterior)) |> 
  tidyr::pivot_longer(posterior,  values_to = "prob")

nd_age

p <- ggplot(nd_age, aes(x = age, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(10))
ggsave(file=paste(FIGURE_PATH,"age_offset.svg", sep="/"), plot=p)

ages = c(c(10,25), c(10, 40), c(25, 40))

for (i in 1:(length(ages)/2)){
  i = i*2
  print (i)
  n = nrow(results)+1
  results[n,] <-  get_diff_row(nd_age, 'age', ages[i-1], ages[i])
  results_ratio[n,] <- get_ratio_row(nd_age, 'age', ages[i-1], ages[i])
}



## lobe  -----

nd_lobe <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = posterior) |> 
  group_by(lobe) |> 
  summarise(posterior = rvar_mean(posterior)) |> 
  tidyr::pivot_longer(posterior, values_to = "prob")

nd_lobe
summary(factor(offset_labels$lobe))
# only frontal and temporal are clearly defined 

p <- ggplot(nd_lobe, aes(x = lobe, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file=paste(FIGURE_PATH,"lobe_offset.svg", sep="/"), plot=p)


n = nrow(results)+1
results[n,] <- get_diff_row(nd_lobe, 'lobe', 'frontal', 'temporal')
results_ratio[n,] <- get_ratio_row(nd_lobe, 'lobe', 'frontal','temporal')


## gender  -----
summary(factor(offset_labels$gender))

nd_gender <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = posterior) |> 
  group_by(gender) |> 
  summarise(posterior = rvar_mean(posterior)) |> 
  tidyr::pivot_longer(posterior, values_to = "prob")

nd_gender

p <- ggplot(nd_gender, aes(x = gender, ydist = prob)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file=paste(FIGURE_PATH, "gender_offset.svg", sep="/"), plot=p)


n = nrow(results)+1
results[n,] <- get_diff_row(nd_gender, 'gender', 'female', 'male')
results_ratio[n,] <- get_ratio_row(nd_gender, 'gender', 'female','male')


## save the data  -----
fwrite(results[,c(1, 4:9)], paste(DATA_PATH, "clinical_measures", "offset_results_summary.csv", sep="/"))
fwrite(results_ratio[,c(1,4:14)], paste(DATA_PATH, "clinical_measures", "offset_results_ratio_summary.csv", sep="/"))

# save.image(paste(DATA_PATH, "clinical_measures","offset_data_summary.RData", sep="/"))

#load(paste(DATA_PATH, "clinical_measures","offset_data_summary.RData", sep="/"))