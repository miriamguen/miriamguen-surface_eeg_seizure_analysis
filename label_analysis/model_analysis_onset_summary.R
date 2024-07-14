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

# Load model -------------------------------------------------------

setwd("D:/Ben Gurion University of Negev Dropbox/Miriam Guen/Miriam Guendelman/Projects/surface_eeg_seizure_analysis")
DATA_PATH = "example_data/labels/analysis_results/"  
FIGURE_PATH = "figures/clinical_analysis"
MODEL_DATA_PATH = "example_data/labels/analysis_results/models"

#load and clean the data
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

#_______________________________________________

summary(onset_labels$bifurcation_start)


mod <- readRDS(paste(MODEL_DATA_PATH, "onset_model.RDS", sep="/"))

# diagnostics: Rhat ESS
summary(mod, waic = TRUE)


# pp check
pp_check(mod, type = "bars") # YAY


# create a posterior draw from the data 
nd <- expand.grid(
  seizure_length = c(0,60),
  age = c(10, 25, 40),
  gender =  unique(onset_labels$gender),
  lobe = unique(onset_labels$lobe),
  etiology = unique(onset_labels$etiology),
  vigilance = unique(onset_labels$vigilance),
  classification =  unique(onset_labels$classification),
  patient = "NA"
)

nd$posterior <- rvar(posterior_epred(mod, newdata = nd, re_formula = NA))

### define help functions -----

rope_rang_diff = c(-0.01, 0.01) 
rope_rang_ratio = c((1/1.2), 0.2)
metrics_diff = c( "rope", "pd", "p_map")
metrics_ratio = c("rope", "pd", "p_map")
rope_ci = 1

get_diff_row <- function(data, name, value.1, value.2){
  data['param'] <- data[name]
  row <- data |>
    group_by(Class) |> 
    summarise(
      parameter = paste(c(name, value.1, value.2), collapse="_"),
      prob_change = prob[param==value.1] - prob[param==value.2],
      describe_posterior(prob_change, centrality='median', test=metrics_diff),
      
    )
 
  return(row)
}

p_map_ratio <- function(x, precision = 4^10, method = "kernel") {
  # Density at MAP
  map_x = map_estimate(x, precision = precision,method = method)$MAP_Estimate
  map <- density_at(x, map_x, precision = precision, method = method)
  print(c('map_estimate: ', map))
  # Density at 0
  d_1 <- density_at(x, 1, precision = precision, method = method)
  print(c('d_1 estimate: ',d_1 ))
  if (is.na(d_1)) d_1 <- 0
  
  # Odds
  p <- d_1 / map
  print(p)

  return(p) 
}

get_ratio_row <- function(data, name, value.1, value.2){
  data['param'] <- data[name]
  row <- data |>
    group_by(Class) |> 
    summarise(
      parameter = paste(c(name, value.1, value.2), collapse="_"),
      ratio_change =  prob[param==value.1] / prob[param==value.2],
      describe_posterior(ratio_change, centrality='map', test=metrics_ratio, 
                         null = 1, rope_range=rope_rang_ratio, rope_ci=rope_ci),
      p_map_sanity = p_map_ratio(ratio_change)
    )
  return(row)
}


## seizure_length -----

nd_sl <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = NULL) |> 
  group_by(seizure_length) |> 
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`),
            SNIC = rvar_mean(SNIC),
            SupH = rvar_mean(SupH)) |> 
  tidyr::pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

nd_sl

p <- ggplot(nd_sl, aes(x = seizure_length, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(10))
ggsave(file=paste(FIGURE_PATH ,"seizure_length_onset.svg", sep="/"), width=4, height=4, plot=p)

results <- get_diff_row(nd_sl, 'seizure_length', 0, 60)
results_ratio <- get_ratio_row(nd_sl, 'seizure_length', 0, 60)


## vigilance 
summary(factor(onset_labels$vigilance))

nd_vig <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = NULL) |> 
  group_by(vigilance) |> 
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`),
            SNIC = rvar_mean(SNIC),
            SupH = rvar_mean(SupH)) |> 
  tidyr::pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

vig_order = factor(nd_vig$vigilance,
                   level = c('REM','awake', 'sleep stage I','sleep stage II','sleep stage III/IV' ))


p <- ggplot(nd_vig, aes(x = vig_order, ydist = prob, color=Class)) + 
  stat_slabinterval(side = "both", position = position_dodge(0.4) ) + 
  labs(x = "Vigilance state", y = "Estimated distribution") + 
  scale_x_discrete(labels = c('awake' = 'Awake',  'sleep stage I' = 'NREM I',      
                              'sleep stage II' = 'NREM II','sleep stage III/IV' = 'NREM III')) +            
  scale_color_manual(values = c("SN/SubH" = "blue", "SNIC" = "gold", "SupH" = "darkgreen")) +
  theme_minimal() +  # Starting with a minimal theme
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # Setting panel background to white
    panel.grid.major = element_line(color = "#dcdcdc"),  # Soft gray for major grid lines
    panel.grid.minor = element_blank(),  # Soft gray, dotted for minor grid lines
    panel.grid.major.x = element_blank(),  # Removing vertical grid lines
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(margin = margin(r = 0 , b=12), size = 11), # Increase space between X-axis tick labels and axis label
    axis.text.y = element_text(margin = margin(t = 0, l=12), size = 11),
    axis.title.x = element_text(size = 14), # Change font size of X-axis label
    axis.title.y = element_text(size = 14)  # Change font size of Y-axis label
  )

ggsave(file=paste(FIGURE_PATH ,"vigilance_onset.svg", sep="/"),  width=6, height=4, plot=p)



p <- ggplot(nd_vig, aes(xdist = prob, y = vig_order,  fill = vig_order)) + 
  theme_bw()+
  facet_wrap(~Class, ncol=1, scales = "free", shrink =TRUE) +
  stat_halfeye() +
  labs(
    title = "Bifurcation Proportion",
    subtitle = "By vigilance state",
    y = "Vigilance State",
    fill = 'Vigilance'
  )

ggsave(file=paste(FIGURE_PATH ,"vigilance_onset_sep.svg", sep="/"), plot=p)


v1 = c('sleep stage I', 'sleep stage II', 'sleep stage III/IV', 'sleep stage III/IV', 'sleep stage I',  'sleep stage III/IV')
v2 = c('awake',         'awake',          'awake',              'sleep stage I',      'sleep stage II', 'sleep stage II')

for (i in 1:length(v1)){
  print(v1[i])
  n = nrow(results)+1
  results[c(n,n+1,n+2),] <- get_diff_row(nd_vig, 'vigilance', v1[i], v2[i])
  results_ratio[c(n,n+1,n+2),] <- get_ratio_row(nd_vig, 'vigilance', v1[i], v2[i])
}


## classification -----

summary(factor(onset_labels$classification))


nd_class <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = NULL) |> 
  group_by(classification) |> 
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`),
            SNIC = rvar_mean(SNIC),
            SupH = rvar_mean(SupH)) |> 
  tidyr::pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

class_order = factor(nd_class$classification,
                     level = c('SP', 'UC', 'CP', 'SG' ))


p <- ggplot(nd_class, aes(x = class_order, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4)) +
  labs(x = "Seizure classification", y = "Estimated distribution") + 
  scale_x_discrete(labels = c('SP'='FAS', 'UC'='UC', 'CP'='FIAS', 'SG'='FTBCS')) +            
  scale_color_manual(values = c("SN/SubH" = "blue", "SNIC" = "gold", "SupH" = "darkgreen")) +
  theme_minimal() +  # Starting with a minimal theme
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # Setting panel background to white
    panel.grid.major = element_line(color = "#dcdcdc"),  # Soft gray for major grid lines
    panel.grid.minor = element_blank(),  # Soft gray, dotted for minor grid lines
    panel.grid.major.x = element_blank(),  # Removing vertical grid lines
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(margin = margin(r = 0 , b=12), size = 11), # Increase space between X-axis tick labels and axis label
    axis.text.y = element_text(margin = margin(t = 0, l=12), size = 11),
    axis.title.x = element_text(size = 14), # Change font size of X-axis label
    axis.title.y = element_text(size = 14)  # Change font size of Y-axis label
  )



ggsave(file=paste(FIGURE_PATH ,"classification_onset.svg", sep="/"), width=6, height=4, plot=p)


p <- ggplot(nd_class, aes(xdist = prob, y = class_order,  color = class_order)) + 
  theme_bw()+
  facet_wrap(~Class, ncol=1, scales = "free", shrink =TRUE) +
  stat_halfeye() +
  labs(
    title = "Bifurcation Proportion",
    subtitle = "By seizure classification",
    y = "Seizure classification",
    fill = 'classification'
  )
ggsave(file=paste(FIGURE_PATH ,"classification_onset.svg_sep.svg", sep="/"), plot=p)

v1 = c('SP', 'SG', 'UC', 'SG', 'SP', 'SG')
v2 = c('CP', 'CP', 'CP', 'UC', 'UC', 'SP' )

for (i in 1:length(v1)){
  n = nrow(results)+1
  results[c(n,n+1,n+2),] <- get_diff_row(nd_class, 'classification', v1[i], v2[i])
  results_ratio[c(n,n+1,n+2),] <- get_ratio_row(nd_class, 'classification', v1[i], v2[i])
}

## etiology ----- 


summary(factor(onset_labels$etiology))


nd_etiology <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = NULL) |> 
  group_by(etiology) |> 
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`),
            SNIC = rvar_mean(SNIC),
            SupH = rvar_mean(SupH)) |> 
  tidyr::pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")


p <- ggplot(nd_etiology, aes(x = etiology, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file=paste(FIGURE_PATH ,"etiology_onset.svg", sep="/" ), width=6, height=4, plot=p)


p <- ggplot(nd_etiology, aes(xdist = prob, y = etiology,  color = etiology)) + 
  theme_bw()+
  facet_wrap(~Class, ncol=1, scales = "free", shrink =TRUE) +
  stat_halfeye() +
  labs(
    title = "Bifurcation Proportion",
    subtitle = "By Etiolog",
    y = "Etiolog",
    fill = 'Etiology'
  )
ggsave(file=paste(FIGURE_PATH ,"Etiology_onset.svg_sep.svg", sep="/"), plot=p)

n = nrow(results)+1
results[c(n,n+1,n+2),] <- get_diff_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis')
results_ratio[c(n,n+1,n+2),] <- get_ratio_row(nd_etiology, 'etiology', 'cortical dysplasia', 'hippocampal sclerosis')

## age -----

nd_age <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = NULL) |> 
  group_by(age) |> 
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`),
            SNIC = rvar_mean(SNIC),
            SupH = rvar_mean(SupH)) |> 
  tidyr::pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")


p <- ggplot(nd_age, aes(x = age, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(10)) + 
  labs(x = "Age group", y = "Estimated distribution") + 
  #scale_x_discrete(labels = c('SP'='FAS', 'UC'='UC', 'CP'='FIAS', 'SG'='FTBCS')) +            
  scale_color_manual(values = c("SN/SubH" = "blue", "SNIC" = "gold", "SupH" = "darkgreen")) +
  theme_minimal() +  # Starting with a minimal theme
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # Setting panel background to white
    panel.grid.major = element_line(color = "#dcdcdc"),  # Soft gray for major grid lines
    panel.grid.minor = element_blank(),  # Soft gray, dotted for minor grid lines
    panel.grid.major.x = element_blank(),  # Removing vertical grid lines
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(margin = margin(r = 0 , b=12), size = 11), # Increase space between X-axis tick labels and axis label
    axis.text.y = element_text(margin = margin(t = 0, l=12), size = 11),
    axis.title.x = element_text(size = 14), # Change font size of X-axis label
    axis.title.y = element_text(size = 14)  # Change font size of Y-axis label
  )

ggsave(file=paste(FIGURE_PATH ,"age_onset.svg", sep="/"), width=6, height=4, plot=p)


ages = c(c(10,25), c(10, 40), c(25, 40))

for (i in 1:(length(ages)/2)){
 i = i*2
 print (i)
 n = nrow(results)+1
 results[c(n,n+1,n+2),] <-  get_diff_row(nd_age, 'age', ages[i-1], ages[i])
 results_ratio[c(n,n+1,n+2),] <- get_ratio_row(nd_age, 'age', ages[i-1], ages[i])
}



## lobe  -----

summary(factor(onset_labels$lobe))

nd_lobe <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = NULL) |> 
  group_by(lobe) |> 
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`),
            SNIC = rvar_mean(SNIC),
            SupH = rvar_mean(SupH)) |> 
  tidyr::pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")


p <- ggplot(nd_lobe, aes(x = lobe, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file=paste(FIGURE_PATH ,"lobe_onset.svg", sep="/"), width=4, height=4, plot=p)

n = nrow(results)+1
results[c(n,n+1,n+2),] <- get_diff_row(nd_lobe, 'lobe', 'frontal', 'temporal')
results_ratio[c(n,n+1,n+2),] <- get_ratio_row(nd_lobe, 'lobe', 'frontal','temporal')



## gender  -----
summary(factor(onset_labels$gender))


nd_gender <- nd |> 
  mutate(as.data.frame(posterior),
         posterior = NULL) |> 
  group_by(gender) |> 
  summarise(`SN/SubH` = rvar_mean(`SN/SubH`),
            SNIC = rvar_mean(SNIC),
            SupH = rvar_mean(SupH)) |> 
  tidyr::pivot_longer(`SN/SubH`:`SupH`, names_to = "Class", values_to = "prob")

nd_gender

p <- ggplot(nd_gender, aes(x = gender, ydist = prob, color = Class)) + 
  stat_slabinterval(position = position_dodge(0.4))
ggsave(file=paste(FIGURE_PATH ,"gender_onset.svg", sep="/"), width=4, height=4, plot=p)

n = nrow(results)+1
results[c(n,n+1,n+2),] <- get_diff_row(nd_gender, 'gender', 'female', 'male')
results_ratio[c(n,n+1,n+2),] <- get_ratio_row(nd_gender, 'gender', 'female','male')

## save the data  -----
fwrite(results[,c(1:2,4:9)], paste(DATA_PATH, "clinical_measures","onset_results_summary.csv", sep="/"))
fwrite(results_ratio[,c(1:2,4:15)], paste(DATA_PATH, "clinical_measures","onset_results_ratio_summary.csv", sep="/"))


# save.image(paste(DATA_PATH, "clinical_measures","onset_data_summary_new.RData", sep="/"))


