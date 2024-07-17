# Install necessary packages if not already installed
if (!requireNamespace("see", quietly = TRUE)) install.packages("see")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

# Load necessary libraries
library(ggplot2)
library(ggdist)
library(bayestestR)
library(dplyr)
library(see)
library(patchwork)


setwd("D:/Ben Gurion University of Negev Dropbox/Miriam Guen/Miriam Guendelman/Projects/surface_eeg_seizure_analysis")
# Generate distributions
posterior <- distribution_normal(10000, 0.5, 0.2)
prior <- distribution_normal(10000, mean = 0, sd = 0.2)


# Create each plot
ci_plot <- plot(ci(posterior, ci = 0.95)) +
  labs(title = "95% Credible interval") +  
  theme_lucid() + 
  scale_fill_see() + 
  scale_colour_see()

pd_plot <- plot(p_direction(posterior)) +
  labs(title = "Probability of Direction") + 
  theme_lucid() + 
  scale_fill_see() + 
  scale_colour_see() 

rope_plot <- plot(rope(posterior, range = c(-0.1, 0.1))) +
  labs(title = "ROPE") +
  theme_lucid()+ scale_fill_see() + 
  scale_colour_see() + scale_fill_see_d()

bf_plot <- plot(bayesfactor_parameters(posterior, prior, direction = "two-sided", null = 0, verbose = FALSE)) +
  labs(title = "Bayes Factor") +
  theme_lucid()+ scale_fill_see() + 
  scale_colour_see()

# Combine plots into one figure
combined_plot <- ci_plot + pd_plot + rope_plot + bf_plot + plot_layout(ncol = 4)

# Print the combined plot
print(combined_plot)

# Save the combined plot
ggsave("bayesian_measures_combined.svg", plot = combined_plot, width = 16, height = 4)
