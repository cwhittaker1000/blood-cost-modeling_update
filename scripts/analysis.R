# Import packages
# TODO: Look into early stopping
# TODO: Improve function naming
# TODO: Make plot style at top
source("functions.R")
library(tidyverse)
library(furrr)
library(scales)
library(future)
library(future.apply)
library(dplyr)

# =============================================================================
# LOAD IN DATA AND ESTIMATE MEW/ALPHA
# =============================================================================

# Load in data
data <- read_csv("../data/new-hiv-data.csv", show_col_types = FALSE)

our_data <- data %>%
  mutate(hiv_ra = hiv_reads / read_depth)

# Use geometric mean to estimate mew
estimate_mew <- our_data %>%
  summarize(geo_mean = exp(mean(log(hiv_ra), na.rm = TRUE))) %>%
  pull(geo_mean)

# TODO: Double check that I can use geometric mean for the mean esitmate, but arithmetic mean for alpha
# Calculate overdispersion using Coefficient of Variation Squared of relative abundances
cv_squared <- our_data %>%
  summarize(
    mean_ra = mean(hiv_ra, na.rm = TRUE),
    sd_ra = sd(hiv_ra, na.rm = TRUE),
    cv_squared = (sd_ra / mean_ra)^2
  ) %>%
  pull(cv_squared)
# Use CV² directly as overdispersion parameter
alpha <- cv_squared

# =============================================================================
# DEFINE OTHER PARAMETERS FOR SIMULATIONS
# =============================================================================

# Constants for simulations
N <- 3e8 # Approximate US population
shedding_weeks <- 12 # The number of weeks an infected person sheds
r <- 0.0155 # weekly growth rate based on HIV
I <- 100 # Number of initial infections

# Varying parameters for simulations
target_cumulative_incidences <- c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1)
batch_sizes <- c(1200, 6000, 12000, 60000, 120000)
mews <- 10^seq(-2, -8, by = -1)
read_thresholds <- c(100, 1000, 10000)

# Setup for running simulations using multiple cores
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)
param_grid <- expand.grid(
  target_cumulative_incidence = target_cumulative_incidences,
  batch_size = batch_sizes,
  mew = mews,
  read_threshold = read_thresholds,
  stringsAsFactors = FALSE
)

# =============================================================================
# RUN EXPENSIVE SIMULATIONS (uncomment to run)
# =============================================================================

#optimal_sequencing_depths_across_param_sweep <- future_apply(
#  X = param_grid, # each row is one job
#  MARGIN = 1,
#  FUN = function(row) {
#    find_optimal_sequencing_depth_with_binary_search(
#      target_cumulative_incidence = row["target_cumulative_incidence"],
#      target_prob = 0.95,
#      r = r,
#      I1 = I1,
#      P = row["batch_size"],
#      theta = row["read_threshold"],
#      mew = row["mew"],
#      alpha = alpha,
#      shedding_weeks = shedding_weeks,
#      N = N,
#      n_sims = 500
#    )
#  },
#  future.seed = TRUE
#)
#
#results_dir <- "../results"
#if (!dir.exists(results_dir)) {
#  dir.create(results_dir)
#}
#saveRDS(
#  optimal_sequencing_depths_across_param_sweep,
#  file = file.path(results_dir, "lenni_final_optimal_results.rds")
#)
#
#rm(optimal_sequencing_depths_across_param_sweep)

# =============================================================================
# LOAD IN RESULTS FROM SIMULATIONS
# =============================================================================

# TODO: Change output file name
opt_seq_depth_results <- readRDS("../results/lenni_final_optimal_results.rds")
# Extract results from each parameter sweep
opt_seq_depth_summary <- param_grid %>%
  mutate(
    optimal_depth = map_dbl(opt_seq_depth_results, "optimal_depth"),
    achieved_prob = map_dbl(opt_seq_depth_results, "final_prob"),
    q25_incidence = map_dbl(
      optimal_results,
      c("final_results", "summary_stats", "q25_incidence")
    ),
    q75_incidence = map_dbl(
      optimal_results,
      c("final_results", "summary_stats", "q75_incidence")
    ),
    annual_reads = optimal_depth * 52
  )
# Turn to a tibble for easier manipulation
opt_seq_depth_summary <- as_tibble(opt_seq_depth_summary)

# =============================================================================
# PLOTS
# =============================================================================

# Extract data for a single mew closest to our estimated mew, and remove any
# parameter combinations that did not converge (hit our ceiling)
select_mew_read_threshold_data <- opt_seq_depth_summary %>%
  filter(mew == 1e-6, read_threshold == 100) %>%
  mutate(
    # Flag values that hit the ceiling
    hit_ceiling = optimal_depth >= 1e12
  ) %>%
  filter(!hit_ceiling)
depth_plot <- ggplot(
  select_mew_read_threshold_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1, linetype = "solid") +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(NA, 1e12)
  ) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  labs(
    title = "",
    y = "Sequencing depth per week to detect HIV-like pathogen",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
depth_plot

ggsave("../output/depth_plot.png", depth_plot, width = 8, height = 6)


## Cost analysis for a specific depth
seq_c <- 2500 / 1e9
proc_c <- 2500 #TODO: Update to new value that we decide on
sample_c <- 15

# Add the cost estimate column
select_mew_read_threshold_data_with_cost <- select_mew_read_threshold_data %>%
  filter(!hit_ceiling) %>% # Exclude points that didn't converge
  mutate(
    total_cost_per_year = calculate_total_cost_per_year(
      seq_depth = optimal_depth,
      cost_of_seq = seq_c,
      cost_of_proc = proc_c,
      cost_of_sample = sample_c,
      num_samples = batch_size
    )
  )
cost_plot <- ggplot(
  select_mew_read_threshold_data_with_cost,
  aes(
    x = target_cumulative_incidence,
    y = total_cost_per_year,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(
    breaks = c(10000000, 100000000),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  labs(
    title = "",
    y = "Total cost per year to detect HIV-like pathogen ($)",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
cost_plot

ggsave("../output/cost_plot.png", cost_plot, width = 8, height = 6)

# Only filter by read threshold so that we can show a plot
# comparing different mews
select_read_threshold_data <- opt_seq_depth_summary %>%
  filter(read_threshold == 100) %>%
  mutate(
    # Flag values that hit the ceiling
    hit_ceiling = optimal_depth >= 1e12
  ) %>%
  filter(!hit_ceiling)

facet_depth_plot <- ggplot(
  select_read_threshold_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_point(
    size = 3,
  ) +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(NA, 1e12)
  ) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  labs(
    title = "",
    y = "Sequencing depth per week to detect HIV-like pathogen",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 15) +
  facet_grid(. ~ mew) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
facet_depth_plot

ggsave(
  "../output/facet_mew_depth_plot.png",
  facet_depth_plot,
  width = 25,
  height = 8
)


## Cost analysis for multiple sample costs
sample_costs <- c(3, 12, 50, 100)

# Add the cost estimate column
select_mew_read_threshold_data_with_varied_costs <- select_read_threshold_data %>%
  filter(!hit_ceiling) %>% # Exclude points that didn't converge
  crossing(sample_cost = sample_costs) %>% # Create all combinations
  mutate(
    total_cost_per_year = find_total_cost(
      seq_depth = optimal_depth,
      cost_of_seq = seq_c,
      cost_of_proc = proc_c,
      cost_of_sample = sample_cost,
      num_samples = batch_size
    ) *
      52
  )

# Compare costs for different mews with different sample costs
facet_cost_plot <- ggplot(
  select_mew_read_threshold_data_with_varied_costs,
  aes(
    x = target_cumulative_incidence,
    y = total_cost_per_year,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  #  facet_wrap(~ paste0("Sample cost: $", sample_cost), scales = "free_y") +
  facet_grid(sample_cost ~ mew) +
  labs(
    title = "Annual cost with varying sample collection costs",
    y = "Total cost per year to detect HIV-like pathogen ($)",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
facet_cost_plot

ggsave(
  "../output/facet_cost_plot.png",
  facet_cost_plot,
  width = 25,
  height = 15
)
