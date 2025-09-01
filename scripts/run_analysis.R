# Import packages
# TODO: Look into early stopping of binary search if we don't think we'll converge.
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

# TODO: Question for Dan: Can I use geometric mean for the mean estimate, but arithmetic mean for alpha? ChatGPT says it's wrong as the NB variance will then be biased; however, I know we talk about using geometric mean since there are few samples (~3), so then how would you suggest I estimate the overdispersion parameter.
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
mews <- c(
  estimate_mew * -1e3,
  estimate_mew * -1e2,
  estimate_mew * -1e1,
  estimate_mew,
  estimate_mew * 1e1,
  estimate_mew * 1e2,
  estimate_mew * 1e3
)
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

optimal_sequencing_depths_across_param_sweep <- future_apply(
  X = param_grid, # each row is one job
  MARGIN = 1,
  FUN = function(row) {
    find_optimal_sequencing_depth_with_binary_search(
      target_cumulative_incidence = row["target_cumulative_incidence"],
      target_prob = 0.95,
      r = r,
      I1 = I1,
      P = row["batch_size"],
      theta = row["read_threshold"],
      mew = row["mew"],
      alpha = alpha,
      shedding_weeks = shedding_weeks,
      N = N,
      n_sims = 500
    )
  },
  future.seed = TRUE
)
results_dir <- "../results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
# Save raw results
saveRDS(
  optimal_sequencing_depths_across_param_sweep,
  file = file.path(results_dir, "raw_results.rds")
)
# Extract results from each parameter sweep to get summary for plotting
opt_seq_depth_summary <- param_grid %>%
  mutate(
    optimal_depth = map_dbl(
      optimal_sequencing_depths_across_param_sweep,
      "optimal_depth"
    ),
    achieved_prob = map_dbl(
      optimal_sequencing_depths_across_param_sweep,
      "final_prob"
    ),
    annual_reads = optimal_depth * 52
  )

# Turn to a tibble for saving
opt_seq_depth_summary <- as_tibble(opt_seq_depth_summary)
write_tsv(opt_seq_depth_summary, "../results/summarized_results.tsv")
