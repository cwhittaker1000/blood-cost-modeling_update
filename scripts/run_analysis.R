# Import packages
# TODO: Look into early stopping of binary search if we don't think we'll converge.

source("functions.R")
library(readr)
library(purrr)
library(future)
library(future.apply)
library(dplyr)

# conditional on me keeping plot
library(ggplot2)

# =============================================================================
# LOAD IN DATA AND ESTIMATE MEW/ALPHA
# =============================================================================

# Load in data
our_data <- read_csv("../data/new-hiv-data.csv", show_col_types = FALSE) %>%
  mutate(hiv_ra = hiv_reads / read_depth)

# Use arithmetic mean to estimate mew
# 1.26e-5
estimate_mew <- our_data %>%
  summarize(arith_mean = mean(hiv_ra)) %>%
  pull(arith_mean)

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
batch_sizes <- c(2000, 10000, 50000, 100000)

# Change color, and keep 1e-3 to 1e-8, include actual estimate
mews <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, estimate_mew)

# Fix at 100
read_thresholds <- c(100)

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
# RUN EXPENSIVE SIMULATIONS
# =============================================================================

optimal_sequencing_depths_across_param_sweep <- future_apply(
  X = param_grid, # each row is one job
  MARGIN = 1, # tells the function the dimension along which to apply the function, in this case each row
  FUN = function(row) {
    find_optimal_sequencing_depth_with_binary_search(
      target_cumulative_incidence = row["target_cumulative_incidence"],
      target_prob = 0.95,
      r = r,
      I = I,
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
  file = file.path(results_dir, "simulation_results.rds")
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
write_tsv(opt_seq_depth_summary, "../results/summarized_simulation_results.tsv")

# =============================================================================
# Data validation
# =============================================================================

our_data %>%
  summarize(
    mean = mean(hiv_reads),
    variance = var(hiv_reads),
    sd = sd(hiv_reads),
    cv = sd / mean
  )

ggplot(aes(y = hiv_reads), data = our_data) +
  geom_histogram()

D <- our_data %>% pull(read_depth) %>% unique()

# TEMP
p <- dnbinom(x = seq(4000), size = 1 / alpha, mu = estimate_mew * D)

ggplot(aes(x = seq(4000), y = p), data = NULL) +
  geom_line() +
  geom_point(aes(x = hiv_reads, y = 0), our_data)


pnbinom(50, size = 1 / alpha, mu = estimate_mew * D)
# TEMP
