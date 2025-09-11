source("functions.R")
library(readr)
library(purrr)
library(future)
library(future.apply)
library(dplyr)

# =============================================================================
# LOAD IN DATA AND ESTIMATE MEW/ALPHA
# =============================================================================

# Load in data from Piantadosi et al. 2019
# Data contains HIV reads from plasma samples of three individuals within two weeks of infection
# For individuals with multi-day samples, we used the last day for each individual
# Sequencing depth was reported as a range (40M-250M reads), but not specified per sample
# We used the geometric mean of this range (100M reads) for all samples
data <- read_csv("../data/new-hiv-data.csv", show_col_types = FALSE) %>%
  mutate(hiv_ra = hiv_reads / read_depth)

# Parameterize pathogen shedding rate (mu) using arithmetic mean relative abundance
# mu represents the expected fraction of sequencing reads that are viral in an infected person's sample
# Calculated as arithmetic mean of (HIV reads / total reads) across three individuals = 1.26 × 10^-5
# This parameter scales the expected number of viral reads based on how many infected people are sampled
mu <- data %>%
  summarize(mean = mean(hiv_ra)) %>%
  pull(mean)

# Parameterize overdispersion parameter using Coefficient of Variation Squared of relative abundances
# cv_squared captures biological and technical variability in viral shedding between individuals
# Calculated as (SD of relative abundances / mean of relative abundances)^2 = 0.69
# Used as 1/size parameter in negative binomial distribution to model count overdispersion
cv_squared <- data %>%
  summarize(
    mean_ra = mean(hiv_ra, na.rm = TRUE),
    sd_ra = sd(hiv_ra, na.rm = TRUE),
    cv_squared = (sd_ra / mean_ra)^2
  ) %>%
  pull(cv_squared)

# =============================================================================
# DEFINE OTHER PARAMETERS FOR SIMULATIONS
# =============================================================================

# Constants for simulations
total_population_size <- 3e8 # Approximate US population
number_of_weeks_shedding <- 12 # The number of weeks an infected person sheds
weekly_growth_rate <- 0.0155 # weekly growth rate based on HIV
initial_infections <- 100 # Number of initial infections
read_threshold <- 100 # Number of viral reads required to detect the pathogen

# Varying parameters for simulations
target_cumulative_incidences <- c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1)
number_of_individuals_sampled <- c(2000, 10000, 50000, 100000)
# Use range of mu to account for variance in sample processing/library prep and include our original estimate
mus <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, mu)

# Setup for running simulations using multiple cores
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)
param_grid <- expand.grid(
  target_cumulative_incidence = target_cumulative_incidences,
  individuals_sampled = number_of_individuals_sampled,
  mu = mus,
  stringsAsFactors = FALSE
)

# =============================================================================
# RUN EXPENSIVE SIMULATIONS
# =============================================================================

minimum_sequencing_depth_for_detection_across_params <- future_apply(
  X = param_grid, # each row is one job
  MARGIN = 1, # tells the function the dimension along which to apply the function, in this case each row
  FUN = function(row) {
    minimize_sequencing_depth_given_detection_constraint_using_binary_search(
      target_cumulative_incidence = row["target_cumulative_incidence"],
      target_prob = 0.95,
      weekly_growth_rate = weekly_growth_rate,
      initial_infections = initial_infections,
      individuals_sampled = row["individuals_sampled"],
      read_detection_threshold = read_threshold,
      mu = row["mu"],
      cv_squared = cv_squared,
      number_of_weeks_shedding = number_of_weeks_shedding,
      total_population_size = total_population_size,
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
  minimum_sequencing_depth_for_detection_across_params,
  file = file.path(results_dir, "simulation_results.rds")
)

# Extract results from each parameter sweep to get summary for plotting
opt_seq_depth_summary <- param_grid %>%
  mutate(
    optimal_depth = map_dbl(
      minimum_sequencing_depth_for_detection_across_params,
      "optimal_depth"
    ),
    achieved_prob = map_dbl(
      minimum_sequencing_depth_for_detection_across_params,
      "final_prob"
    ),
    converged = map_lgl(
      minimum_sequencing_depth_for_detection_across_params,
      "converged"
    ),
    relative_error = map_dbl(
      minimum_sequencing_depth_for_detection_across_params,
      "relative_error"
    ),
    annual_reads = optimal_depth * 52
  )

# Turn to a tibble for saving
opt_seq_depth_summary <- as_tibble(opt_seq_depth_summary)
write_tsv(opt_seq_depth_summary, "../results/summarized_simulation_results.tsv")
