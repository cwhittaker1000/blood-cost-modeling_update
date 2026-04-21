# Load required libraries
library(readr); library(purrr); library(future); library(future.apply); library(dplyr)
library(progressr); library(cli)

# Sourcing required function
source("scripts/functions.R")

# Optional: pick a handler. "txtprogressbar" works everywhere;
# "cli" is prettier if you have the cli package.
handlers(global = TRUE)
# handlers("cli")
handlers("txtprogressbar")

# =============================================================================
# LOAD IN DATA AND ESTIMATE MEW/ALPHA
# =============================================================================

# Load in data from Piantadosi et al. 2019
# Data contains HIV reads from plasma samples of three individuals within two weeks of infection
# For individuals with multi-day samples, we used the last day for each individual
# Sequencing depth was reported as a range (40M-250M reads), but not specified per sample
# We used the geometric mean of this range (100M reads) for all samples
data <- read_csv("data/new-hiv-data.csv", show_col_types = FALSE) %>%
  mutate(hiv_ra = hiv_reads / read_depth)

# Parameterize pathogen shedding rate (mu) using arithmetic mean relative abundance
# mu represents the expected fraction of sequencing reads that are viral in an infected person's sample
# Calculated as arithmetic mean of (HIV reads / total reads) across three individuals = 1.26 × 10^-5
# This parameter scales the expected number of viral reads based on how many infected people are sampled
# mu <- data %>%
#   summarize(mean = mean(hiv_ra)) %>%
#   pull(mean)
mu <- 1.26 * 10^-5 

# Parameterize overdispersion parameter using Coefficient of Variation Squared of relative abundances
# cv_squared captures biological and technical variability in viral shedding between individuals
# Calculated as (SD of relative abundances / mean of relative abundances)^2 = 0.69
# Note: This ignores Poisson noise (which contributes to CV in the data but shouldn't be in the size parameter);
# however, in our low-abundance regime, this approximation should be fine
# Used as individuals_shedding/cv_squared as size parameter in negative binomial distribution
# This scaling accounts for variance reduction when sampling multiple infected individuals
# cv_squared <- data %>%
#   summarize(
#     mean_ra = mean(hiv_ra, na.rm = TRUE),
#     sd_ra = sd(hiv_ra, na.rm = TRUE),
#     cv_squared = (sd_ra / mean_ra)^2
#   ) %>%
#   pull(cv_squared)
cv_squared <- 0.69

# =============================================================================
# DEFINE OTHER PARAMETERS FOR SIMULATIONS
# =============================================================================

# Constants for simulations
total_population_size <- 3.4e8 # Approximate US population
number_of_weeks_shedding <- 12 # The number of weeks an infected person sheds
weekly_growth_rate <- 0.0155 # weekly growth rate based on HIV
initial_infections <- 100 # Number of initial infections
read_threshold <- 100 # Number of viral reads required to detect the pathogen

# The new functions.R interface parameterises growth via (R0, infectious_weeks)
# rather than weekly_growth_rate directly. For the "deterministic_simple" model,
# the tracker derives r internally from (R0, infectious_weeks) via the
# Euler-Lotka equation, so to exactly replicate the old r = 0.0155/week we
# convert it to a corresponding R0. Note: for the simple exponential model,
# the choice of infectious_weeks only affects the R0 <-> r conversion — the
# growth trajectory itself is purely I0 * exp(r*(t-1)) and doesn't depend on
# infectious_weeks. So any infectious_weeks value will reproduce the old
# dynamics as long as we compute R0 consistently with it. We use
# number_of_weeks_shedding here as a reasonable default.
infectious_weeks <- number_of_weeks_shedding
doubling_time <- log(2) / weekly_growth_rate
R0 <- calc_R0_from_doubling_time(doubling_time, infectious_weeks)
# Sanity check: tracker should derive back the same r we started with
stopifnot(abs(calc_r_from_R0(R0, infectious_weeks) - weekly_growth_rate) < 1e-6)

# Varying parameters for simulations
target_cumulative_incidences <- c(
  0.000001,
  0.00001,
  0.00005,
  0.0001,
  0.0005,
  0.001
)
number_of_individuals_sampled <- c(2000, 10000, 50000, 100000)
# Use range of mu to account for variance in sample processing/library prep and include our original estimate
mus <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, mu)

# Setup for running simulations using multiple cores
workers <- max(1, parallel::detectCores() - 2)
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

# minimum_sequencing_depth_for_detection_across_params <- future_apply(
#   X = param_grid, # each row is one job
#   MARGIN = 1, # tells the function the dimension along which to apply the function, in this case each row
#   FUN = function(row) {
#     minimize_sequencing_depth_given_detection_constraint_using_binary_search(
#       target_cumulative_incidence = row[["target_cumulative_incidence"]],
#       target_prob = 0.95,
#       R0 = R0,
#       initial_infections = initial_infections,
#       individuals_sampled = row[["individuals_sampled"]],
#       read_detection_threshold = read_threshold,
#       mu = row[["mu"]],
#       cv_squared = cv_squared,
#       number_of_weeks_shedding = number_of_weeks_shedding,
#       infectious_weeks = infectious_weeks,
#       total_population_size = total_population_size,
#       growth_model = "deterministic_simple",
#       n_sims = 500
#     )
#   },
#   future.seed = TRUE
# )
minimum_sequencing_depth_for_detection_across_params <- with_progress({
  p <- progressor(steps = nrow(param_grid))
  
  future_apply(
    X = param_grid,
    MARGIN = 1,
    FUN = function(row) {
      result <- minimize_sequencing_depth_given_detection_constraint_using_binary_search(
        target_cumulative_incidence = row[["target_cumulative_incidence"]],
        target_prob = 0.95,
        R0 = R0,
        initial_infections = initial_infections,
        individuals_sampled = row[["individuals_sampled"]],
        read_detection_threshold = read_threshold,
        mu = row[["mu"]],
        cv_squared = cv_squared,
        number_of_weeks_shedding = number_of_weeks_shedding,
        infectious_weeks = infectious_weeks,
        total_population_size = total_population_size,
        growth_model = "deterministic_simple",
        n_sims = 100 # 500
      )
      p()  # signal that one job is done
      result
    },
    future.seed = TRUE
  )
})

results_dir <- "results"
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
write_tsv(opt_seq_depth_summary, "results/summarized_simulation_results.tsv")
