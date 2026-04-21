# Load required libraries
library(readr); library(purrr); library(future); library(future.apply); library(dplyr)
library(progressr); library(cli)

# Sourcing required function
source("scripts/functions.R")

handlers(global = TRUE)
handlers("txtprogressbar")

# =============================================================================
# LOAD IN DATA AND ESTIMATE MU/CV
# =============================================================================
# Same as the simple-exponential national run for consistency.

mu <- 1.26 * 10^-5
cv_squared <- 0.69

# =============================================================================
# DEFINE PARAMETERS FOR SIMULATIONS
# =============================================================================
# All parameters held identical to the original national run except for the
# growth model. This isolates the effect of switching from simple exponential
# to deterministic renewal.

total_population_size <- 3.4e8 # Approximate US population
number_of_weeks_shedding <- 12 # The number of weeks an infected person sheds
weekly_growth_rate <- 0.0155   # weekly growth rate based on HIV
initial_infections <- 100      # Number of initial infections
read_threshold <- 100          # Number of viral reads required to detect the pathogen

# Renewal model uses infectious_weeks for transmission dynamics directly.
# Setting it equal to number_of_weeks_shedding as the baseline; this can be
# swept later as a sensitivity analysis. Unlike the simple-exponential model,
# this choice DOES affect the trajectory of the renewal model: it controls the
# length of the "cliff" during which seed infections build up to exponential
# growth.
infectious_weeks <- 52 # number_of_weeks_shedding
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
mus <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, mu)

# Range of values assuming different probability of donating conditional on being infected
prob_infected_donate_values <- c(0.1, 0.2, 0.25, 0.3333333, 0.5, 1)

# Setup for running simulations using multiple cores
workers <- max(1, parallel::detectCores() - 2)
plan(multisession, workers = workers)
param_grid <- expand.grid(
  target_cumulative_incidence = target_cumulative_incidences,
  individuals_sampled         = number_of_individuals_sampled,
  mu                          = mus,
  prob_infected_donates       = prob_infected_donate_values,
  stringsAsFactors            = FALSE
)

# =============================================================================
# RUN EXPENSIVE SIMULATIONS
# =============================================================================
# Same as the simple-exponential run but with growth_model = "deterministic_renewal".

minimum_sequencing_depth_for_detection_across_params <- with_progress({
  p <- progressor(steps = nrow(param_grid))
  
  future_apply(
    X = param_grid,
    MARGIN = 1,
    FUN = function(row) {
      result <- minimize_sequencing_depth_given_detection_constraint_using_binary_search(
        target_cumulative_incidence = row[["target_cumulative_incidence"]],
        target_prob                 = 0.95,
        R0                          = R0,
        initial_infections          = initial_infections,
        individuals_sampled         = row[["individuals_sampled"]],
        read_detection_threshold    = read_threshold,
        mu                          = row[["mu"]],
        cv_squared                  = cv_squared,
        number_of_weeks_shedding    = number_of_weeks_shedding,
        infectious_weeks            = infectious_weeks,
        total_population_size       = total_population_size,
        prob_infected_donates       = row[["prob_infected_donates"]],
        growth_model                = "deterministic_renewal",
        n_sims                      = 250 # 500 for production
      )
      p()
      result
    },
    future.seed = TRUE
  )
})

results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

saveRDS(
  minimum_sequencing_depth_for_detection_across_params,
  file = file.path(results_dir, "simulation_results_renewal_probdonates_52weekSensitivity.rds")
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

opt_seq_depth_summary <- as_tibble(opt_seq_depth_summary)
write_tsv(opt_seq_depth_summary, "results/summarized_simulation_results_renewal_probdonates_52weekSensitivity.tsv")