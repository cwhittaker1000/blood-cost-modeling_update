# Load required libraries
library(readr); library(purrr); library(future); library(future.apply); library(dplyr)
library(progressr); library(cli)

# Sourcing required function
source("scripts/functions.R")

handlers(global = TRUE)
handlers("txtprogressbar")

# =============================================================================
# PARAMETERS (shared across all sensitivity runs)
# =============================================================================
mu <- 1.26 * 10^-5
cv_squared <- 0.69

total_population_size <- 3.4e8
number_of_weeks_shedding <- 12
weekly_growth_rate <- 0.0155
initial_infections <- 100
read_threshold <- 100

target_cumulative_incidences <- c(
  0.000001, 0.00001, 0.00005, 0.0001, 0.0005, 0.001
)
number_of_individuals_sampled <- c(2000, 10000, 50000, 100000)
mus <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, mu)
prob_infected_donate_values <- c(0.1, 0.2, 0.25, 0.3333333, 0.5, 1)

workers <- max(1, parallel::detectCores() - 2)
plan(multisession, workers = workers)

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir)

# =============================================================================
# SENSITIVITY SWEEP OVER infectious_weeks
# =============================================================================
# The baseline (12 weeks) produces the "main" results; additional values
# produce sensitivity files with a suffix.
infectious_weeks_values <- c(12, 52)

for (iw in infectious_weeks_values) {
  
  cat(sprintf("\n=== Running infectious_weeks = %d ===\n", iw))
  
  # Derive R0 consistent with the target growth rate for this D
  doubling_time <- log(2) / weekly_growth_rate
  R0 <- calc_R0_from_doubling_time(doubling_time, iw)
  stopifnot(abs(calc_r_from_R0(R0, iw) - weekly_growth_rate) < 1e-6)
  cat(sprintf("  R0 = %.4f (doubling time = %.1f weeks)\n", R0, doubling_time))
  
  # Build file suffixes: baseline (12 wk) gets no suffix, others get one
  if (iw == number_of_weeks_shedding) {
    suffix <- ""
  } else {
    suffix <- sprintf("_%dwkSensitivity", iw)
  }
  
  rds_file <- file.path(results_dir, sprintf("simulation_results_renewal_probdonates%s.rds", suffix))
  tsv_file <- file.path(results_dir, sprintf("summarized_simulation_results_renewal_probdonates%s.tsv", suffix))
  
  # Skip if output already exists (remove files to force re-run)
  if (file.exists(tsv_file)) {
    cat(sprintf("  Output already exists: %s — skipping.\n", tsv_file))
    next
  }
  
  param_grid <- expand.grid(
    target_cumulative_incidence = target_cumulative_incidences,
    individuals_sampled         = number_of_individuals_sampled,
    mu                          = mus,
    prob_infected_donates       = prob_infected_donate_values,
    stringsAsFactors            = FALSE
  )
  
  raw_results <- with_progress({
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
          infectious_weeks            = iw,
          total_population_size       = total_population_size,
          prob_infected_donates       = row[["prob_infected_donates"]],
          growth_model                = "deterministic_renewal",
          n_sims                      = 250
        )
        p()
        result
      },
      future.seed = TRUE
    )
  })
  
  saveRDS(raw_results, file = rds_file)
  
  opt_seq_depth_summary <- param_grid %>%
    mutate(
      optimal_depth  = map_dbl(raw_results, "optimal_depth"),
      achieved_prob  = map_dbl(raw_results, "final_prob"),
      converged      = map_lgl(raw_results, "converged"),
      relative_error = map_dbl(raw_results, "relative_error"),
      annual_reads   = optimal_depth * 52
    ) %>%
    as_tibble()
  
  write_tsv(opt_seq_depth_summary, tsv_file)
  cat(sprintf("  Saved: %s\n", tsv_file))
}

cat("\nAll runs complete.\n")