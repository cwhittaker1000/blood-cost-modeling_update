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
hiv_mu <- 1.26 * 10^-5
cv_squared <- 0.69

total_population_size <- 3.4e8
number_of_weeks_shedding <- 12
weekly_growth_rate <- 0.0155
initial_infections <- 100
read_threshold <- 100

# Reduced parameter grid relative to run_analysis_renewal.R:
#   - mu: HIV central estimate and one order of magnitude either side (3 values)
#   - individuals_sampled: 10,000 only
#   - prob_infected_donates: 1 only
#   - target_cumulative_incidences: full sweep (unchanged)
target_cumulative_incidences <- c(
  0.000001, 0.00001, 0.00005, 0.0001, 0.0005, 0.001
)
number_of_individuals_sampled <- 10000
mus <- c(hiv_mu / 10, hiv_mu, hiv_mu * 10)   # 1.26e-6, 1.26e-5, 1.26e-4
prob_infected_donate_values <- 1

workers <- max(1, parallel::detectCores() - 8)
plan(multisession, workers = workers)

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir)

# =============================================================================
# SENSITIVITY SWEEP OVER infectious_weeks
# =============================================================================
# The baseline (12 weeks) produces the "main" results; additional values
# produce sensitivity files with a suffix.
infectious_weeks_values <- c(12)

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
  
  rds_file <- file.path(results_dir, sprintf("simulation_results_renewal_hivmuscan%s.rds", suffix))
  tsv_file <- file.path(results_dir, sprintf("summarized_simulation_results_renewal_hivmuscan%s.tsv", suffix))
  
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
  
  cat(sprintf("  Parameter grid: %d combinations\n", nrow(param_grid)))
  
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

# =============================================================================
# Compute annual cost and component breakdown for each mu in the hivmuscan
# analysis (3 mus x 6 cumulative incidences).
# =============================================================================

library(readr); library(dplyr); library(scales)

source("scripts/functions.R")

# -----------------------------------------------------------------------------
# Cost parameters (match plot_figures.R)
# -----------------------------------------------------------------------------
cost_per_individual_sample              <- 5.00
cost_sample_processing_per_week         <- 500
cost_per_billion_reads_sequencing       <- 2500
cost_per_billion_reads_data_processing  <- 1
cost_per_billion_reads_storage_per_week <- 0.25
annual_labor_cost                       <- 650000

hiv_mu <- 1.26e-5

# -----------------------------------------------------------------------------
# Load results (baseline infectious_weeks = 12; swap suffix for sensitivity)
# -----------------------------------------------------------------------------
input_tsv  <- "results/summarized_simulation_results_renewal_hivmuscan.tsv"
output_tsv <- "results/cost_breakdown_hivmuscan.tsv"

results <- read_tsv(input_tsv, show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Component-by-component annual cost
#   - sample_acquisition: per-donor cost * donors/wk * 52
#   - sample_processing : fixed weekly cost * 52
#   - sequencing        : reads/wk * $/B reads * 52
#   - data_processing   : reads/wk * $/B reads * 52
#   - data_storage      : cumulative (triangular) storage over 52 weeks
#   - labor             : annual flat cost
# -----------------------------------------------------------------------------
cost_breakdown <- results %>%
  mutate(
    depth_billion      = optimal_depth / 1e9,
    sample_acquisition = cost_per_individual_sample * individuals_sampled * 52,
    sample_processing  = cost_sample_processing_per_week * 52,
    sequencing         = depth_billion * cost_per_billion_reads_sequencing * 52,
    data_processing    = depth_billion * cost_per_billion_reads_data_processing * 52,
    data_storage       = calc_cumulative_storage_cost(
      calc_storage_unit_cost(depth_billion, cost_per_billion_reads_storage_per_week),
      num_weeks_in_period = 52
    ),
    labor       = annual_labor_cost,
    annual_cost = sample_acquisition + sample_processing + sequencing +
      data_processing + data_storage + labor,
    # Label mus by their ratio to the HIV central estimate (robust to
    # round-trip float noise from TSV write/read)
    mu_label = case_when(
      abs(mu / hiv_mu - 0.1) < 0.01 ~ "1 OoM below HIV (1.26e-6)",
      abs(mu / hiv_mu - 1)   < 0.01 ~ "HIV central    (1.26e-5)",
      abs(mu / hiv_mu - 10)  < 0.01 ~ "1 OoM above HIV (1.26e-4)",
      TRUE                          ~ sprintf("%.2e", mu)
    )
  ) %>%
  arrange(mu, target_cumulative_incidence) %>%
  select(mu, mu_label, target_cumulative_incidence,
         individuals_sampled, prob_infected_donates,
         optimal_depth, depth_billion, annual_reads,
         sample_acquisition, sample_processing, sequencing,
         data_processing, data_storage, labor, annual_cost,
         converged, achieved_prob, relative_error)

write_tsv(cost_breakdown, output_tsv)
cat(sprintf("Saved full cost breakdown to: %s\n\n", output_tsv))

# -----------------------------------------------------------------------------
# Readable console summary: one table per mu
# -----------------------------------------------------------------------------
mu_labels_ordered <- c(
  "1 OoM below HIV (1.26e-6)",
  "HIV central    (1.26e-5)",
  "1 OoM above HIV (1.26e-4)"
)

for (lbl in mu_labels_ordered) {
  sub <- cost_breakdown %>% filter(mu_label == lbl)
  if (nrow(sub) == 0) next
  
  cat(sprintf("=== mu = %s ===\n", lbl))
  
  display <- sub %>%
    transmute(
      CI            = percent(target_cumulative_incidence, accuracy = 0.0001),
      `Seq/wk (B)`  = comma(depth_billion, accuracy = 0.01),
      `Annual cost` = paste0("$", comma(round(annual_cost))),
      `Acq %`       = percent(sample_acquisition / annual_cost, accuracy = 1),
      `Seq %`       = percent(sequencing         / annual_cost, accuracy = 1),
      `Labor %`     = percent(labor              / annual_cost, accuracy = 1),
      `Other %`     = percent((sample_processing + data_processing + data_storage)
                              / annual_cost, accuracy = 1),
      Converged     = converged
    )
  print(display, n = Inf)
  cat("\n")
}

# -----------------------------------------------------------------------------
# Cross-mu comparison: annual cost at each CI, one column per mu
# -----------------------------------------------------------------------------
cat("=== Annual cost by CI x mu ===\n")
cost_breakdown %>%
  mutate(`Annual cost` = paste0("$", comma(round(annual_cost)))) %>%
  select(target_cumulative_incidence, mu_label, `Annual cost`) %>%
  tidyr::pivot_wider(names_from = mu_label, values_from = `Annual cost`) %>%
  mutate(CI = percent(target_cumulative_incidence, accuracy = 0.0001)) %>%
  select(CI, everything(), -target_cumulative_incidence) %>%
  print(n = Inf)

