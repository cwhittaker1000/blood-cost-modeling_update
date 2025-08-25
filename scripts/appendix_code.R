# Import packages
library(tidyverse)
library(furrr)
library(scales)
library(future)
library(future.apply)
library(dplyr)

############################
# Start here if from scratch
############################

# Function to solve for T given r, I1, and X
solve_time_period <- function(r, I1, X_pct, N) {
  target_total <- (X_pct / 100) * N
  if (r == 0) {
    # Constant infections case
    return(target_total / I1)
  } else {
    # Exponential growth case
    # T = (1/r) * ln(1 + (X*N*(e^r - 1))/(100*I1))
    exp_r <- exp(r)
    ratio <- (X_pct * N * (exp_r - 1)) / (100 * I1)
    T <- (1 / r) * log(1 + ratio)
    return(T)
  }
}

# Function to calculate infections at time t
calc_infections <- function(t, I1, r) {
  I1 * exp(r * (t - 1))
}

# Function to calculate shedding population at time t
calc_shedding <- function(t, I1, r, shedding_weeks) {
  start_week <- max(1, t - shedding_weeks + 1)
  end_week <- t
  weeks <- start_week:end_week
  infections <- sapply(weeks, function(w) calc_infections(w, I1, r))
  sum(infections)
}

# Function to calculate cumulative incidence as percentage of population
calc_cumulative_incidence <- function(week, I1, r, N) {
  # Calculate total infections from week 1 to target week
  if (r == 0) {
    total_infections <- I1 * week
  } else {
    total_infections <- I1 * (exp(r * week) - 1) / (exp(r) - 1)
  }
  # Return as percentage of population
  return(100 * total_infections / N)
}

# Function to simulate detection timing distribution
simulate_detection_distribution <- function(
  r,
  I1,
  P,
  D,
  theta,
  mew,
  alpha,
  shedding_weeks,
  N,
  max_weeks = 1000
) {
  # Run epidemic week by week until detection
  week <- 1
  total_reads <- 0
  while (week <= max_weeks && total_reads < theta) {
    # Calculate shedding population for current week
    S_t <- calc_shedding(week, I1, r, shedding_weeks)
    S_t <- round(S_t)
    S_t <- min(S_t, N) # Cap at population size
    if (S_t > 0 && P < N) {
      # For large populations, approximate with binomial
      sampled_shedders <- rbinom(1, P, S_t / N)
      # Generate reads from negative binomial distribution
      expected_reads <- sampled_shedders / P * mew * D
      weekly_reads <- rnbinom(1, size = 1 / alpha, mu = expected_reads)
      total_reads <- total_reads + weekly_reads
      if (total_reads >= theta) {
        # Detection occurred!
        cumulative_incidence <- calc_cumulative_incidence(week, I1, r, N)
        return(list(
          detected = TRUE,
          detection_week = week,
          cumulative_incidence = cumulative_incidence,
          total_reads = total_reads
        ))
      }
    }
    week <- week + 1
  }
  # No detection within max_weeks
  return(list(
    detected = FALSE,
    detection_week = NA,
    cumulative_incidence = NA,
    total_reads = total_reads
  ))
}


# Function to analyze detection distribution across many simulations
analyze_detection_distribution <- function(
  r,
  I1,
  P,
  D,
  theta,
  mew,
  alpha,
  shedding_weeks,
  N,
  n_sims = 1000
) {
  # Run many simulations
  results <- replicate(
    n_sims,
    {
      simulate_detection_distribution(
        r,
        I1,
        P,
        D,
        theta,
        mew,
        alpha,
        shedding_weeks,
        N
      )
    },
    simplify = FALSE
  )
  # Extract detection data
  detected_sims <- results[sapply(results, function(x) x$detected)]
  if (length(detected_sims) == 0) {
    return(list(detection_rate = 0, summary_stats = NA))
  }
  # Extract cumulative incidence at detection
  detection_incidences <- sapply(detected_sims, function(x) {
    x$cumulative_incidence
  })
  # Calculate summary statistics
  summary_stats <- list(
    detection_rate = length(detected_sims) / n_sims,
    mean_incidence = mean(detection_incidences),
    median_incidence = median(detection_incidences),
    q25_incidence = quantile(detection_incidences, 0.25),
    q75_incidence = quantile(detection_incidences, 0.75),
    min_incidence = min(detection_incidences),
    max_incidence = max(detection_incidences)
  )
  return(list(
    detection_rate = summary_stats$detection_rate,
    summary_stats = summary_stats,
    raw_data = detection_incidences
  ))
}


# Function to calculate detection probability from saved results
calc_detection_probability <- function(detection_results, threshold_incidence) {
  # Calculate P(detect | X% incidence) from saved results
  # detection_results: output from analyze_detection_distribution()
  # threshold_incidence: X% incidence threshold
  if (detection_results$detection_rate == 0) {
    return(0)
  }
  # Count simulations that detected at or before threshold
  detected_early <- sum(detection_results$raw_data <= threshold_incidence)
  total_sims <- length(detection_results$raw_data) /
    detection_results$detection_rate
  return(detected_early / total_sims)
}

# Function to find optimal sequencing depth using binary search
find_optimal_depth <- function(
  target_incidence,
  target_prob = 0.95,
  r,
  I1,
  P,
  theta,
  mew,
  alpha,
  shedding_weeks,
  N,
  n_sims = 500,
  tolerance = 1e4
) {
  cat(sprintf(
    "Finding optimal depth for %.1f%% cumulative incidence at %.0f%% detection probability\n",
    target_incidence,
    target_prob * 100
  ))
  # Binary search bounds
  low <- 1e4 # Start with very low depth
  high <- 1e12 # Maximum reasonable depth
  # Track results for visualization
  tested_depths <- c()
  tested_probs <- c()
  while (high - low > tolerance) {
    mid <- round((low + high) / 2)
    cat(sprintf("Testing depth %.2e... ", mid))
    # Run simulations at this depth
    results <- analyze_detection_distribution(
      r,
      I1,
      P,
      mid,
      theta,
      mew,
      alpha,
      shedding_weeks,
      N,
      n_sims
    )
    p_detect <- calc_detection_probability(results, target_incidence)
    cat(sprintf("P(detect) = %.2f%%\n", p_detect * 100))
    # Store for visualization
    tested_depths <- c(tested_depths, mid)
    tested_probs <- c(tested_probs, p_detect)
    if (p_detect >= target_prob) {
      high <- mid # Can use less sequencing
    } else {
      low <- mid # Need more sequencing
    }
  }
  # Run final simulation at selected depth to confirm
  final_results <- analyze_detection_distribution(
    r,
    I1,
    P,
    high,
    theta,
    mew,
    alpha,
    shedding_weeks,
    N,
    n_sims * 2
  )
  final_p_detect <- calc_detection_probability(final_results, target_incidence)
  cat(sprintf("\nOptimal depth: %.2e reads/week\n", high))
  cat(sprintf(
    "Final P(detect | %.1f%% incidence) = %.2f%%\n",
    target_incidence,
    final_p_detect * 100
  ))
  return(list(
    optimal_depth = high,
    final_prob = final_p_detect,
    tested_depths = tested_depths,
    tested_probs = tested_probs,
    final_results = final_results
  ))
}


############################
# Start here if already run
############################

find_total_cost <- function(
  seq_depth,
  cost_of_seq,
  cost_of_proc,
  cost_of_sample,
  num_samples
) {
  return(
    seq_depth * (cost_of_seq + cost_of_proc) + cost_of_sample * num_samples
  )
}


# load in data
data <- read_csv("../data/new-hiv-data.csv", show_col_types = FALSE)

our_data <- data %>%
  mutate(hiv_ra = hiv_reads / read_depth)

our_data

estimate_mew <- our_data %>%
  summarize(geo_mean = exp(mean(log(hiv_ra), na.rm = TRUE))) %>%
  pull(geo_mean)

# Calculate overdispersion using Coefficient of Variation Squared of relative abundances
cv_squared <- our_data %>%
  summarize(
    mean_ra = mean(hiv_ra, na.rm = TRUE),
    sd_ra = sd(hiv_ra, na.rm = TRUE),
    cv_squared = (sd_ra / mean_ra)^2
  ) %>%
  pull(cv_squared)


## Model parameters

# CONSTANTS
N <- 3e8 # US population
#theta <- 100  # detection threshold (reads)
shedding_weeks <- 12
r <- 0.0155 # weekly growth based on HIV
I1 <- 100 # 100 initial infections
alpha <- cv_squared # Use CV² directly as overdispersion parameter

# Find optimal sequencing depths for different target incidences
target_incidences <- c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1) # Different cumulative incidence targets
pool_sizes <- c(1200, 6000, 12000, 60000, 120000) # Different cumulative incidence targets

#mews <- c(
#  estimate_mew * 1e-2,
#  estimate_mew * 1e-1,
#  estimate_mew,
#  estimate_mew * 1e1,
#  estimate_mew * 1e2
#)
mews <- 10^seq(-2, -8, by = -1)

read_thresholds <- c(100, 1000, 10000)
optimal_results <- list()


# Parallel
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)

param_grid <- expand.grid(
  target_incidence = target_incidences,
  pool_size = pool_sizes,
  mew = mews,
  read_threshold = read_thresholds,
  stringsAsFactors = FALSE
)

new_optimal_results <- future_apply(
  X = param_grid, # each row is one job
  MARGIN = 1,
  FUN = function(row) {
    find_optimal_depth(
      target_incidence = row["target_incidence"],
      target_prob = 0.95,
      r = r,
      I1 = I1,
      P = row["pool_size"],
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

# After the future_apply() call finishes
results_dir <- "../results" # adjust as needed

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
saveRDS(
  new_optimal_results,
  file = file.path(results_dir, "lenni_final_optimal_results.rds")
)

rm(new_optimal_results)

optimal_results <- readRDS("../results/lenni_final_optimal_results.rds")

optimal_summary <- param_grid %>%
  mutate(
    optimal_depth = map_dbl(optimal_results, "optimal_depth"),
    achieved_prob = map_dbl(optimal_results, "final_prob"),
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

optimal_summary <- as_tibble(optimal_summary)

normal_data <- optimal_summary %>%
  filter(mew == 1e-6, read_threshold == 100) %>%
  mutate(
    # Flag values that hit the ceiling
    hit_ceiling = optimal_depth >= 1e12
  ) %>%
  filter(!hit_ceiling)

depth_plot <- ggplot(
  normal_data,
  aes(
    x = target_incidence,
    y = optimal_depth,
    colour = factor(pool_size),
    group = pool_size
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
proc_c <- 2500 / 1e9
sample_c <- 15

# 1. add a new column with the cost for each depth
# Filter out non-converged results for cost plot
normal_data_with_cost <- normal_data %>%
  filter(!hit_ceiling) %>% # Exclude points that didn't converge
  mutate(
    total_cost_per_year = find_total_cost(
      seq_depth = optimal_depth,
      cost_of_seq = seq_c,
      cost_of_proc = proc_c,
      cost_of_sample = sample_c,
      num_samples = pool_size
    ) *
      52
  )

cost_plot <- ggplot(
  normal_data_with_cost,
  aes(
    x = target_incidence,
    y = total_cost_per_year,
    colour = factor(pool_size),
    group = pool_size
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

#optimal_summary %>% select(target_incidence, pool_size, total_cost_per_year, q25_incidence,q75_incidence) %>%
#  filter(pool_size == 1024)

normal_all_data <- optimal_summary %>%
  filter(read_threshold == 100) %>%
  mutate(
    # Flag values that hit the ceiling
    hit_ceiling = optimal_depth >= 1e12
  ) %>%
  filter(!hit_ceiling)

facet_depth_plot <- ggplot(
  normal_all_data,
  aes(
    x = target_incidence,
    y = optimal_depth,
    colour = factor(pool_size),
    group = pool_size
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
seq_c <- 2500 / 1e9
proc_c <- 2500 / 1e9
sample_costs <- c(3, 12, 50, 100)

# Create data frame with all combinations of sample costs
normal_data_with_varied_costs <- normal_all_data %>%
  filter(!hit_ceiling) %>% # Exclude points that didn't converge
  crossing(sample_cost = sample_costs) %>% # Create all combinations
  mutate(
    total_cost_per_year = find_total_cost(
      seq_depth = optimal_depth,
      cost_of_seq = seq_c,
      cost_of_proc = proc_c,
      cost_of_sample = sample_cost,
      num_samples = pool_size
    ) *
      52
  )

# Create faceted plot by sample cost
facet_cost_plot <- ggplot(
  normal_data_with_varied_costs,
  aes(
    x = target_incidence,
    y = total_cost_per_year,
    colour = factor(pool_size),
    group = pool_size
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
