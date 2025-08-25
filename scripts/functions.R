# Import packages
library(tidyverse)
library(furrr)
library(dplyr)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

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
      sampled_shedders <- rbinom(1, P, S_t / N)
      expected_reads <- sampled_shedders / P * mew * D
      weekly_reads <- rnbinom(1, size = 1 / alpha, mu = expected_reads)
      total_reads <- total_reads + weekly_reads
      # Detection occurred
      if (total_reads >= theta) {
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

# =============================================================================
# MAIN FUNCTIONS FOR RUNNING SIMULATIONS
# =============================================================================

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
