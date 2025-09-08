# Import packages
library(docstring)

# =============================================================================
# SINGLE SIMULATION FUNCTIONS
# =============================================================================

calc_infections <- function(t, I1, r) {
  #' Calculate the number of infected people at time t
  I1 * exp(r * (t - 1))
}

calc_shedding <- function(t, I1, r, shedding_weeks) {
  #' Calculate how many people are shedding viral reads at any time t
  start_week <- max(1, t - shedding_weeks + 1)
  end_week <- t
  weeks <- start_week:end_week
  infections <- sapply(weeks, function(w) calc_infections(w, I1, r))
  sum(infections)
}

calc_cumulative_incidence <- function(week, I1, r, N) {
  #' Function to calculate cumulative incidence as percentage of population

  # Calculate total infections from week 1 to target week
  if (r == 0) {
    total_infections <- I1 * week
  } else {
    total_infections <- I1 * (exp(r * week) - 1) / (exp(r) - 1)
  }
  # Return as percentage of population
  return(100 * total_infections / N)
}

run_single_simulation_of_outbreak <- function(
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
  #' Function to simulate one outbreak

  # Run epidemic week by week until detection
  week <- 1
  total_pathogen_reads <- 0
  while (week <= max_weeks && total_pathogen_reads < theta) {
    # Calculate shedding population for current week
    S_t <- calc_shedding(week, I1, r, shedding_weeks)
    S_t <- min(S_t, N / 4) # Cap at 1/4 population size since exponential growth slows down much earlier
    if (S_t > 0 && P < N) {
      sampled_shedders <- rbinom(1, P, S_t / N)
      expected_reads <- sampled_shedders / P * mew * D
      weekly_reads <- rnbinom(1, size = 1 / alpha, mu = expected_reads)
      total_pathogen_reads <- total_pathogen_reads + weekly_reads
      # Detection occurred
      if (total_pathogen_reads >= theta) {
        cumulative_incidence <- calc_cumulative_incidence(week, I1, r, N)
        return(list(
          detected = TRUE,
          detection_week = week,
          cumulative_incidence = cumulative_incidence,
          total_pathogen_reads = total_pathogen_reads
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
    total_pathogen_reads = total_pathogen_reads
  ))
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

run_simulations_at_given_sequencing_depth <- function(
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
  #' Function to simulate many outbreaks, then calculate detection probability and summary statistics

  # Run many simulations
  results <- replicate(
    n_sims,
    {
      run_single_simulation_of_outbreak(
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


calc_detection_probability_from_all_simulations <- function(
  detection_results,
  threshold_incidence
) {
  #' Function to calculate detection probability from saved results

  # Calculate P(detect | X% incidence) from saved results
  # detection_results: output from run_simulations_at_given_sequencing_depth()
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
# MAIN FUNCTIONS
# =============================================================================

find_optimal_sequencing_depth_with_binary_search <- function(
  target_cumulative_incidence,
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
    target_cumulative_incidence,
    target_prob * 100
  ))
  #' Using binary search, select a sequencing depth, then run simulations
  #' check if the target cumulative incidence is reached, if so
  #' adjust sequencing depth. Continue until bounds for
  #' binary search are below tolerance threshold.

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
    results <- run_simulations_at_given_sequencing_depth(
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
    p_detect <- calc_detection_probability_from_all_simulations(
      results,
      target_cumulative_incidence
    )
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
  final_results <- run_simulations_at_given_sequencing_depth(
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
  final_p_detect <- calc_detection_probability_from_all_simulations(
    final_results,
    target_cumulative_incidence
  )
  cat(sprintf("\nOptimal depth: %.2e reads/week\n", high))
  cat(sprintf(
    "Final P(detect | %.1f%% incidence) = %.2f%%\n",
    target_cumulative_incidence,
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

calc_weekly_operational_cost <- function(
  P,
  D,
  C_samp,
  C_samp_proc,
  C_seq,
  C_data_proc
) {
  #' Calculate weekly operational costs (C_weekly_ops)
  #' @param P Number of sampled individuals
  #' @param D Sequencing depth in billions of reads
  #' @param C_samp Cost per sample acquisition
  #' @param C_samp_proc Cost of sample processing (fixed)
  #' @param C_seq Cost per billion reads of sequencing
  #' @param C_data_proc Cost per billion reads of data processing

  return((C_samp * P) + C_samp_proc + D * (C_seq + C_data_proc))
}

calc_storage_unit_cost <- function(D, C_data_store) {
  #' Calculate cost to store one week's data for one week
  #' @param D Sequencing depth in billions of reads
  #' @param C_data_store Cost per billion reads to store for one week

  return(D * C_data_store)
}

calc_cumulative_storage_cost <- function(C_storage_unit, weeks = 52) {
  #' Calculate cumulative storage cost over the year
  #' Storage accumulates: week 1 data stored for 52 weeks,
  #' week 2 data for 51 weeks, etc.
  #' @param C_storage_unit Cost to store one week's data for one week
  #' @param weeks Number of weeks (default 52)

  # Sum of arithmetic series: sum(1 to 52) = 52 * 53 / 2
  total_storage_weeks <- weeks * (weeks + 1) / 2
  return(C_storage_unit * total_storage_weeks)
}

calc_total_cost_per_year <- function(
  seq_depth,
  cost_of_seq,
  cost_of_proc,
  cost_of_sample,
  num_samples,
  cost_of_data_proc,
  cost_of_data_store,
  cost_of_labor
) {
  #' Calculate the total cost per year of detecting this HIV-like pathogen
  #'
  #' @param seq_depth Sequencing depth (D in billions of reads)
  #' @param cost_of_seq Cost per billion reads of sequencing (C_seq)
  #' @param cost_of_proc Cost of sample processing per week (C_samp_proc)
  #' @param cost_of_sample Cost per sample acquisition (C_samp)
  #' @param num_samples Number of samples per week (P)
  #' @param cost_of_data_proc Cost per billion reads of data processing (C_data_proc)
  #' @param cost_of_data_store Cost per billion reads to store for one week (C_data_store)
  #' @param cost_of_labor Annual labor cost (C_labor)

  # Calculate weekly operational cost (constant each week)
  C_weekly_ops <- calc_weekly_operational_cost(
    P = num_samples,
    D = seq_depth,
    C_samp = cost_of_sample,
    C_samp_proc = cost_of_proc,
    C_seq = cost_of_seq,
    C_data_proc = cost_of_data_proc
  )

  # Calculate storage unit cost (cost to store one week's data for one week)
  C_storage_unit <- calc_storage_unit_cost(
    D = seq_depth,
    C_data_store = cost_of_data_store
  )

  # Calculate cumulative storage cost over 52 weeks
  cumulative_storage <- calc_cumulative_storage_cost(
    C_storage_unit = C_storage_unit,
    weeks = 52
  )

  # Total annual cost = 52 weeks of operations + cumulative storage + labor
  C_annual <- (C_weekly_ops * 52) + cumulative_storage + cost_of_labor

  return(C_annual)
}
