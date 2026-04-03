library(docstring)

create_outbreak_tracker <- function(
    initial_infections,
    weekly_growth_rate,
    max_weeks,
    stochastic_growth = FALSE
) {
  #' Create an outbreak tracker that manages infection dynamics
  #' @param initial_infections Number of infected people at week 1
  #' @param weekly_growth_rate Growth rate of the pathogen
  #' @param max_weeks Maximum number of weeks to track
  #' @param stochastic_growth If TRUE, use Poisson branching process; if FALSE, use deterministic exponential growth
  #' @return List of functions: advance_to_week, get_shedding, get_cumulative
  
  new_infections <- numeric(max_weeks)
  new_infections[1] <- initial_infections
  
  list(
    advance_to_week = function(week_t) {
      if (week_t == 1) return(invisible(NULL))
      if (new_infections[week_t] != 0) return(invisible(NULL))
      if (stochastic_growth) {
        new_infections[week_t] <<- rpois(
          1, new_infections[week_t - 1] * exp(weekly_growth_rate)
        )
      } else {
        new_infections[week_t] <<- initial_infections *
          exp(weekly_growth_rate * (week_t - 1))
      }
    },
    get_shedding = function(week_t, number_of_weeks_shedding) {
      start <- max(1, week_t - number_of_weeks_shedding + 1)
      sum(new_infections[start:week_t])
    },
    get_cumulative = function(week_t) {
      sum(new_infections[1:week_t])
    }
  )
}

# =============================================================================
# SINGLE SIMULATION FUNCTIONS
# =============================================================================

calc_total_infections_between_weeks <- function(
  week_start,
  week_end,
  initial_infections,
  weekly_growth_rate
) {
  #' Analytically calculate total infections between two weeks (inclusive)
  #' @param week_start Starting week (inclusive)
  #' @param week_end Ending week (inclusive)
  #' @param initial_infections Number of infected people at week 1
  #' @param weekly_growth_rate Growth rate of the pathogen (r)
  #' @return Total number of infections from week_start to week_end

  if (abs(weekly_growth_rate) < 1e-10) {
    # Handle weekly_growth_rate ≈ 0 case to avoid numerical issues
    return(initial_infections * (week_end - week_start + 1))
  }
  # Geometric series sum of infections from week_start to week_end
  n_terms <- week_end - week_start + 1
  first_term <- initial_infections * exp(weekly_growth_rate * (week_start - 1))
  series_sum <- (exp(weekly_growth_rate * n_terms) - 1) /
    (exp(weekly_growth_rate) - 1)
  return(first_term * series_sum)
}

calc_shedding_in_week_t <- function(
  week_t,
  initial_infections,
  weekly_growth_rate,
  number_of_weeks_shedding
) {
  #' Calculate how many people are shedding viral reads at week t
  #' @param week_t Week number during outbreak
  #' @param initial_infections Number of infected people at week 1
  #' @param weekly_growth_rate Growth rate of the pathogen
  #' @param number_of_weeks_shedding Number of weeks an infected person sheds viral reads
  #' @details People infected in weeks [t-shedding_weeks+1, t] are still shedding at week t

  start_week <- max(1, week_t - number_of_weeks_shedding + 1)
  end_week <- week_t

  # Use analytical formula instead of numerical summation
  calc_total_infections_between_weeks(
    start_week,
    end_week,
    initial_infections,
    weekly_growth_rate
  )
}

simulate_single_outbreak_on_weekly_basis <- function(
    weekly_growth_rate,
    initial_infections,
    individuals_sampled,
    sequencing_depth,
    read_detection_threshold,
    mu,
    cv_squared,
    number_of_weeks_shedding,
    total_population_size,
    prop_pop_donating = 1,
    prob_infected_donates = 1,
    stochastic_growth = FALSE,
    max_weeks = 1000
) {
  #' Simulate a single outbreak on a weekly basis
  #' @param weekly_growth_rate Growth rate of the pathogen
  #' @param initial_infections Number of infected people at week 1
  #' @param individuals_sampled Number of people sampled each week
  #' @param sequencing_depth Sequencing depth
  #' @param read_detection_threshold Number of viral reads required to detect the pathogen
  #' @param mu Infected person to relative abundance ratio calculated using the arithmetic mean from untargeted MGS data for this pathogen
  #' @param cv_squared Coefficient of Variation Squared of relative abundance from untargeted MGS data for this pathogen (using the same data as mu)
  #' @param number_of_weeks_shedding Number of weeks an infected person sheds viral reads
  #' @param total_population_size Total population size
  #' @param prop_pop_donating Proportion of the total population donating blood
  #' @param prob_infected_donates Probability that someone who is infected actually donates
  #' @param stochastic_growth If TRUE, use Poisson branching process for infection dynamics; if FALSE, use deterministic exponential growth
  #' @param max_weeks Number of weeks to run simulation for before terminating if the pathogen is not detected
  
  # Calculate effective donor pool and validate sampling feasibility
  effective_donor_pool <- total_population_size * prop_pop_donating
  if (individuals_sampled >= effective_donor_pool) {
    stop(sprintf(
      "individuals_sampled (%s) >= effective donor pool (%s). Reduce individuals_sampled or increase prop_pop_donating.",
      format(individuals_sampled, big.mark = ","),
      format(effective_donor_pool, big.mark = ",")
    ))
  }
  
  # Initialize outbreak tracker
  tracker <- create_outbreak_tracker(
    initial_infections = initial_infections,
    weekly_growth_rate = weekly_growth_rate,
    max_weeks = max_weeks,
    stochastic_growth = stochastic_growth
  )
  
  # Run epidemic week by week until detection
  week_t <- 1
  total_pathogen_reads <- 0
  while (
    week_t <= max_weeks && total_pathogen_reads < read_detection_threshold
  ) {
    tracker$advance_to_week(week_t)
    shedding_population <- tracker$get_shedding(week_t, number_of_weeks_shedding)
    
    # Stop at 1/4 population size since exponential growth slows down around then
    if (shedding_population > total_population_size / 4) {
      break
    }
    
    # Adjust shedding population to reflect donor pool dynamics
    shedding_donors <- shedding_population * prop_pop_donating * prob_infected_donates
    if (shedding_donors > 0) {
      individuals_shedding_viral_reads <- rbinom(
        1,
        individuals_sampled,
        shedding_donors / effective_donor_pool
      )
      if (individuals_shedding_viral_reads > 0) {
        expected_viral_reads <-
          (individuals_shedding_viral_reads / individuals_sampled) *
          mu *
          sequencing_depth
        weekly_viral_reads <- rnbinom(
          1,
          size = individuals_shedding_viral_reads / cv_squared,
          mu = expected_viral_reads
        )
        total_pathogen_reads <- total_pathogen_reads + weekly_viral_reads
        if (total_pathogen_reads >= read_detection_threshold) {
          total_infections <- tracker$get_cumulative(week_t)
          return(list(
            detected = TRUE,
            detection_week = week_t,
            cumulative_incidence = total_infections / total_population_size,
            total_pathogen_reads = total_pathogen_reads
          ))
        }
      }
    }
    week_t <- week_t + 1
  }
  
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
  weekly_growth_rate,
  initial_infections,
  individuals_sampled,
  sequencing_depth,
  read_detection_threshold,
  mu,
  cv_squared,
  number_of_weeks_shedding,
  total_population_size,
  prop_pop_donating = 1,
  prob_infected_donates = 1,
  stochastic_growth = FALSE,
  n_sims
) {
  #' Simulate many outbreaks, then calculate detection probability and summary statistics
  #' @param weekly_growth_rate Growth rate of the pathogen
  #' @param initial_infections Number of infected people at week 1
  #' @param individuals_sampled Number of people sampled each week
  #' @param sequencing_depth Sequencing depth
  #' @param read_detection_threshold Number of viral reads required to detect the pathogen
  #' @param mu Infected person to relative abundance ratio calculated using the arithmetic mean from untargeted MGS data for this pathogen
  #' @param cv_squared Coefficient of Variation Squared of relative abundance from untargeted MGS data for this pathogen (using the same data as mu)
  #' @param number_of_weeks_shedding Number of weeks an infected person sheds viral reads
  #' @param total_population_size Total population size
  #' @param prop_pop_donating Proportion of the total population donating blood (default 1)
  #' @param prob_infected_donates Probability that an infected individual donates, relative to the general donation probability (default 1)
  #' @param stochastic_growth If TRUE, use Poisson branching process for infection dynamics; if FALSE, use deterministic exponential growth (default FALSE)
  #' @param n_sims Number of simulations to run at this sequencing depth

  # Run many simulations
  results <- replicate(
    n_sims,
    {
      simulate_single_outbreak_on_weekly_basis(
        weekly_growth_rate = weekly_growth_rate,
        initial_infections = initial_infections,
        individuals_sampled = individuals_sampled,
        sequencing_depth = sequencing_depth,
        read_detection_threshold = read_detection_threshold,
        mu = mu,
        cv_squared = cv_squared,
        number_of_weeks_shedding = number_of_weeks_shedding,
        total_population_size = total_population_size,
        prop_pop_donating = prop_pop_donating,
        prob_infected_donates = prob_infected_donates,
        stochastic_growth = stochastic_growth
      )
    },
    simplify = FALSE
  )
  sims_with_pathogen_detection <- results[sapply(results, function(x) {
    x$detected
  })]
  # Treat non-detections as 100% incidence for unconditional statistics
  all_incidences <- sapply(results, function(x) {
    ifelse(x$detected, x$cumulative_incidence, 100)
  })
  cumulative_incidences_at_detection <- if (
    length(sims_with_pathogen_detection) > 0
  ) {
    sapply(sims_with_pathogen_detection, function(x) x$cumulative_incidence)
  } else {
    numeric(0)
  }
  summary_stats <- list(
    detection_rate = length(sims_with_pathogen_detection) / n_sims,
    mean_incidence = mean(all_incidences),
    median_incidence = median(all_incidences),
    q25_incidence = quantile(all_incidences, 0.25),
    q75_incidence = quantile(all_incidences, 0.75),
    min_incidence = min(all_incidences),
    max_incidence = max(all_incidences)
  )
  return(list(
    detection_rate = summary_stats$detection_rate,
    summary_stats = summary_stats,
    cumulative_incidences_at_detection = cumulative_incidences_at_detection,
    n_sims = n_sims
  ))
}


calc_prob_detection_before_threshold <- function(
  detection_results,
  threshold_incidence
) {
  #' Calculate probability of detecting pathogen before cumulative incidence exceeds threshold
  #' @param detection_results Output from run_simulations_at_given_sequencing_depth containing detection data and n_sims
  #' @param threshold_incidence Cumulative incidence threshold (as proportion 0-1)
  #' @return Numeric value between 0 and 1 representing the probability of early detection

  if (detection_results$detection_rate == 0) {
    return(0)
  }
  detected_early <- sum(
    detection_results$cumulative_incidences_at_detection <= threshold_incidence
  )
  return(detected_early / detection_results$n_sims)
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================
minimize_sequencing_depth_given_detection_constraint_using_binary_search <- function(
  target_cumulative_incidence,
  target_prob,
  weekly_growth_rate,
  initial_infections,
  individuals_sampled,
  read_detection_threshold,
  mu,
  cv_squared,
  number_of_weeks_shedding,
  total_population_size,
  prop_pop_donating = 1,
  prob_infected_donates = 1,
  stochastic_growth = FALSE,
  n_sims,
  depth_precision = 1e4,
  relative_error_tolerance = 0.025,
  max_attempts = 3
) {
  #' Using binary search, find the minimum sequencing depth that achieves the target
  #' detection probability before the cumulative incidence threshold. Includes retry
  #' logic with relative error tolerance to handle simulation stochasticity. Will retry
  #' up to max_attempts times with increasing simulation counts to ensure convergence.
  #' @param target_cumulative_incidence Cumulative incidence threshold before which we want to detect the pathogen (proportion 0-1)
  #' @param target_prob Minimum probability of detecting the pathogen before cumulative incidence reaches target_cumulative_incidence (default 0.95 for 95% detection probability)
  #' @param weekly_growth_rate Growth rate of the pathogen
  #' @param initial_infections Number of infected people at week 1
  #' @param individuals_sampled Number of people sampled each week
  #' @param read_detection_threshold Number of viral reads required to detect the pathogen
  #' @param mu Infected person to relative abundance ratio calculated using the arithmetic mean from untargeted MGS data for this pathogen
  #' @param cv_squared Coefficient of Variation Squared of relative abundance from untargeted MGS data for this pathogen (using the same data as mu)
  #' @param number_of_weeks_shedding Number of weeks an infected person sheds viral reads
  #' @param total_population_size Total population size
  #' @param prop_pop_donating Proportion of the total population donating blood (default 1)
  #' @param prob_infected_donates Probability that an infected individual donates, relative to the general donation probability (default 1)
  #' @param stochastic_growth If TRUE, use Poisson branching process for infection dynamics; if FALSE, use deterministic exponential growth (default FALSE)
  #' @param n_sims Number of simulations to run at this sequencing depth
  #' @param relative_error_tolerance Relative error tolerance for accepting final probability (default 0.025 for ±2.5%)
  #' @param max_attempts Maximum number of attempts to achieve target probability within tolerance (default 3)

  # Store achieved probabilities for error message if needed
  achieved_probs <- numeric(max_attempts)
  for (attempt in 1:max_attempts) {
    # Increase simulations on retry for better accuracy
    current_n_sims <- n_sims * (2^(attempt - 1))
    # Run binary search
    binary_search_lower_bound <- 1e4 # Start with very low depth
    binary_search_upper_bound <- 1e12 # Maximum reasonable depth
    curr_seq_depth <- NULL # Initialize outside loop to maintain scope
    while (
      binary_search_upper_bound - binary_search_lower_bound > depth_precision
    ) {
      curr_seq_depth <- round(
        (binary_search_lower_bound + binary_search_upper_bound) / 2
      )
      # Run simulations at this depth
      results <- run_simulations_at_given_sequencing_depth(
        weekly_growth_rate = weekly_growth_rate,
        initial_infections = initial_infections,
        individuals_sampled = individuals_sampled,
        sequencing_depth = curr_seq_depth,
        read_detection_threshold = read_detection_threshold,
        mu = mu,
        cv_squared = cv_squared,
        number_of_weeks_shedding = number_of_weeks_shedding,
        total_population_size = total_population_size,
        prop_pop_donating = prop_pop_donating,
        prob_infected_donates = prob_infected_donates,
        stochastic_growth = stochastic_growth,
        n_sims = current_n_sims
      )
      p_detect <- calc_prob_detection_before_threshold(
        results,
        target_cumulative_incidence
      )
      if (p_detect >= target_prob) {
        binary_search_upper_bound <- curr_seq_depth # Can use less sequencing
      } else {
        binary_search_lower_bound <- curr_seq_depth # Need more sequencing
      }
    }
    # Run final simulation at selected depth to confirm
    final_results <- run_simulations_at_given_sequencing_depth(
      weekly_growth_rate = weekly_growth_rate,
      initial_infections = initial_infections,
      individuals_sampled = individuals_sampled,
      sequencing_depth = curr_seq_depth,
      read_detection_threshold = read_detection_threshold,
      mu = mu,
      cv_squared = cv_squared,
      number_of_weeks_shedding = number_of_weeks_shedding,
      total_population_size = total_population_size,
      prop_pop_donating = prop_pop_donating,
      prob_infected_donates = prob_infected_donates,
      stochastic_growth = stochastic_growth,
      n_sims = current_n_sims * 2
    )
    final_p_detect <- calc_prob_detection_before_threshold(
      final_results,
      target_cumulative_incidence
    )
    achieved_probs[attempt] <- final_p_detect
    # Check if within tolerance
    relative_error <- abs(final_p_detect - target_prob) / target_prob
    if (relative_error <= relative_error_tolerance) {
      return(list(
        optimal_depth = curr_seq_depth,
        final_prob = final_p_detect,
        final_results = final_results,
        attempt_number = attempt,
        relative_error = relative_error,
        converged = TRUE,
        achieved_probs = achieved_probs[1:attempt]
      ))
    }
  }
  # Return best result even if not converged
  # Use the last attempt's results as the best we could achieve
  return(list(
    optimal_depth = curr_seq_depth,
    final_prob = final_p_detect,
    final_results = final_results,
    attempt_number = max_attempts,
    relative_error = relative_error,
    converged = FALSE,
    achieved_probs = achieved_probs
  ))
}

# =============================================================================
# COST FUNCTIONS
# =============================================================================

calc_weekly_operational_cost <- function(
  num_individuals_sampled_per_week,
  sequencing_depth_billion_reads_per_week,
  cost_per_individual_sample,
  cost_sample_processing_per_week,
  cost_per_billion_reads_sequencing,
  cost_per_billion_reads_data_processing
) {
  #' Calculate weekly operational costs
  #' @param num_individuals_sampled_per_week Number of individuals sampled each week
  #' @param sequencing_depth_billion_reads_per_week Total sequencing depth per week (in billions of reads)
  #' @param cost_per_individual_sample Cost to acquire one individual sample (unit cost)
  #' @param cost_sample_processing_per_week Fixed weekly cost for sample processing
  #' @param cost_per_billion_reads_sequencing Unit cost per billion reads of sequencing
  #' @param cost_per_billion_reads_data_processing Unit cost per billion reads of data processing

  return(
    (cost_per_individual_sample * num_individuals_sampled_per_week) +
      cost_sample_processing_per_week +
      sequencing_depth_billion_reads_per_week *
        (cost_per_billion_reads_sequencing +
          cost_per_billion_reads_data_processing)
  )
}

calc_storage_unit_cost <- function(
  sequencing_depth_billion_reads_per_week,
  cost_per_billion_reads_storage_per_week
) {
  #' Calculate cost to store one week's data for one week
  #' @param sequencing_depth_billion_reads_per_week Weekly sequencing depth (in billions of reads)
  #' @param cost_per_billion_reads_storage_per_week Unit cost to store one billion reads for one week

  return(
    sequencing_depth_billion_reads_per_week *
      cost_per_billion_reads_storage_per_week
  )
}

calc_cumulative_storage_cost <- function(
  cost_to_store_one_week_data_for_one_week,
  num_weeks_in_period
) {
  #' Calculate cumulative storage cost over the year
  #' @param cost_to_store_one_week_data_for_one_week Cost to store one week's worth of data for one week
  #' @param num_weeks_in_period Number of weeks in the period

  # Sum of arithmetic series
  total_storage_weeks <- num_weeks_in_period * (num_weeks_in_period + 1) / 2
  return(cost_to_store_one_week_data_for_one_week * total_storage_weeks)
}

calc_total_cost_per_year <- function(
  sequencing_depth_billion_reads_per_week,
  cost_per_billion_reads_sequencing,
  cost_sample_processing_per_week,
  cost_per_individual_sample,
  num_individuals_sampled_per_week,
  cost_per_billion_reads_data_processing,
  cost_per_billion_reads_storage_per_week,
  annual_labor_cost
) {
  #' Calculate the total cost per year of detecting this HIV-like pathogen
  #'
  #' @param sequencing_depth_billion_reads_per_week Weekly sequencing depth (in billions of reads)
  #' @param cost_per_billion_reads_sequencing Unit cost per billion reads of sequencing
  #' @param cost_sample_processing_per_week Fixed weekly cost for sample processing
  #' @param cost_per_individual_sample Unit cost to acquire one individual sample
  #' @param num_individuals_sampled_per_week Number of individuals sampled each week
  #' @param cost_per_billion_reads_data_processing Unit cost per billion reads of data processing
  #' @param cost_per_billion_reads_storage_per_week Unit cost to store one billion reads for one week
  #' @param annual_labor_cost Total annual labor cost

  cost_weekly_ops <- calc_weekly_operational_cost(
    num_individuals_sampled_per_week = num_individuals_sampled_per_week,
    sequencing_depth_billion_reads_per_week = sequencing_depth_billion_reads_per_week,
    cost_per_individual_sample = cost_per_individual_sample,
    cost_sample_processing_per_week = cost_sample_processing_per_week,
    cost_per_billion_reads_sequencing = cost_per_billion_reads_sequencing,
    cost_per_billion_reads_data_processing = cost_per_billion_reads_data_processing
  )
  cost_storage_unit <- calc_storage_unit_cost(
    sequencing_depth_billion_reads_per_week = sequencing_depth_billion_reads_per_week,
    cost_per_billion_reads_storage_per_week = cost_per_billion_reads_storage_per_week
  )
  cumulative_storage <- calc_cumulative_storage_cost(
    cost_to_store_one_week_data_for_one_week = cost_storage_unit,
    num_weeks_in_period = 52
  )
  cost_annual <- (cost_weekly_ops * 52) + cumulative_storage + annual_labor_cost
  return(cost_annual)
}
