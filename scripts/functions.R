library(docstring)

# =============================================================================
# HELPER FUNCTION
# =============================================================================

calc_R0_from_doubling_time <- function(doubling_time, infectious_weeks) {
  #' Calculate R0 from doubling time and duration of infectiousness
  #' using the Euler-Lotka equation with a uniform generation interval.
  #' @param doubling_time Doubling time of the epidemic in weeks
  #' @param infectious_weeks Duration of infectiousness in weeks (D)
  #' @return R0: basic reproduction number
  r <- log(2) / doubling_time
  R0 <- infectious_weeks / sum(exp(-r * (1:infectious_weeks)))
  return(R0)
}

calc_r_from_R0 <- function(R0, infectious_weeks, tol = 1e-10) {
  #' Calculate the exponential growth rate r from R0 and duration of
  #' infectiousness by numerically solving the Euler-Lotka equation.
  #' @param R0 Basic reproduction number
  #' @param infectious_weeks Duration of infectiousness in weeks (D)
  #' @param tol Numerical tolerance for root finding
  #' @return r: exponential growth rate per week
  f <- function(r) {
    R0 - infectious_weeks / sum(exp(-r * (1:infectious_weeks)))
  }
  result <- uniroot(f, interval = c(0, 2), tol = tol)
  return(result$root)
}

# =============================================================================
# OUTBREAK TRACKER
# =============================================================================

create_outbreak_tracker <- function(
    initial_infections,
    R0,
    max_weeks,
    growth_model = "deterministic_renewal",
    infectious_weeks = 12
) {
  #' Create an outbreak tracker that manages infection dynamics.
  #' Three growth models are available:
  #'   "deterministic_simple" - Simple exponential growth: I_new(t) = I0 * exp(r*(t-1)),
  #'     where r is derived from R0 and infectious_weeks via the Euler-Lotka equation.
  #'     No assumptions about infectiousness duration in the growth dynamics themselves.
  #'   "deterministic_renewal" - Deterministic renewal equation with finite infectiousness:
  #'     I_new(t) = beta * sum(I_new(t-k) for k in 1:infectious_weeks), where
  #'     beta = R0 / infectious_weeks.
  #'   "stochastic_renewal" - Poisson renewal process: I_new(t) ~ Poisson(beta * infectious_pop).
  #'     Mathematically equivalent to an age-dependent (Bellman-Harris) branching process.
  #'     Matches deterministic_renewal in expectation.
  #' @param initial_infections Number of infected people at week 1
  #' @param R0 Basic reproduction number
  #' @param max_weeks Maximum number of weeks to track
  #' @param growth_model One of "deterministic_simple", "deterministic_renewal",
  #'   or "stochastic_renewal" (default "deterministic_renewal")
  #' @param infectious_weeks Number of weeks an infected person remains infectious (D).
  #'   Used by renewal models for transmission dynamics and by deterministic_simple
  #'   to derive r from R0 via the Euler-Lotka equation. (default 12)
  #' @return List of functions: advance_to_week, get_shedding, get_cumulative
  
  valid_models <- c("deterministic_simple", "deterministic_renewal", "stochastic_renewal")
  if (!growth_model %in% valid_models) {
    stop(sprintf(
      "growth_model must be one of: %s. Got '%s'.",
      paste(valid_models, collapse = ", "), growth_model
    ))
  }
  
  new_infections <- numeric(max_weeks)
  new_infections[1] <- initial_infections
  
  # Per-week transmission rate for renewal models
  beta <- R0 / infectious_weeks
  
  # Growth rate for simple exponential model, derived from R0 and D
  # via the Euler-Lotka equation so that all models are parameterised
  # consistently from the same (R0, D) pair
  r <- if (growth_model == "deterministic_simple") {
    calc_r_from_R0(R0, infectious_weeks)
  } else {
    NULL
  }
  
  list(
    advance_to_week = function(week_t) {
      if (week_t == 1) return(invisible(NULL))
      if (new_infections[week_t] != 0) return(invisible(NULL))
      
      if (growth_model == "deterministic_simple") {
        new_infections[week_t] <<- initial_infections *
          exp(r * (week_t - 1))
        
      } else if (growth_model == "deterministic_renewal") {
        infectious_start <- max(1, week_t - infectious_weeks)
        infectious_pop <- sum(new_infections[infectious_start:(week_t - 1)])
        new_infections[week_t] <<- beta * infectious_pop
        
      } else if (growth_model == "stochastic_renewal") {
        infectious_start <- max(1, week_t - infectious_weeks)
        infectious_pop <- sum(new_infections[infectious_start:(week_t - 1)])
        new_infections[week_t] <<- rpois(1, beta * infectious_pop)
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
# SINGLE SIMULATION FUNCTION
# =============================================================================

simulate_single_outbreak_on_weekly_basis <- function(
    R0,
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
    growth_model = "deterministic_renewal",
    max_weeks = 1000,
    max_extinction_retries = 100
) {
  #' Simulate a single outbreak on a weekly basis
  #' @param R0 Basic reproduction number
  #' @param initial_infections Number of infected people at week 1
  #' @param individuals_sampled Number of people sampled each week
  #' @param sequencing_depth Sequencing depth
  #' @param read_detection_threshold Number of viral reads required to detect the pathogen
  #' @param mu Infected person to relative abundance ratio
  #' @param cv_squared Coefficient of Variation Squared of relative abundance
  #' @param number_of_weeks_shedding Number of weeks an infected person sheds viral reads
  #' @param total_population_size Total population size
  #' @param prop_pop_donating Proportion of the total population donating blood (default 1)
  #' @param prob_infected_donates Probability that an infected individual donates (default 1)
  #' @param growth_model One of "deterministic_simple", "deterministic_renewal",
  #'   or "stochastic_renewal" (default "deterministic_renewal")
  #' @param max_weeks Number of weeks to run simulation before terminating (default 1000)
  #' @param max_extinction_retries Maximum retries if stochastic outbreak goes
  #'   extinct (default 100). Only applies when growth_model = "stochastic_renewal".
  
  # Calculate effective donor pool and validate sampling feasibility
  effective_donor_pool <- total_population_size * prop_pop_donating
  if (individuals_sampled >= effective_donor_pool) {
    stop(sprintf(
      "individuals_sampled (%s) >= effective donor pool (%s). Reduce individuals_sampled or increase prop_pop_donating.",
      format(individuals_sampled, big.mark = ","),
      format(effective_donor_pool, big.mark = ",")
    ))
  }
  
  retries <- 0
  repeat {
    # Initialize outbreak tracker
    tracker <- create_outbreak_tracker(
      initial_infections = initial_infections,
      R0 = R0,
      max_weeks = max_weeks,
      growth_model = growth_model,
      infectious_weeks = number_of_weeks_shedding
    )
    
    # Run epidemic week by week until detection
    week_t <- 1
    total_pathogen_reads <- 0
    outbreak_extinct <- FALSE
    
    while (
      week_t <= max_weeks && total_pathogen_reads < read_detection_threshold
    ) {
      tracker$advance_to_week(week_t)
      shedding_population <- tracker$get_shedding(week_t, number_of_weeks_shedding)
      
      # Check for stochastic extinction: if no one is shedding after week 1,
      # the outbreak has died out and cannot recover
      if (growth_model == "stochastic_renewal" && week_t > 1 && shedding_population == 0) {
        outbreak_extinct <- TRUE
        break
      }
      
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
              total_pathogen_reads = total_pathogen_reads,
              extinction_retries = retries
            ))
          }
        }
      }
      week_t <- week_t + 1
    }
    
    # If outbreak went extinct and we're in stochastic mode, retry
    if (outbreak_extinct && retries < max_extinction_retries) {
      retries <- retries + 1
      next
    }
    
    # If we hit the retry cap, warn
    if (outbreak_extinct && retries >= max_extinction_retries) {
      warning(sprintf(
        "Outbreak went extinct %d times in a row. Returning non-detection.",
        max_extinction_retries
      ))
    }
    
    # No detection
    return(list(
      detected = FALSE,
      detection_week = NA,
      cumulative_incidence = NA,
      total_pathogen_reads = total_pathogen_reads,
      extinction_retries = retries
    ))
  }
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

run_simulations_at_given_sequencing_depth <- function(
    R0,
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
    growth_model = "deterministic_renewal",
    n_sims
) {
  #' Simulate many outbreaks, then calculate detection probability and summary statistics
  #' @param R0 Basic reproduction number
  #' @param initial_infections Number of infected people at week 1
  #' @param individuals_sampled Number of people sampled each week
  #' @param sequencing_depth Sequencing depth
  #' @param read_detection_threshold Number of viral reads required to detect the pathogen
  #' @param mu Infected person to relative abundance ratio
  #' @param cv_squared Coefficient of Variation Squared of relative abundance
  #' @param number_of_weeks_shedding Number of weeks an infected person sheds viral reads
  #' @param total_population_size Total population size
  #' @param prop_pop_donating Proportion of the total population donating blood (default 1)
  #' @param prob_infected_donates Probability that an infected individual donates (default 1)
  #' @param growth_model One of "deterministic_simple", "deterministic_renewal",
  #'   or "stochastic_renewal" (default "deterministic_renewal")
  #' @param n_sims Number of simulations to run at this sequencing depth
  
  results <- replicate(
    n_sims,
    {
      simulate_single_outbreak_on_weekly_basis(
        R0 = R0,
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
        growth_model = growth_model
      )
    },
    simplify = FALSE
  )
  sims_with_pathogen_detection <- results[sapply(results, function(x) {
    x$detected
  })]
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
  #' @param detection_results Output from run_simulations_at_given_sequencing_depth
  #' @param threshold_incidence Cumulative incidence threshold (as proportion 0-1)
  #' @return Numeric value between 0 and 1
  
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
    R0,
    initial_infections,
    individuals_sampled,
    read_detection_threshold,
    mu,
    cv_squared,
    number_of_weeks_shedding,
    total_population_size,
    prop_pop_donating = 1,
    prob_infected_donates = 1,
    growth_model = "deterministic_renewal",
    n_sims,
    depth_precision = 1e4,
    relative_error_tolerance = 0.025,
    max_attempts = 3
) {
  #' Using binary search, find the minimum sequencing depth that achieves the target
  #' detection probability before the cumulative incidence threshold.
  #' @param target_cumulative_incidence Cumulative incidence threshold (proportion 0-1)
  #' @param target_prob Minimum detection probability before threshold
  #' @param R0 Basic reproduction number
  #' @param initial_infections Number of infected people at week 1
  #' @param individuals_sampled Number of people sampled each week
  #' @param read_detection_threshold Number of viral reads required to detect the pathogen
  #' @param mu Infected person to relative abundance ratio
  #' @param cv_squared Coefficient of Variation Squared of relative abundance
  #' @param number_of_weeks_shedding Number of weeks an infected person sheds viral reads
  #' @param total_population_size Total population size
  #' @param prop_pop_donating Proportion of the total population donating blood (default 1)
  #' @param prob_infected_donates Probability that an infected individual donates (default 1)
  #' @param growth_model One of "deterministic_simple", "deterministic_renewal",
  #'   or "stochastic_renewal" (default "deterministic_renewal")
  #' @param n_sims Number of simulations to run at each sequencing depth
  #' @param depth_precision Precision of binary search (default 1e4)
  #' @param relative_error_tolerance Relative error tolerance for final probability (default 0.025)
  #' @param max_attempts Maximum attempts to converge (default 3)
  
  achieved_probs <- numeric(max_attempts)
  for (attempt in 1:max_attempts) {
    current_n_sims <- n_sims * (2^(attempt - 1))
    binary_search_lower_bound <- 1e4
    binary_search_upper_bound <- 1e12
    curr_seq_depth <- NULL
    while (
      binary_search_upper_bound - binary_search_lower_bound > depth_precision
    ) {
      curr_seq_depth <- round(
        (binary_search_lower_bound + binary_search_upper_bound) / 2
      )
      results <- run_simulations_at_given_sequencing_depth(
        R0 = R0,
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
        growth_model = growth_model,
        n_sims = current_n_sims
      )
      p_detect <- calc_prob_detection_before_threshold(
        results,
        target_cumulative_incidence
      )
      if (p_detect >= target_prob) {
        binary_search_upper_bound <- curr_seq_depth
      } else {
        binary_search_lower_bound <- curr_seq_depth
      }
    }
    final_results <- run_simulations_at_given_sequencing_depth(
      R0 = R0,
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
      growth_model = growth_model,
      n_sims = current_n_sims * 2
    )
    final_p_detect <- calc_prob_detection_before_threshold(
      final_results,
      target_cumulative_incidence
    )
    achieved_probs[attempt] <- final_p_detect
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
  #' Calculate the total cost per year
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