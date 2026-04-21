# =============================================================================
# COMBINED COMPARISON: national + Boston on a single plot
# Colour = individuals_sampled (absolute count, shared across datasets)
# Linetype = scenario (national vs Boston)
# =============================================================================

library(readr); library(dplyr); library(ggplot2); library(scales)

hiv_mu <- 1.26e-5

# -----------------------------------------------------------------------------
# Load and tag both result sets
# -----------------------------------------------------------------------------
national_results <- read_tsv(
  "results/summarized_simulation_results.tsv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "National (US, ~340M)")

boston_results <- read_tsv(
  "results/boston_summarized_simulation_results.tsv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "Boston MSA (~5M)")

# Boston has prob_infected_donates; national doesn't. Filter Boston to == 1
# so we're comparing apples to apples, and add the column to national so the
# bind works.
boston_filtered <- boston_results %>%
  filter(prob_infected_donates == 1)

national_filtered <- national_results %>%
  mutate(prob_infected_donates = 1)

# Common columns only, then bind
common_cols <- intersect(names(national_filtered), names(boston_filtered))

combined <- bind_rows(
  national_filtered[, common_cols],
  boston_filtered[, common_cols]
) %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         converged == TRUE) %>%
  mutate(
    individuals_sampled = factor(
      individuals_sampled,
      levels = sort(unique(individuals_sampled))
    ),
    scenario = factor(scenario, levels = c("National (US, ~340M)",
                                           "Boston MSA (~5M)"))
  )

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
p_combined <- ggplot(
  combined,
  aes(
    x        = target_cumulative_incidence,
    y        = optimal_depth,
    colour   = individuals_sampled,
    linetype = scenario,
    shape    = scenario,
    group    = interaction(individuals_sampled, scenario)
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_x_log10(
    labels = function(x) paste0(signif(x * 100, 2), "%"),
    breaks = c(1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3)
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(7:13)
  ) +
  scale_linetype_manual(values = c("National (US, ~340M)" = "solid",
                                   "Boston MSA (~5M)"      = "dashed")) +
  scale_shape_manual(values = c("National (US, ~340M)" = 16,
                                "Boston MSA (~5M)"      = 17)) +
  labs(
    x        = "Cumulative incidence",
    y        = "Weekly sequencing depth",
    colour   = "# sampled / wk",
    linetype = "Scenario",
    shape    = "Scenario",
    title    = "National vs Boston: required sequencing depth at matched sample size"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "right",
    legend.background  = element_rect(fill = "white", colour = "grey80"),
    plot.title         = element_text(face = "bold", hjust = 0),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

print(p_combined)

# =============================================================================
# COST PARAMETERS
# =============================================================================
# See supplementary text for sources. Central estimates used throughout.
cost_per_individual_sample              <- 5.00      # $ per sample (CDC/Vitalant/ARC, central)
cost_sample_processing_per_week         <- 500       # $ per week, fixed (MIT BioMicroCenter)
cost_per_billion_reads_sequencing       <- 2500      # $ per 1B read pairs (Illumina NovaSeq X, Grimm et al. 2025)
cost_per_billion_reads_data_processing  <- 1         # $ per 1B read pairs (NAO)
cost_per_billion_reads_storage_per_week <- 0.25      # $ per 1B read pairs per week (NAO)
annual_labor_cost                       <- 650000    # $ per year (NAO: lab manager + tech + 2 bioinformaticians)

# Requires source("scripts/functions.R") earlier in the session for calc_total_cost_per_year
combined <- combined %>%
  mutate(
    annual_cost = calc_total_cost_per_year(
      sequencing_depth_billion_reads_per_week   = optimal_depth / 1e9,
      cost_per_billion_reads_sequencing         = cost_per_billion_reads_sequencing,
      cost_sample_processing_per_week           = cost_sample_processing_per_week,
      cost_per_individual_sample                = cost_per_individual_sample,
      num_individuals_sampled_per_week          = as.numeric(as.character(individuals_sampled)),
      cost_per_billion_reads_data_processing    = cost_per_billion_reads_data_processing,
      cost_per_billion_reads_storage_per_week   = cost_per_billion_reads_storage_per_week,
      annual_labor_cost                         = annual_labor_cost
    )
  )

# =============================================================================
# Companion plot: annual cost vs cumulative incidence
# =============================================================================
p_combined_cost <- ggplot(
  combined,
  aes(
    x        = target_cumulative_incidence,
    y        = annual_cost,
    colour   = individuals_sampled,
    linetype = scenario,
    shape    = scenario,
    group    = interaction(individuals_sampled, scenario)
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_x_log10(
    labels = function(x) paste0(signif(x * 100, 2), "%"),
    breaks = c(1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3)
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(5:9)
  ) +
  scale_linetype_manual(values = c("National (US, ~340M)" = "solid",
                                   "Boston MSA (~5M)"      = "dashed")) +
  scale_shape_manual(values = c("National (US, ~340M)" = 16,
                                "Boston MSA (~5M)"      = 17)) +
  labs(
    x        = "Cumulative incidence",
    y        = "Annual cost ($)",
    colour   = "# sampled / wk",
    linetype = "Scenario",
    shape    = "Scenario",
    title    = "National vs Boston: annual cost at matched sample size"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "right",
    legend.background  = element_rect(fill = "white", colour = "grey80"),
    plot.title         = element_text(face = "bold", hjust = 0),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

print(p_combined_cost)
