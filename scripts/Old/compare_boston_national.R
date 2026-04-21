library(readr); library(dplyr); library(ggplot2); library(scales)

hiv_mu <- 1.26e-5

# -----------------------------------------------------------------------------
# Load and tag both result sets
# -----------------------------------------------------------------------------
national_results <- read_tsv(
  "results/summarized_simulation_results.tsv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "National (uniform spread, 340M)")

boston_results <- read_tsv(
  "results/boston_summarized_simulation_results.tsv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "Boston (concentrated, 5M)")

# -----------------------------------------------------------------------------
# Filter and rescale CI to a common US-population basis
# -----------------------------------------------------------------------------
# Boston's internal CI is computed against 5M. Re-express it as a fraction of
# the US population (340M) so the two scenarios share a common x-axis.
# This represents a "Boston-first emergence" framing: the absolute infection
# counts are real, the population denominator is the full US.
boston_filtered <- boston_results %>%
  filter(prob_infected_donates == 1) %>%
  mutate(
    target_cumulative_incidence_us = target_cumulative_incidence * (5e6 / 3.4e8)
  )

national_filtered <- national_results %>%
  mutate(
    prob_infected_donates = 1,
    target_cumulative_incidence_us = target_cumulative_incidence
  )

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
    scenario = factor(scenario, levels = c("National (uniform spread, 340M)",
                                           "Boston (concentrated, 5M)"))
  )

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
p_combined <- ggplot(
  combined,
  aes(
    x        = target_cumulative_incidence_us,
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
    labels = function(x) {
      ifelse(x >= 1e-4,
             paste0(signif(x * 100, 2), "%"),
             sprintf("%.0e", x))
    }
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(7:13)
  ) +
  scale_linetype_manual(values = c("National (uniform spread, 340M)" = "solid",
                                   "Boston (concentrated, 5M)"        = "dashed")) +
  scale_shape_manual(values = c("National (uniform spread, 340M)" = 16,
                                "Boston (concentrated, 5M)"        = 17)) +
  labs(
    x        = "Cumulative incidence (% of US population)",
    y        = "Weekly sequencing depth",
    colour   = "# sampled / wk",
    linetype = "Scenario",
    shape    = "Scenario",
    title    = "Detection threshold under uniform vs concentrated outbreak",
    subtitle = "Boston scenario assumes all infections are contained in a 5M metro;\nCI re-expressed against 340M US population for direct comparison"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "right",
    legend.background  = element_rect(fill = "white", colour = "grey80"),
    plot.title         = element_text(face = "bold", hjust = 0),
    plot.subtitle      = element_text(size = 10, colour = "grey30"),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

print(p_combined)