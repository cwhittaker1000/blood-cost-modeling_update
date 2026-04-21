## Renewal-model version of reproducing_old_figures.R
## Reads the renewal-model results and produces the same a1/a2/b panels,
## plus a sensitivity panel on prob_infected_donates.

library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)

source("scripts/functions.R")

# Load the results produced by run_analysis_renewal.R
opt_seq_depth_summary <- read_tsv(
  "results/summarized_simulation_results_renewal_probdonates.tsv",
  show_col_types = FALSE
)

hiv_mu <- 1.26e-5

# =============================================================================
# STANDARDISED COLOUR PALETTE FOR individuals_sampled
# =============================================================================
# Define one palette keyed by individuals_sampled value, used across a1/a2/b
# so the same n always maps to the same colour. Built from scales::hue_pal()
# (ggplot's default) so we don't change the existing colour scheme — we just
# anchor it.
n_levels  <- as.character(sort(unique(opt_seq_depth_summary$individuals_sampled)))
n_palette <- scales::hue_pal()(length(n_levels))
names(n_palette) <- n_levels

# =============================================================================
# PANEL A1: weekly sequencing depth vs cumulative incidence
# =============================================================================
a1_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7) %>%
  filter(prob_infected_donates == 1) %>%
  filter(converged == TRUE) %>%
  mutate(individuals_sampled = factor(
    individuals_sampled,
    levels = n_levels
  ))

p_a1 <- ggplot(
  a1_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = individuals_sampled,
    group = individuals_sampled
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
    breaks = 10^(9:13)
  ) +
  scale_colour_manual(values = n_palette) +
  labs(
    x = "Cumulative incidence",
    y = "Weekly sequencing depth",
    colour = "# of individuals"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.75, 0.78),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_a1

# =============================================================================
# PANEL A3: prob_infected_donates sensitivity at fixed n = 50,000
# Matched colour gradient: grey (at 0.25) -> n=50,000's colour (at 1.00)
# =============================================================================
fixed_n <- 50000

a3_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7) %>%
  filter(individuals_sampled == fixed_n) %>%
  filter(converged == TRUE) %>%
  mutate(prob_infected_donates = factor(
    prob_infected_donates,
    levels = sort(unique(prob_infected_donates))
  ))

# Build discrete palette: grey at lowest prob, n=50000's colour at highest
prob_levels      <- levels(a3_data$prob_infected_donates)
n_prob_levels    <- length(prob_levels)
fixed_n_colour   <- n_palette[as.character(fixed_n)]
prob_palette     <- colorRampPalette(c("grey85", fixed_n_colour))(n_prob_levels)
names(prob_palette) <- prob_levels

p_a3_line2 <- ggplot(
  a3_data,
  aes(
    x      = target_cumulative_incidence,
    y      = optimal_depth,
    colour = prob_infected_donates,
    group  = prob_infected_donates
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
    breaks = 10^(8:13)
  ) +
  scale_colour_manual(
    values = prob_palette,
    guide  = guide_legend(
      override.aes = list(size = 4, shape = 15)  # filled squares
    )
  ) +
  labs(
    x      = "Cumulative incidence",
    y      = "Weekly sequencing depth",
    colour = "P(infected\ndonates)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position   = c(0.85, 0.72),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    plot.title        = element_text(face = "bold", hjust = 0),
    axis.text.x       = element_text(angle = 45, hjust = 1)
  )

p_a3_line2

# =============================================================================
# COST PARAMETERS
# =============================================================================
cost_per_individual_sample              <- 5.00
cost_sample_processing_per_week         <- 500
cost_per_billion_reads_sequencing       <- 2500
cost_per_billion_reads_data_processing  <- 1
cost_per_billion_reads_storage_per_week <- 0.25
annual_labor_cost                       <- 650000

# =============================================================================
# ADD ANNUAL COST TO THE SUMMARY TABLE
# =============================================================================
opt_seq_depth_summary <- opt_seq_depth_summary %>%
  mutate(
    annual_cost = calc_total_cost_per_year(
      sequencing_depth_billion_reads_per_week   = optimal_depth / 1e9,
      cost_per_billion_reads_sequencing         = cost_per_billion_reads_sequencing,
      cost_sample_processing_per_week           = cost_sample_processing_per_week,
      cost_per_individual_sample                = cost_per_individual_sample,
      num_individuals_sampled_per_week          = individuals_sampled,
      cost_per_billion_reads_data_processing    = cost_per_billion_reads_data_processing,
      cost_per_billion_reads_storage_per_week   = cost_per_billion_reads_storage_per_week,
      annual_labor_cost                         = annual_labor_cost
    )
  )

# =============================================================================
# PANEL A2: annual cost vs cumulative incidence, filtered to HIV mu
# =============================================================================
a2_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7) %>%
  filter(prob_infected_donates == 1) %>%
  filter(converged == TRUE) %>%
  mutate(individuals_sampled = factor(
    individuals_sampled,
    levels = n_levels
  ))

p_a2 <- ggplot(
  a2_data,
  aes(
    x = target_cumulative_incidence,
    y = annual_cost,
    colour = individuals_sampled,
    group = individuals_sampled
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
    breaks = 10^(6:9)
  ) +
  scale_colour_manual(values = n_palette) +
  labs(
    x = "Cumulative incidence",
    y = "Annual cost ($)",
    colour = "# of individuals"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_a2

# =============================================================================
# PANEL B: annual cost vs mu, filtered to cumulative incidence = 0.01% (1e-4)
# =============================================================================
b_ci_target <- 1e-4

b_data <- opt_seq_depth_summary %>%
  filter(target_cumulative_incidence == b_ci_target) %>%
  filter(prob_infected_donates == 1) %>%
  filter(converged == TRUE) %>%
  mutate(
    parameterization = ifelse(
      abs(mu - hiv_mu) < 1e-7, "HIV estimate", "Simulated"
    ),
    individuals_sampled = factor(
      individuals_sampled,
      levels = n_levels
    )
  ) %>%
  filter(individuals_sampled != 2e+03)

p_b <- ggplot(
  b_data,
  aes(
    x = mu,
    y = annual_cost,
    colour = individuals_sampled,
    group = individuals_sampled
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(shape = parameterization), size = 2.5) +
  scale_shape_manual(values = c("HIV estimate" = 16, "Simulated" = 15)) +
  scale_colour_manual(values = n_palette, drop = FALSE) +
  scale_x_log10(
    labels = label_log(base = 10),
    breaks = 10^(-8:-3)
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(6:9)
  ) +
  labs(
    x = "RA of target pathogen",
    y = sprintf("Annual cost ($), CI = %s%%", signif(b_ci_target * 100, 2)),
    colour = "# of individuals",
    shape = "Parameterization"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.78, 0.78),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    legend.box = "vertical",
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(colour = "none")

p_b

cowplot::plot_grid(p_a1, p_a3_line2,
                   p_a2, p_b, nrow = 2,
                   labels = c("a", "b", "c", "d"))