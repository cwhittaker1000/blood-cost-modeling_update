# =============================================================================
# COMBINED COMPARISON: simple exponential vs deterministic renewal (national)
# Colour    = individuals_sampled
# Linetype  = growth model (simple vs renewal)
# Shape     = growth model
# =============================================================================

library(readr); library(dplyr); library(ggplot2); library(scales); library(cowplot)
source("scripts/functions.R")

hiv_mu <- 1.26e-5

# -----------------------------------------------------------------------------
# Load and tag both result sets
# -----------------------------------------------------------------------------
simple_results <- read_tsv(
  "results/summarized_simulation_results.tsv",
  show_col_types = FALSE
) %>%
  mutate(growth_model = "Simple exponential")

renewal_results <- read_tsv(
  "results/summarized_simulation_results_renewal.tsv",
  show_col_types = FALSE
) %>%
  mutate(growth_model = "Deterministic renewal")

# Bind on common columns so we're robust to either dataset having extras
common_cols <- intersect(names(simple_results), names(renewal_results))

combined <- bind_rows(
  simple_results[, common_cols],
  renewal_results[, common_cols]
) %>%
  filter(converged == TRUE) %>%
  mutate(
    individuals_sampled = factor(
      individuals_sampled,
      levels = sort(unique(individuals_sampled))
    ),
    growth_model = factor(growth_model,
                          levels = c("Simple exponential",
                                     "Deterministic renewal"))
  )

# =============================================================================
# COST PARAMETERS
# =============================================================================
cost_per_individual_sample              <- 5.00
cost_sample_processing_per_week         <- 500
cost_per_billion_reads_sequencing       <- 2500
cost_per_billion_reads_data_processing  <- 1
cost_per_billion_reads_storage_per_week <- 0.25
annual_labor_cost                       <- 650000

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

# Shared scales for linetype/shape so all three panels stay consistent
model_linetypes <- c("Simple exponential"    = "solid",
                     "Deterministic renewal" = "dashed")
model_shapes    <- c("Simple exponential"    = 16,
                     "Deterministic renewal" = 17)

# =============================================================================
# PANEL A1: weekly sequencing depth vs cumulative incidence (HIV mu)
# =============================================================================
a1_data <- combined %>%
  filter(abs(mu - hiv_mu) < 1e-7)

p_a1 <- ggplot(
  a1_data,
  aes(
    x        = target_cumulative_incidence,
    y        = optimal_depth,
    colour   = individuals_sampled,
    linetype = growth_model,
    shape    = growth_model,
    group    = interaction(individuals_sampled, growth_model)
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
  scale_linetype_manual(values = model_linetypes) +
  scale_shape_manual(values = model_shapes) +
  labs(
    x        = "Cumulative incidence",
    y        = "Weekly sequencing depth",
    colour   = "# of individuals",
    linetype = "Growth model",
    shape    = "Growth model",
    title    = "a1"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "right",
    legend.background  = element_rect(fill = "white", colour = "grey80"),
    plot.title         = element_text(face = "bold", hjust = 0),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

p_a1

# =============================================================================
# PANEL A2: annual cost vs cumulative incidence (HIV mu)
# =============================================================================
a2_data <- combined %>%
  filter(abs(mu - hiv_mu) < 1e-7)

p_a2 <- ggplot(
  a2_data,
  aes(
    x        = target_cumulative_incidence,
    y        = annual_cost,
    colour   = individuals_sampled,
    linetype = growth_model,
    shape    = growth_model,
    group    = interaction(individuals_sampled, growth_model)
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
  scale_linetype_manual(values = model_linetypes) +
  scale_shape_manual(values = model_shapes) +
  labs(
    x        = "Cumulative incidence",
    y        = "Annual cost ($)",
    colour   = "# of individuals",
    linetype = "Growth model",
    shape    = "Growth model",
    title    = "a2"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", hjust = 0),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

p_a2

# =============================================================================
# PANEL B: annual cost vs mu, fixed cumulative incidence = 0.01% (1e-4)
# =============================================================================
b_ci_target <- 1e-4

b_data <- combined %>%
  filter(target_cumulative_incidence == b_ci_target) %>%
  filter(individuals_sampled != 2e+03)

p_b <- ggplot(
  b_data,
  aes(
    x        = mu,
    y        = annual_cost,
    colour   = individuals_sampled,
    linetype = growth_model,
    shape    = growth_model,
    group    = interaction(individuals_sampled, growth_model)
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_linetype_manual(values = model_linetypes) +
  scale_shape_manual(values = model_shapes) +
  scale_x_log10(
    labels = label_log(base = 10),
    breaks = 10^(-8:-3)
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(6:9)
  ) +
  labs(
    x        = "RA of target pathogen",
    y        = sprintf("Annual cost ($), CI = %s%%", signif(b_ci_target * 100, 2)),
    colour   = "# of individuals",
    linetype = "Growth model",
    shape    = "Growth model",
    title    = "b"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "right",
    legend.background  = element_rect(fill = "white", colour = "grey80"),
    legend.box         = "vertical",
    plot.title         = element_text(face = "bold", hjust = 0),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  ) +
  guides(colour = "none")  # colour legend already shown in a1

p_b

# =============================================================================
# Compose
# =============================================================================
cowplot::plot_grid(p_a1, p_a2, p_b, nrow = 1)