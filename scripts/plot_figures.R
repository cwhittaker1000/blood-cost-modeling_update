# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)

# Sourcing required functions
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
# Image-matched palette inferred from the reference figure.
image_matched_palette <- c(
  "#4077AB", # 2,000
  "#EF6677", # 10,000
  "#1F8942", # 50,000
  "#CCBB44"  # 100,000
)

n_levels <- as.character(sort(unique(opt_seq_depth_summary$individuals_sampled)))

if (length(n_levels) != 4) {
  stop(
    sprintf(
      "Expected 4 individuals_sampled levels, but found %d: %s",
      length(n_levels),
      paste(n_levels, collapse = ", ")
    )
  )
}

n_palette <- image_matched_palette
names(n_palette) <- n_levels
n_labels <- comma(as.numeric(n_levels))

# =============================================================================
# SHARED FORMATTING HELPERS
# =============================================================================
ci_breaks <- c(1e-5, 5e-5, 1e-4, 5e-4, 1e-3)
ci_labels <- percent(ci_breaks, accuracy = 0.001)

ref_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
      panel.grid.minor = element_line(colour = "grey95", linewidth = 0.25),
      axis.line = element_line(colour = "grey20", linewidth = 0.6),
      axis.ticks = element_line(colour = "grey20", linewidth = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, colour = "grey20"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold", hjust = 0, size = 20,
                                margin = margin(b = -2)),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = alpha("white", 0.95),
                                       colour = "grey70", linewidth = 0.35),
      legend.key = element_blank(),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
}

# =============================================================================
# PANEL A1: weekly sequencing depth vs cumulative incidence
# =============================================================================
a1_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7) %>%
  filter(prob_infected_donates == 1) %>%
  filter(converged == TRUE) %>%
  mutate(individuals_sampled = factor(individuals_sampled, levels = n_levels))

p_a1 <- ggplot(
  a1_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = individuals_sampled,
    group = individuals_sampled
  )
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_x_log10(
    breaks = ci_breaks,
    labels = ci_labels,
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(9:12),
    limits = c(7e8, 1.3e12),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = n_palette,
    labels = n_labels,
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(linewidth = 1.0, size = 3.0)
    )
  ) +
  labs(
    title = "a",
    x = "Cumulative incidence",
    y = "Weekly sequencing depth",
    colour = "# of individuals"
  ) +
  ref_theme() +
  theme(
    legend.position = c(0.73, 0.80)
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
  mutate(individuals_sampled = factor(individuals_sampled, levels = n_levels))

x <- a2_data %>%
  filter(target_cumulative_incidence == 1e-04,
         individuals_sampled == 10000)

p_a2 <- ggplot(
  a2_data,
  aes(
    x = target_cumulative_incidence,
    y = annual_cost,
    colour = individuals_sampled,
    group = individuals_sampled
  )
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_x_log10(
    breaks = ci_breaks,
    labels = ci_labels,
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(6:8),
    limits = c(8e5, 1.2e8),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(values = n_palette, labels = n_labels, drop = FALSE) +
  labs(
    title = "c",
    x = "Cumulative incidence",
    y = "Annual cost ($)",
    colour = "# of individuals"
  ) +
  ref_theme() +
  theme(
    legend.position = "none"
  )

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
    individuals_sampled = factor(individuals_sampled, levels = n_levels)
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
  geom_line(linewidth = 0.95, show.legend = FALSE) +
  geom_point(aes(shape = parameterization), size = 2.8) +
  scale_shape_manual(
    values = c("HIV estimate" = 16, "Simulated" = 15),
    guide = guide_legend(override.aes = list(colour = "black", size = 3.0))
  ) +
  scale_colour_manual(values = n_palette, labels = n_labels, drop = FALSE) +
  scale_x_log10(
    labels = label_log(base = 10),
    breaks = 10^(-6:-3),
    limits = c(8e-7, 1.2e-3),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(6:8),
    limits = c(8e5, 1.2e8),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "d",
    x = "RA of target pathogen",
    y = sprintf("Annual cost ($): CI = %s", percent_format(accuracy = 0.01)(b_ci_target)),
    shape = "Parameterization"
  ) +
  ref_theme() +
  theme(
    legend.position = c(0.76, 0.84),
  ) +
  guides(colour = "none")

# =============================================================================
# OPTIONAL SENSITIVITY PANEL (kept available, but not included in the final
# 3-panel layout because the reference figure only includes a1 / a2 / b)
# =============================================================================
fixed_n <- 50000

a3_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7) %>%
  filter(!(prob_infected_donates %in% c(0.01, 0.25, 0.3333333))) %>%
  filter(individuals_sampled == fixed_n) %>%
  filter(converged == TRUE) %>%
  mutate(prob_infected_donates = factor(
    prob_infected_donates,
    levels = sort(unique(prob_infected_donates))
  ))

prob_levels <- levels(a3_data$prob_infected_donates)
n_prob_levels <- length(prob_levels)
fixed_n_colour <- n_palette[as.character(fixed_n)]
prob_palette <- colorRampPalette(c("grey85", fixed_n_colour))(n_prob_levels)
names(prob_palette) <- prob_levels

p_a3_line2 <- ggplot(
  a3_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = prob_infected_donates,
    group = prob_infected_donates
  )
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_x_log10(
    breaks = ci_breaks,
    labels = ci_labels,
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(9:12),
    limits = c(7e8, 1.3e12),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = unname(prob_palette),
    breaks = c("0.1", "0.2", "0.5", "1"),
    labels = paste(c("1/10", "1/5", "1/2", "1")),
    guide = guide_legend(
      title = "RR of Infected\nPerson Donating",
      override.aes = list(size = 4, shape = 15)
    )
  ) +
  labs(
    title = "b",
    x = "Cumulative incidence",
    y = "Weekly sequencing depth",
    colour = "Prob. "
  ) +
  ref_theme() +
  theme(
    legend.position = c(0.78, 0.78)
  )

# =============================================================================
# FINAL 3-PANEL FIGURE MATCHING THE REFERENCE LAYOUT
# =============================================================================
final_plot <- cowplot::plot_grid(
  p_a1,
  p_a2,
  p_b,
  nrow = 1,
  align = "h",
  axis = "tb"
)

final_plot

# =============================================================================
# FINAL 4-PANEL FIGURE IN A 2 x 2 LAYOUT
# =============================================================================
final_plot <- cowplot::plot_grid(
  p_a1,
  p_a3_line2,
  p_a2,
  p_b,
  nrow = 2,
  align = "hv",
  axis = "tblr"
)

final_plot

# =============================================================================
# NEW PANEL: weekly sequencing depth vs prob_infected_donates at CI = 0.01%
# =============================================================================
new_ci_target <- 1e-4

new_data <- opt_seq_depth_summary %>%
  filter(target_cumulative_incidence == new_ci_target) %>%
  filter(abs(mu - hiv_mu) < 1e-7) %>%
  filter(converged == TRUE) %>%
  mutate(individuals_sampled = factor(individuals_sampled, levels = n_levels))

# Build breaks at the actual data values and compact diagonal-fraction labels
prob_breaks <- sort(unique(new_data$prob_infected_donates))
prob_labels <- do.call(expression, lapply(prob_breaks, function(v) {
  if (abs(v - 1) < 1e-9) {
    quote(1)
  } else {
    denom <- 1 / v
    denom <- if (abs(denom - round(denom)) < 1e-6) as.integer(round(denom)) else round(denom, 2)
    bquote(scriptstyle(""^1 * "/" * ""[.(denom)]))
  }
}))

p_new <- ggplot(
  new_data,
  aes(
    x = prob_infected_donates,
    y = optimal_depth,
    colour = individuals_sampled,
    group = individuals_sampled
  )
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = prob_breaks,
    labels = prob_labels,
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(9:12),
    limits = c(7e8, 1.3e12),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = n_palette,
    labels = n_labels,
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(linewidth = 1.0, size = 3.0)
    )
  ) +
  labs(
    title = "b",
    x = "Chance of infected individual donating relative to uninfected person",
    y = sprintf("Weekly sequencing depth: CI = %s", percent_format(accuracy = 0.01)(new_ci_target)),
    colour = "# of individuals"
  ) +
  ref_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 15)
  )

p_new

final_plot2 <- cowplot::plot_grid(
  p_a1,
  p_a2 + labs(title = "b"),
  p_new + labs(title = "c"),
  p_b,
  nrow = 2,
  align = "hv",
  axis = "tblr"
)

final_plot2

# Extract the row for 10,000 individuals, HIV mu, 0.01% CI, prob_donates = 1
row <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         individuals_sampled == 10000,
         target_cumulative_incidence == 1e-4,
         prob_infected_donates == 1,
         converged == TRUE)

depth_billion <- row$optimal_depth / 1e9
n_sampled     <- 10000

# Decompose into components
sample_acquisition <- cost_per_individual_sample * n_sampled * 52
sequencing         <- depth_billion * cost_per_billion_reads_sequencing * 52
sample_processing  <- cost_sample_processing_per_week * 52
labor              <- annual_labor_cost
data_processing    <- depth_billion * cost_per_billion_reads_data_processing * 52
storage            <- calc_cumulative_storage_cost(
  calc_storage_unit_cost(depth_billion, cost_per_billion_reads_storage_per_week),
  num_weeks_in_period = 52
)

# Build the table
cost_table <- tibble(
  Component = c("Sample acquisition", "Sequencing", "Labor",
                "Sample processing", "Data storage & processing"),
  Annual_Cost = c(sample_acquisition, sequencing, labor,
                  sample_processing, data_processing + storage)
) %>%
  mutate(
    Total = sum(Annual_Cost),
    Percentage = Annual_Cost / Total * 100
  )

# Add total row
cost_table <- bind_rows(
  cost_table %>% select(Component, Annual_Cost, Percentage),
  tibble(Component = "Total",
         Annual_Cost = sum(cost_table$Annual_Cost),
         Percentage = 100)
)

# Format for display
cost_table %>%
  mutate(
    Annual_Cost = paste0("$", comma(round(Annual_Cost))),
    Percentage  = paste0(round(Percentage, 1), "%")
  )
