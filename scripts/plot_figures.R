# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)

source("scripts/functions.R")

hiv_mu <- 1.26e-5

# =============================================================================
# EPIDEMIC TRAJECTORY CHECK
# =============================================================================
max_weeks_sim <- 1000
I0 <- 100
weekly_growth_rate <- 0.0155
doubling_time <- log(2) / weekly_growth_rate
pop <- 3.4e8
D_vals <- c(12, 52)

traj_df <- do.call(rbind, lapply(D_vals, function(D) {
  R0 <- calc_R0_from_doubling_time(doubling_time, D)
  tr <- create_outbreak_tracker(I0, R0, max_weeks_sim,
                                growth_model = "deterministic_renewal",
                                infectious_weeks = D)
  for (t in 1:max_weeks_sim) tr$advance_to_week(t)
  data.frame(
    week       = 1:max_weeks_sim,
    cumulative = sapply(1:max_weeks_sim, tr$get_cumulative),
    D          = sprintf("D = %d weeks", D)
  )
}))

ci_targets <- c(1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3)

p_trajectory <- ggplot(traj_df, aes(x = week, y = cumulative / pop, colour = D)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = ci_targets, linetype = "dashed", colour = "grey50") +
  scale_y_log10(labels = percent_format(accuracy = 0.0001)) +
  labs(x = "Week", y = "Cumulative incidence",
       title = "Renewal model: time to reach CI targets",
       subtitle = sprintf("I0 = %d, r = %.4f/wk, max_weeks = %d",
                          I0, weekly_growth_rate, max_weeks_sim)) +
  theme_minimal(base_size = 13)

p_trajectory

# =============================================================================
# LOAD RESULTS
# =============================================================================
opt_seq_depth_summary <- read_tsv(
  "results/summarized_simulation_results_renewal_probdonates.tsv",
  show_col_types = FALSE
)

# =============================================================================
# COLOUR PALETTE
# =============================================================================
image_matched_palette <- c(
  "#4077AB", "#EF6677", "#1F8942", "#CCBB44"
)

n_levels <- as.character(sort(unique(opt_seq_depth_summary$individuals_sampled)))

if (length(n_levels) != 4) {
  stop(sprintf(
    "Expected 4 individuals_sampled levels, but found %d: %s",
    length(n_levels), paste(n_levels, collapse = ", ")
  ))
}

n_palette <- image_matched_palette
names(n_palette) <- n_levels
n_labels <- comma(as.numeric(n_levels))

# =============================================================================
# SHARED FORMATTING
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
# COST PARAMETERS & ANNUAL COST
# =============================================================================
cost_per_individual_sample              <- 5.00
cost_sample_processing_per_week         <- 500
cost_per_billion_reads_sequencing       <- 2500
cost_per_billion_reads_data_processing  <- 1
cost_per_billion_reads_storage_per_week <- 0.25
annual_labor_cost                       <- 650000

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
# PANEL A: weekly sequencing depth vs cumulative incidence
# =============================================================================
a1_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         prob_infected_donates == 1,
         converged == TRUE) %>%
  mutate(individuals_sampled = factor(individuals_sampled, levels = n_levels))

p_a <- ggplot(
  a1_data,
  aes(x = target_cumulative_incidence, y = optimal_depth,
      colour = individuals_sampled, group = individuals_sampled)
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_x_log10(breaks = ci_breaks, labels = ci_labels,
                expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_log10(labels = label_log(base = 10), breaks = 10^(9:12),
                limits = c(7e8, 1.3e12), expand = expansion(mult = c(0, 0.02))) +
  scale_colour_manual(values = n_palette, labels = n_labels, drop = FALSE,
                      guide = guide_legend(override.aes = list(linewidth = 1.0, size = 3.0))) +
  labs(title = "a", x = "Cumulative incidence",
       y = "Weekly sequencing depth", colour = "# of individuals") +
  ref_theme() +
  theme(legend.position = c(0.73, 0.80))

# =============================================================================
# PANEL B: annual cost vs cumulative incidence
# =============================================================================
a2_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         prob_infected_donates == 1,
         converged == TRUE) %>%
  mutate(individuals_sampled = factor(individuals_sampled, levels = n_levels))

p_b <- ggplot(
  a2_data,
  aes(x = target_cumulative_incidence, y = annual_cost,
      colour = individuals_sampled, group = individuals_sampled)
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_x_log10(breaks = ci_breaks, labels = ci_labels,
                expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_log10(labels = label_log(base = 10), breaks = 10^(6:8),
                limits = c(8e5, 1.2e8), expand = expansion(mult = c(0, 0.02))) +
  scale_colour_manual(values = n_palette, labels = n_labels, drop = FALSE) +
  labs(title = "b", x = "Cumulative incidence",
       y = "Annual cost ($)", colour = "# of individuals") +
  ref_theme() +
  theme(legend.position = "none")

# =============================================================================
# PANEL C: weekly sequencing depth vs prob_infected_donates at CI = 0.01%
# =============================================================================
new_ci_target <- 1e-4

new_data <- opt_seq_depth_summary %>%
  filter(target_cumulative_incidence == new_ci_target,
         abs(mu - hiv_mu) < 1e-7,
         converged == TRUE) %>%
  mutate(individuals_sampled = factor(individuals_sampled, levels = n_levels))

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

p_c <- ggplot(
  new_data,
  aes(x = prob_infected_donates, y = optimal_depth,
      colour = individuals_sampled, group = individuals_sampled)
) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = prob_breaks, labels = prob_labels,
                     expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_log10(labels = label_log(base = 10), breaks = 10^(9:12),
                limits = c(7e8, 1.3e12), expand = expansion(mult = c(0, 0.02))) +
  scale_colour_manual(values = n_palette, labels = n_labels, drop = FALSE,
                      guide = guide_legend(override.aes = list(linewidth = 1.0, size = 3.0))) +
  labs(title = "c",
       x = "Chance of infected individual donating relative to uninfected person",
       y = sprintf("Weekly sequencing depth: CI = %s", percent_format(accuracy = 0.01)(new_ci_target)),
       colour = "# of individuals") +
  ref_theme() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 15))

# =============================================================================
# PANEL D: annual cost vs mu at CI = 0.01%
# =============================================================================
b_ci_target <- 1e-4

b_data <- opt_seq_depth_summary %>%
  filter(target_cumulative_incidence == b_ci_target,
         prob_infected_donates == 1,
         converged == TRUE) %>%
  mutate(
    parameterization = ifelse(abs(mu - hiv_mu) < 1e-7, "HIV estimate", "Simulated"),
    individuals_sampled = factor(individuals_sampled, levels = n_levels)
  ) %>%
  filter(individuals_sampled != 2e+03)

p_d <- ggplot(
  b_data,
  aes(x = mu, y = annual_cost,
      colour = individuals_sampled, group = individuals_sampled)
) +
  geom_line(linewidth = 0.95, show.legend = FALSE) +
  geom_point(aes(shape = parameterization), size = 2.8) +
  scale_shape_manual(values = c("HIV estimate" = 16, "Simulated" = 15),
                     guide = guide_legend(override.aes = list(colour = "black", size = 3.0))) +
  scale_colour_manual(values = n_palette, labels = n_labels, drop = FALSE) +
  scale_x_log10(labels = label_log(base = 10), breaks = 10^(-6:-3),
                limits = c(8e-7, 1.2e-3), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_log10(labels = label_log(base = 10), breaks = 10^(6:8),
                limits = c(8e5, 1.2e8), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "d", x = "RA of target pathogen",
       y = sprintf("Annual cost ($): CI = %s", percent_format(accuracy = 0.01)(b_ci_target)),
       shape = "Parameterization") +
  ref_theme() +
  theme(legend.position = c(0.76, 0.84)) +
  guides(colour = "none")

# =============================================================================
# FINAL 4-PANEL FIGURE
# =============================================================================
final_plot <- cowplot::plot_grid(
  p_a, p_b, p_c, p_d,
  nrow = 2, align = "hv", axis = "tblr"
)

final_plot

# =============================================================================
# COST TABLE: 10,000 individuals, HIV mu, 0.01% CI
# =============================================================================
row <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         individuals_sampled == 10000,
         target_cumulative_incidence == 1e-4,
         prob_infected_donates == 1,
         converged == TRUE)

depth_billion <- row$optimal_depth / 1e9
n_sampled     <- 10000

sample_acquisition <- cost_per_individual_sample * n_sampled * 52
sequencing         <- depth_billion * cost_per_billion_reads_sequencing * 52
sample_processing  <- cost_sample_processing_per_week * 52
labor              <- annual_labor_cost
data_processing    <- depth_billion * cost_per_billion_reads_data_processing * 52
storage            <- calc_cumulative_storage_cost(
  calc_storage_unit_cost(depth_billion, cost_per_billion_reads_storage_per_week),
  num_weeks_in_period = 52
)

cost_table <- tibble(
  Component = c("Sample acquisition", "Sequencing", "Labor",
                "Sample processing", "Data storage & processing"),
  Annual_Cost = c(sample_acquisition, sequencing, labor,
                  sample_processing, data_processing + storage)
) %>%
  mutate(Total = sum(Annual_Cost),
         Percentage = Annual_Cost / Total * 100)

cost_table <- bind_rows(
  cost_table %>% select(Component, Annual_Cost, Percentage),
  tibble(Component = "Total",
         Annual_Cost = sum(cost_table$Annual_Cost),
         Percentage = 100)
)

cost_table %>%
  mutate(
    Annual_Cost = paste0("$", comma(round(Annual_Cost))),
    Percentage  = paste0(round(Percentage, 1), "%")
  )
