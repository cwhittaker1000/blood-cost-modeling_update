library(tidyverse)
library(furrr)
library(scales)
library(dplyr)

# =============================================================================
# Load data
# =============================================================================

opt_seq_depth_summary <- read_tsv("../results/summarized_results.tsv")

# =============================================================================
# PLOTS
# =============================================================================

# Extract data for a single mew closest to our estimated mew, and remove any
# parameter combinations that did not converge (hit our ceiling)
select_mew_read_threshold_data <- opt_seq_depth_summary %>%
  filter(mew == 1e-6, read_threshold == 100) %>%
  mutate(
    # Flag values that hit the ceiling
    hit_ceiling = optimal_depth >= 1e12
  ) %>%
  filter(!hit_ceiling)
depth_plot <- ggplot(
  select_mew_read_threshold_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1, linetype = "solid") +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(NA, 1e12)
  ) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  labs(
    title = "",
    y = "Sequencing depth per week to detect HIV-like pathogen",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
depth_plot

ggsave("../output/depth_plot.png", depth_plot, width = 8, height = 6)


## Cost analysis for a specific depth
seq_c <- 2500 / 1e9
proc_c <- 170
sample_c <- 15

# Add the cost estimate column
select_mew_read_threshold_data_with_cost <- select_mew_read_threshold_data %>%
  filter(!hit_ceiling) %>% # Exclude points that didn't converge
  mutate(
    total_cost_per_year = calculate_total_cost_per_year(
      seq_depth = optimal_depth,
      cost_of_seq = seq_c,
      cost_of_proc = proc_c,
      cost_of_sample = sample_c,
      num_samples = batch_size
    )
  )
cost_plot <- ggplot(
  select_mew_read_threshold_data_with_cost,
  aes(
    x = target_cumulative_incidence,
    y = total_cost_per_year,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(
    breaks = c(10000000, 100000000),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  labs(
    title = "",
    y = "Total cost per year to detect HIV-like pathogen ($)",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
cost_plot

ggsave("../output/cost_plot.png", cost_plot, width = 8, height = 6)

# Only filter by read threshold so that we can show a plot
# comparing different mews
select_read_threshold_data <- opt_seq_depth_summary %>%
  filter(read_threshold == 100) %>%
  mutate(
    # Flag values that hit the ceiling
    hit_ceiling = optimal_depth >= 1e12
  ) %>%
  filter(!hit_ceiling)

facet_depth_plot <- ggplot(
  select_read_threshold_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_point(
    size = 3,
  ) +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(NA, 1e12)
  ) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  labs(
    title = "",
    y = "Sequencing depth per week to detect HIV-like pathogen",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 15) +
  facet_grid(. ~ mew) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
facet_depth_plot

ggsave(
  "../output/facet_mew_depth_plot.png",
  facet_depth_plot,
  width = 25,
  height = 8
)


## Cost analysis for multiple sample costs
sample_costs <- c(3, 12, 50, 100)

# Add the cost estimate column
select_mew_read_threshold_data_with_varied_costs <- select_read_threshold_data %>%
  filter(!hit_ceiling) %>% # Exclude points that didn't converge
  crossing(sample_cost = sample_costs) %>% # Create all combinations
  mutate(
    total_cost_per_year = find_total_cost(
      seq_depth = optimal_depth,
      cost_of_seq = seq_c,
      cost_of_proc = proc_c,
      cost_of_sample = sample_cost,
      num_samples = batch_size
    ) *
      52
  )

# Compare costs for different mews with different sample costs
facet_cost_plot <- ggplot(
  select_mew_read_threshold_data_with_varied_costs,
  aes(
    x = target_cumulative_incidence,
    y = total_cost_per_year,
    colour = factor(batch_size),
    group = batch_size
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_brewer(
    palette = "Accent",
    name = "Number of individuals \nin weekly sequencing\nbatch"
  ) +
  #  facet_wrap(~ paste0("Sample cost: $", sample_cost), scales = "free_y") +
  facet_grid(sample_cost ~ mew) +
  labs(
    title = "Annual cost with varying sample collection costs",
    y = "Total cost per year to detect HIV-like pathogen ($)",
    x = "Cumulative incidence (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    legend.position = "right"
  )
facet_cost_plot

ggsave(
  "../output/facet_cost_plot.png",
  facet_cost_plot,
  width = 25,
  height = 15
)
