# Loading required libraries
library(readr); library(dplyr); library(ggplot2); library(scales); library(cowplot)
source("scripts/functions.R")

# =============================================================================
# Load Boston results
# =============================================================================
opt_seq_depth_summary <- read_tsv(
  "results/boston_summarized_simulation_results.tsv",
  show_col_types = FALSE
)

# The HIV-derived mu value (arithmetic mean relative abundance from Piantadosi data).
# Used to filter rows that correspond to the HIV estimate as opposed to the
# synthetic sweep values. Tolerance-based match since mu is a float.
hiv_mu <- 1.26e-5

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

# =============================================================================
# ADD ANNUAL COST TO THE SUMMARY TABLE
# =============================================================================
# optimal_depth is in reads per week. The cost function takes billions of reads
# per week, so divide by 1e9.
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
# PANEL A1: weekly sequencing depth vs cumulative incidence
# Filter: HIV mu, prob_infected_donates == 1
# =============================================================================
a1_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         prob_infected_donates == 1,
         converged == TRUE) %>%
  mutate(individuals_sampled = factor(
    individuals_sampled,
    levels = sort(unique(individuals_sampled))
  )) %>%
  filter(individuals_sampled != 500)

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
    breaks = 10^(7:13)
  ) +
  labs(
    x = "Cumulative incidence",
    y = "Weekly sequencing depth",
    colour = "# of individuals",
    title = "a1"
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
# PANEL A2: annual cost vs cumulative incidence
# Filter: HIV mu, prob_infected_donates == 1
# =============================================================================
a2_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         prob_infected_donates == 1,
         converged == TRUE) %>%
  mutate(individuals_sampled = factor(
    individuals_sampled,
    levels = sort(unique(individuals_sampled))
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
    breaks = 10^(5:9)
  ) +
  labs(
    x = "Cumulative incidence",
    y = "Annual cost ($)",
    colour = "# of individuals",
    title = "a2"
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
# Filter: prob_infected_donates == 1
# =============================================================================
# Points are shaped by "Parameterization": the HIV estimate (1.26e-5) vs the
# synthetic sweep values. In the current grid the HIV mu sits alongside the
# synthetic points on the same curve.
b_ci_target <- 5e-4

b_data <- opt_seq_depth_summary %>%
  filter(target_cumulative_incidence == b_ci_target,
         prob_infected_donates == 1,
         converged == TRUE) %>%
  mutate(
    parameterization = ifelse(
      abs(mu - hiv_mu) < 1e-7, "HIV estimate", "Simulated"
    ),
    individuals_sampled = factor(
      individuals_sampled,
      levels = sort(unique(individuals_sampled))
    )
  ) %>%
  filter(individuals_sampled != 500)

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
  scale_x_log10(
    labels = label_log(base = 10),
    breaks = 10^(-8:-3)
  ) +
  scale_y_log10(
    labels = label_log(base = 10),
    breaks = 10^(5:9)
  ) +
  labs(
    x = "RA of target pathogen",
    y = sprintf("Annual cost ($), CI = %s%%", signif(b_ci_target * 100, 2)),
    colour = "# of individuals",
    shape = "Parameterization",
    title = "b"
  ) # +
  # theme_minimal(base_size = 13) +
  # theme(
  #   legend.position = c(0.78, 0.78),
  #   legend.background = element_rect(fill = "white", colour = "grey80"),
  #   legend.box = "vertical",
  #   plot.title = element_text(face = "bold", hjust = 0),
  #   axis.text.x = element_text(angle = 45, hjust = 1)
  # ) +
  # guides(colour = "none")  # colour legend already shown in a1

p_b

# =============================================================================
# PANEL C: weekly sequencing depth vs cumulative incidence,
# fixed to 4,000 individuals/wk, coloured by P(infected donates).
# Filter: HIV mu only.
# =============================================================================
# Pull the exact colour ggplot assigns to "4000" in panel a1 so the gradient
# in panel c is anchored to it. Discrete palette (scale_colour_hue) assigns
# evenly-spaced hues; we just grab the right index.
a1_levels   <- levels(a1_data$individuals_sampled)
a1_palette  <- scales::hue_pal()(length(a1_levels))
colour_4000 <- a1_palette[which(a1_levels == "5000")]

c_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         individuals_sampled == 5000,
         converged == TRUE)

p_c <- ggplot(
  c_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = prob_infected_donates,
    group = prob_infected_donates
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
  scale_colour_gradient(
    low    = scales::alpha(colour_4000, 0.25),
    high   = colour_4000,
    limits = c(0, 1),
    breaks = seq(0.25, 1, 0.25)
  ) +
  labs(
    x = "Cumulative incidence",
    y = "Weekly sequencing depth",
    colour = "P(infected donates)",
    title = "c: 4,000 sampled/wk"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.78, 0.78),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_c

c_data <- opt_seq_depth_summary %>%
  filter(abs(mu - hiv_mu) < 1e-7,
         individuals_sampled == 5000,
         converged == TRUE) %>%
  mutate(prob_infected_donates = factor(prob_infected_donates,
                                        levels = sort(unique(prob_infected_donates))))

# Build the discrete palette by interpolating between grey and colour_4000
n_levels <- length(levels(c_data$prob_infected_donates))
discrete_palette <- colorRampPalette(c("grey85", colour_4000))(n_levels)
names(discrete_palette) <- levels(c_data$prob_infected_donates)


p_c <- ggplot(
  c_data,
  aes(
    x = target_cumulative_incidence,
    y = optimal_depth,
    colour = prob_infected_donates,
    group = prob_infected_donates
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
  scale_colour_manual(
    values = discrete_palette,
    guide  = guide_legend(
      override.aes = list(size = 4, shape = 15)  # filled squares
    )
  ) +
  labs(
    x = "Cumulative incidence",
    y = "Weekly sequencing depth",
    colour = "P(infected donates)",
    title = "b"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.78, 0.78),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_c

# =============================================================================
# Compose
# =============================================================================
cowplot::plot_grid(p_a1, p_c, p_a2, p_b, nrow = 1)
