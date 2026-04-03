library(ggplot2)
library(patchwork)
source("scripts/functions.R")

set.seed(123)
r <- 0.0155
I0 <- 1
max_weeks <- 500
n_realisations <- 250
infectious_weeks = 12

# Run many stochastic realisations
stoch_list <- lapply(1:n_realisations, function(i) {
  tr <- create_outbreak_tracker(I0, r, max_weeks, stochastic_growth = TRUE, infectious_weeks)
  for (t in 1:max_weeks) tr$advance_to_week(t)
  data.frame(
    week = 1:max_weeks,
    new_infections = sapply(1:max_weeks, function(t) tr$get_shedding(t, 1)),
    cumulative = sapply(1:max_weeks, tr$get_cumulative),
    realisation = i
  )
})
stoch_df <- do.call(rbind, stoch_list)

# Deterministic expectation
det <- create_outbreak_tracker(I0, r, max_weeks, stochastic_growth = FALSE, infectious_weeks)
for (t in 1:max_weeks) det$advance_to_week(t)
det_df <- data.frame(
  week = 1:max_weeks,
  new_infections = sapply(1:max_weeks, function(t) det$get_shedding(t, 1)),
  cumulative = sapply(1:max_weeks, det$get_cumulative)
)

# Calculate stochastic mean and percentile bands per week
stoch_summary <- do.call(rbind, lapply(split(stoch_df, stoch_df$week), function(d) {
  data.frame(
    week = d$week[1],
    mean_new = mean(d$new_infections),
    median_new = median(d$new_infections),
    q05_new = quantile(d$new_infections, 0.05),
    q95_new = quantile(d$new_infections, 0.95),
    mean_cum = mean(d$cumulative),
    median_cum = median(d$cumulative),
    q05_cum = quantile(d$cumulative, 0.05),
    q95_cum = quantile(d$cumulative, 0.95)
  )
}))

# Row 1: Summary plots
p1 <- ggplot(stoch_summary, aes(x = week)) +
  geom_ribbon(aes(ymin = q05_new, ymax = q95_new), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = mean_new, colour = "Stochastic mean"), linewidth = 0.8) +
  geom_line(data = det_df, aes(y = new_infections, colour = "Deterministic"), linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Deterministic" = "steelblue", "Stochastic mean" = "coral")) +
  labs(x = "Week", y = "New infections (log scale)", colour = NULL,
       title = sprintf("New infections: summary of %d realisations", n_realisations),
       subtitle = sprintf("I₀ = %d, r = %.4f/week. Ribbon = 5th–95th percentile", I0, r)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

p2 <- ggplot(stoch_summary, aes(x = week)) +
  geom_ribbon(aes(ymin = q05_cum, ymax = q95_cum), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = mean_cum, colour = "Stochastic mean"), linewidth = 0.8) +
  geom_line(data = det_df, aes(y = cumulative, colour = "Deterministic"), linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Deterministic" = "steelblue", "Stochastic mean" = "coral")) +
  labs(x = "Week", y = "Cumulative infections (log scale)", colour = NULL,
       title = "Cumulative infections: summary",
       subtitle = sprintf("I₀ = %d, r = %.4f/week. Ribbon = 5th–95th percentile", I0, r)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# Row 2: Individual trajectories
p3 <- ggplot(stoch_df, aes(x = week, y = new_infections + 1, group = realisation)) +
  geom_line(alpha = 0.05, colour = "coral") +
  geom_line(data = det_df, aes(x = week, y = new_infections + 1, group = NULL),
            colour = "steelblue", linewidth = 1) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Week", y = "New infections + 1 (log scale)",
       title = sprintf("New infections: %d individual trajectories", n_realisations)) +
  theme_minimal(base_size = 13)

p4 <- ggplot(stoch_df, aes(x = week, y = cumulative, group = realisation)) +
  geom_line(alpha = 0.05, colour = "coral") +
  geom_line(data = det_df, aes(x = week, y = cumulative, group = NULL),
            colour = "steelblue", linewidth = 1) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Week", y = "Cumulative infections (log scale)",
       title = "Cumulative infections: individual trajectories") +
  theme_minimal(base_size = 13)

# Row 3: Ratio check
p5 <- ggplot(ratio_df, aes(x = week)) +
  geom_line(aes(y = ratio_new, colour = "New infections"), linewidth = 0.6) +
  geom_line(aes(y = ratio_cum, colour = "Cumulative"), linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(x = "Week", y = "Stochastic mean / Deterministic",
       colour = NULL, title = "Ratio check (should be ~1)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ratio_df <- data.frame(
  week = stoch_summary$week,
  ratio_new = stoch_summary$mean_new / det_df$new_infections,
  ratio_cum = stoch_summary$mean_cum / det_df$cumulative
)

p5 <- ggplot(ratio_df, aes(x = week)) +
  geom_line(aes(y = ratio_new, colour = "New infections"), linewidth = 0.6) +
  geom_line(aes(y = ratio_cum, colour = "Cumulative"), linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(x = "Week", y = "Stochastic mean / Deterministic",
       colour = NULL, title = "Ratio check (should be ~1)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

(p1 + p2) / (p3 + p4) / p5

# Calculate extinction rate
extinct <- sapply(stoch_list, function(d) {
  d$new_infections[d$week == max_weeks] == 0
})
cat(sprintf(
  "\nExtinction rate: %d/%d (%.1f%%)\n",
  sum(extinct), n_realisations, 100 * mean(extinct)
))
