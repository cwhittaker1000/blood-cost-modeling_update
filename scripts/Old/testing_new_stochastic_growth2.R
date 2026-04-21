library(ggplot2)
library(patchwork)
source("scripts/functions.R")

set.seed(123)
I0 <- 1
max_weeks <- 1000
n_surviving_target <- 100
infectious_weeks <- 52 * 5
shedding_weeks <- 12
doubling_time <- 52 # log(2) / 0.0155

# Derive R0 from doubling time
R0 <- calc_R0_from_doubling_time(doubling_time, infectious_weeks)
r <- calc_r_from_R0(R0, infectious_weeks)
cat(sprintf("R0 = %.4f, r = %.5f, doubling time = %.1f weeks, D = %d weeks\n",
            R0, r, doubling_time, infectious_weeks))

# =========================================================================
# Deterministic simple exponential
# =========================================================================
det_simple <- create_outbreak_tracker(I0, R0, max_weeks,
                                      growth_model = "deterministic_simple",
                                      infectious_weeks = infectious_weeks)
for (t in 1:max_weeks) det_simple$advance_to_week(t)
det_simple_df <- data.frame(
  week = 1:max_weeks,
  new_infections = sapply(1:max_weeks, function(t) det_simple$get_shedding(t, 1)),
  cumulative = sapply(1:max_weeks, det_simple$get_cumulative),
  model = "Deterministic simple"
)

# =========================================================================
# Deterministic renewal
# =========================================================================
det_renewal <- create_outbreak_tracker(I0, R0, max_weeks,
                                       growth_model = "deterministic_renewal",
                                       infectious_weeks = infectious_weeks)
for (t in 1:max_weeks) det_renewal$advance_to_week(t)
det_renewal_df <- data.frame(
  week = 1:max_weeks,
  new_infections = sapply(1:max_weeks, function(t) det_renewal$get_shedding(t, 1)),
  cumulative = sapply(1:max_weeks, det_renewal$get_cumulative),
  model = "Deterministic renewal"
)

# =========================================================================
# Stochastic renewal (collect surviving realisations)
# =========================================================================
surviving_list <- list()
total_attempts <- 0
while (length(surviving_list) < n_surviving_target) {
  total_attempts <- total_attempts + 1
  tr <- create_outbreak_tracker(I0, R0, max_weeks,
                                growth_model = "stochastic_renewal",
                                infectious_weeks = infectious_weeks)
  for (t in 1:max_weeks) tr$advance_to_week(t)
  if (tr$get_shedding(max_weeks, 1) > 0) {
    surviving_list[[length(surviving_list) + 1]] <- data.frame(
      week = 1:max_weeks,
      new_infections = sapply(1:max_weeks, function(t) tr$get_shedding(t, 1)),
      cumulative = sapply(1:max_weeks, tr$get_cumulative),
      realisation = length(surviving_list)
    )
  }
}
extinction_rate <- 1 - n_surviving_target / total_attempts
cat(sprintf(
  "Stochastic: %d surviving from %d attempts (extinction rate: %.1f%%)\n",
  n_surviving_target, total_attempts, 100 * extinction_rate
))

stoch_df <- do.call(rbind, surviving_list)

# Stochastic summary stats
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

# =========================================================================
# Combine deterministic models for plotting
# =========================================================================
det_df <- rbind(det_simple_df, det_renewal_df)
det_df$model <- factor(det_df$model, levels = c("Deterministic simple", "Deterministic renewal"))

# =========================================================================
# Row 1: New infections — all three models
# =========================================================================
p1 <- ggplot() +
  geom_ribbon(data = stoch_summary, aes(x = week, ymin = q05_new, ymax = q95_new),
              fill = "coral", alpha = 0.2) +
  geom_line(data = stoch_summary, aes(x = week, y = mean_new, colour = "Stochastic renewal (mean)"),
            linewidth = 0.8) +
  geom_line(data = det_df, aes(x = week, y = new_infections, colour = model),
            linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Deterministic simple" = "grey50",
                                 "Deterministic renewal" = "steelblue",
                                 "Stochastic renewal (mean)" = "coral")) +
  labs(x = "Week", y = "New infections + 1 (log scale)", colour = NULL,
       title = "Weekly new infections: three models compared",
       subtitle = sprintf("I₀ = %d, R₀ = %.3f, D = %d weeks, T₂ = %.1f weeks",
                          I0, R0, infectious_weeks, doubling_time)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# =========================================================================
# Row 1: Cumulative infections — all three models
# =========================================================================
p2 <- ggplot() +
  geom_ribbon(data = stoch_summary, aes(x = week, ymin = q05_cum, ymax = q95_cum),
              fill = "coral", alpha = 0.2) +
  geom_line(data = stoch_summary, aes(x = week, y = mean_cum, colour = "Stochastic renewal (mean)"),
            linewidth = 0.8) +
  geom_line(data = det_df, aes(x = week, y = cumulative, colour = model),
            linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Deterministic simple" = "grey50",
                                 "Deterministic renewal" = "steelblue",
                                 "Stochastic renewal (mean)" = "coral")) +
  labs(x = "Week", y = "Cumulative infections (log scale)", colour = NULL,
       title = "Cumulative infections: three models compared") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# =========================================================================
# Row 2: Individual stochastic trajectories with both deterministic lines
# =========================================================================
p3 <- ggplot() +
  geom_line(data = stoch_df,
            aes(x = week, y = new_infections, group = realisation),
            alpha = 0.05, colour = "coral") +
  geom_line(data = det_df,
            aes(x = week, y = new_infections, colour = model),
            linewidth = 1) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Deterministic simple" = "grey50",
                                 "Deterministic renewal" = "steelblue")) +
  labs(x = "Week", y = "New infections + 1 (log scale)", colour = NULL,
       title = sprintf("New infections: %d surviving stochastic trajectories", n_surviving_target)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

p4 <- ggplot() +
  geom_line(data = stoch_df,
            aes(x = week, y = cumulative, group = realisation),
            alpha = 0.05, colour = "coral") +
  geom_line(data = det_df,
            aes(x = week, y = cumulative, colour = model),
            linewidth = 1) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Deterministic simple" = "grey50",
                                 "Deterministic renewal" = "steelblue")) +
  labs(x = "Week", y = "Cumulative infections + 1 (log scale)", colour = NULL,
       title = "Cumulative infections: surviving stochastic trajectories") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# =========================================================================
# Row 3: Ratio of stochastic mean to each deterministic model
# =========================================================================
ratio_df <- data.frame(
  week = stoch_summary$week,
  ratio_vs_simple = stoch_summary$mean_new / (det_simple_df$new_infections + 1),
  ratio_vs_renewal = stoch_summary$mean_new / (det_renewal_df$new_infections + 1),
  ratio_cum_vs_simple = stoch_summary$mean_cum / det_simple_df$cumulative,
  ratio_cum_vs_renewal = stoch_summary$mean_cum / det_renewal_df$cumulative
)

p5 <- ggplot(ratio_df, aes(x = week)) +
  geom_line(aes(y = ratio_vs_simple, colour = "vs Deterministic simple"), linewidth = 0.6) +
  geom_line(aes(y = ratio_vs_renewal, colour = "vs Deterministic renewal"), linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(x = "Week", y = "Stochastic mean / Deterministic",
       colour = NULL,
       title = "Ratio check: new infections (stochastic mean vs each deterministic)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

p6 <- ggplot(ratio_df, aes(x = week)) +
  geom_line(aes(y = ratio_cum_vs_simple, colour = "vs Deterministic simple"), linewidth = 0.6) +
  geom_line(aes(y = ratio_cum_vs_renewal, colour = "vs Deterministic renewal"), linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(x = "Week", y = "Stochastic mean / Deterministic",
       colour = NULL,
       title = "Ratio check: cumulative infections") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

(p1 + p2)
(p3 + p4)

# =========================================================================
# Rolling shedding population (infections in last shedding_weeks weeks)
# =========================================================================

# Deterministic simple
det_simple_shedding <- data.frame(
  week = 1:max_weeks,
  shedding = sapply(1:max_weeks, function(t) det_simple$get_shedding(t, shedding_weeks)),
  model = "Deterministic simple"
)

# Deterministic renewal
det_renewal_shedding <- data.frame(
  week = 1:max_weeks,
  shedding = sapply(1:max_weeks, function(t) det_renewal$get_shedding(t, shedding_weeks)),
  model = "Deterministic renewal"
)

det_shedding_df <- rbind(det_simple_shedding, det_renewal_shedding)
det_shedding_df$model <- factor(det_shedding_df$model,
                                levels = c("Deterministic simple", "Deterministic renewal"))

# Stochastic: compute shedding per realisation
stoch_shedding_df <- do.call(rbind, lapply(surviving_list, function(d) {
  tr <- create_outbreak_tracker(I0, R0, max_weeks,
                                growth_model = "stochastic_renewal",
                                infectious_weeks = infectious_weeks)
  # We already have the data in d, so compute rolling sum directly
  new_inf <- d$new_infections
  shedding <- sapply(1:max_weeks, function(t) {
    start <- max(1, t - shedding_weeks + 1)
    sum(new_inf[start:t])
  })
  data.frame(week = 1:max_weeks, shedding = shedding, realisation = d$realisation[1])
}))

# Stochastic summary
stoch_shedding_summary <- do.call(rbind, lapply(split(stoch_shedding_df, stoch_shedding_df$week), function(d) {
  data.frame(
    week = d$week[1],
    mean_shedding = mean(d$shedding),
    median_shedding = quantile(d$shedding, 0.5),
    q05_shedding = quantile(d$shedding, 0.05),
    q95_shedding = quantile(d$shedding, 0.95)
  )
}))

# Plot: shedding population over time
p_shedding <- ggplot() +
  geom_ribbon(data = stoch_shedding_summary,
              aes(x = week, ymin = q05_shedding, ymax = q95_shedding),
              fill = "coral", alpha = 0.2) +
  geom_line(data = stoch_shedding_summary,
            aes(x = week, y = median_shedding, colour = "Stochastic renewal (mean)"),
            linewidth = 0.8) +
  geom_line(data = det_shedding_df,
            aes(x = week, y = shedding, colour = model),
            linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Deterministic simple" = "grey50",
                                 "Deterministic renewal" = "steelblue",
                                 "Stochastic renewal (mean)" = "coral")) +
  labs(x = "Week",
       y = "Shedding population (log scale)",
       colour = NULL,
       title = sprintf("Rolling shedding population (last %d weeks of infections)", shedding_weeks),
       subtitle = sprintf("I₀ = %d, R₀ = %.3f, D = %d weeks, shedding = %d weeks",
                          I0, R0, infectious_weeks, shedding_weeks)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

p_shedding
