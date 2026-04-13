library(ggplot2)
library(tidyr)
source("scripts/functions.R")

# =============================================================================
# QUICK COMPARISON: simple vs renewal deterministic across parameter sweep
# Anchor: doubling time = 46 weeks (r ≈ 0.0155/week), held constant.
# For each infectious_weeks D, derive the R0 that yields that r via Euler-Lotka.
# This way the simple model produces an identical trajectory across all D
# (since it only depends on r), and we can isolate how D changes the renewal
# model's behaviour.
# =============================================================================

target_doubling_time <- 46  # weeks
target_r             <- log(2) / target_doubling_time

# Parameter grid
init_inf_vals  <- c(1, 10, 100)
inf_weeks_vals <- c(12, 52, 104, 260)
max_weeks_sim  <- 1000
models         <- c("deterministic_simple", "deterministic_renewal")

# For each D, compute the R0 that produces the target r
R0_for_D <- sapply(inf_weeks_vals, function(D) {
  calc_R0_from_doubling_time(target_doubling_time, D)
})
names(R0_for_D) <- inf_weeks_vals

cat("Anchor: doubling time =", target_doubling_time, "weeks (r =",
    round(target_r, 5), "/week)\n")
cat("Implied R0 for each D:\n")
print(round(R0_for_D, 4))

sweep_grid <- expand.grid(
  initial_infections = init_inf_vals,
  infectious_weeks   = inf_weeks_vals,
  growth_model       = models,
  stringsAsFactors   = FALSE
)
# Look up the corresponding R0 for each row's D
sweep_grid$R0 <- R0_for_D[as.character(sweep_grid$infectious_weeks)]

# Run all trackers and collect weekly new infections + cumulative
sweep_results <- do.call(rbind, lapply(seq_len(nrow(sweep_grid)), function(i) {
  row <- sweep_grid[i, ]
  tracker <- create_outbreak_tracker(
    initial_infections = row$initial_infections,
    R0                 = row$R0,
    max_weeks          = max_weeks_sim,
    growth_model       = row$growth_model,
    infectious_weeks   = row$infectious_weeks
  )
  for (t in 1:max_weeks_sim) tracker$advance_to_week(t)
  
  data.frame(
    week               = 1:max_weeks_sim,
    new_infections     = sapply(1:max_weeks_sim, function(t) tracker$get_shedding(t, 1)),
    cumulative         = sapply(1:max_weeks_sim, tracker$get_cumulative),
    initial_infections = row$initial_infections,
    infectious_weeks   = row$infectious_weeks,
    R0                 = row$R0,
    growth_model       = row$growth_model
  )
}))

# Pretty facet labels
sweep_results$inf_weeks_lab <- factor(
  sprintf("D = %d wk (R0 = %.3f)",
          sweep_results$infectious_weeks,
          sweep_results$R0),
  levels = sprintf("D = %d wk (R0 = %.3f)",
                   inf_weeks_vals,
                   R0_for_D)
)
sweep_results$init_lab <- factor(
  paste0("I0 = ", sweep_results$initial_infections),
  levels = paste0("I0 = ", init_inf_vals)
)

# Plot: new infections, faceted by I0 (rows) × D (cols), coloured by model
p_new <- ggplot(sweep_results,
                aes(x = week, y = new_infections,
                    colour = growth_model)) +
  geom_line(linewidth = 0.7) +
  facet_grid(init_lab ~ inf_weeks_lab, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c(
    "deterministic_simple"  = "grey40",
    "deterministic_renewal" = "steelblue"
  )) +
  labs(x = "Week", y = "New infections (log)",
       colour = "Model",
       title = "Simple vs renewal deterministic, anchored on r = 0.0155/week",
       subtitle = sprintf("Doubling time = %d weeks; R0 derived per D via Euler-Lotka",
                          target_doubling_time)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_new)

# Same plot but for cumulative infections
p_cum <- p_new %+% sweep_results +
  aes(y = cumulative) +
  labs(y = "Cumulative infections (log)")

print(p_cum)


x <- ggplot(sweep_results, aes(x = week, y = new_infections, colour = inf_weeks_lab)) +
  geom_line(linewidth = 0.7) +
  facet_grid(growth_model ~ init_lab, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  # scale_colour_manual(values = c(
  #   "deterministic_simple"  = "grey40",
  #   "deterministic_renewal" = "steelblue"
  # )) +
  labs(x = "Week", y = "New infections (log)",
       colour = "Model",
       title = "Simple vs renewal deterministic, anchored on r = 0.0155/week",
       subtitle = sprintf("Doubling time = %d weeks; R0 derived per D via Euler-Lotka",
                          target_doubling_time)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p_cum <- x %+% sweep_results +
  aes(y = cumulative) +
  labs(y = "Cumulative infections (log)")

x
print(p_cum)

# =============================================================================
# Rolling shedding population: infections in the last 12 weeks
# =============================================================================

shedding_window <- 12

# Recompute trackers and pull the rolling shedding sum directly from each one.
# This uses the tracker's get_shedding(t, window) method, which is the same
# function used downstream in the simulation pipeline.
shedding_results <- do.call(rbind, lapply(seq_len(nrow(sweep_grid)), function(i) {
  row <- sweep_grid[i, ]
  tracker <- create_outbreak_tracker(
    initial_infections = row$initial_infections,
    R0                 = row$R0,
    max_weeks          = max_weeks_sim,
    growth_model       = row$growth_model,
    infectious_weeks   = row$infectious_weeks
  )
  for (t in 1:max_weeks_sim) tracker$advance_to_week(t)
  
  data.frame(
    week               = 1:max_weeks_sim,
    shedding           = sapply(1:max_weeks_sim,
                                function(t) tracker$get_shedding(t, shedding_window)),
    initial_infections = row$initial_infections,
    infectious_weeks   = row$infectious_weeks,
    R0                 = row$R0,
    growth_model       = row$growth_model
  )
}))

# Reuse the same facet labels
shedding_results$inf_weeks_lab <- factor(
  sprintf("D = %d wk (R0 = %.3f)",
          shedding_results$infectious_weeks,
          shedding_results$R0),
  levels = sprintf("D = %d wk (R0 = %.3f)",
                   inf_weeks_vals,
                   R0_for_D)
)
shedding_results$init_lab <- factor(
  paste0("I0 = ", shedding_results$initial_infections),
  levels = paste0("I0 = ", init_inf_vals)
)

p_shedding <- ggplot(shedding_results,
                     aes(x = week, y = shedding,
                         colour = inf_weeks_lab)) +
  geom_line(linewidth = 0.7) +
  facet_grid(growth_model ~ init_lab, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Week",
       y = sprintf("Shedding population (last %d wk, log)", shedding_window),
       colour = "Infectious window",
       title = sprintf("Rolling shedding population (window = %d weeks)",
                       shedding_window),
       subtitle = sprintf("Doubling time = %d weeks; R0 derived per D via Euler-Lotka",
                          target_doubling_time)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_shedding)

