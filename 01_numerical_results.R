################################################################################
# 01_numerical_results.R
#
# Reproduces: Main text Section 3.2 — Numerical results (Figure 5)
#
# Computes the weighted average hazard difference A(tau) as a function of
# follow-up percentile for three scenarios (all with HR_u = 0.5):
#   Solid line:  No cure in either group (pi0 = 0, pi1 = 0)
#   Dashed line: Cure in treatment only (pi0 = 0, pi1 = 0.2)
#   Dotted line: Cure in both groups (pi0 = 0.2, OR = 1.5)
#
# Control uncured distribution: Weibull(shape = 2, scale = 1)
#
# Outputs:
#   figure_weighted_diff_threecurves_marked.pdf  (Figure 5)
#   figure_weighted_diff_threecurves_marked.png
#
# Required packages: ggplot2
################################################################################

library(ggplot2)

################################################################################
# Helper functions
################################################################################

# Weibull survival function: S(t) = exp(-(t/scale)^shape)
S_weibull <- function(t, shape, scale) {
  exp(-(t / scale)^shape)
}

# Weibull density function
f_weibull <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1) * exp(-(t / scale)^shape)
}

# Mixture cure model survival: S(t) = pi + (1 - pi) * S_u(t)
S_cure <- function(t, pi_g, shape, scale) {
  pi_g + (1 - pi_g) * S_weibull(t, shape, scale)
}

# Mixture cure model hazard: h(t) = (1 - pi) * f_u(t) / S(t)
h_cure <- function(t, pi_g, shape, scale) {
  S_uc <- S_weibull(t, shape, scale)
  f_uc <- f_weibull(t, shape, scale)
  S_mix <- pi_g + (1 - pi_g) * S_uc
  (1 - pi_g) * f_uc / S_mix
}

# Compute A(tau) = integral_0^tau w(t) * {h_1(t) - h_0(t)} dt
# via left Riemann sum, where
#   w(t) = G(t) * 2 * S_0(t) * S_1(t) / (S_0(t) + S_1(t))
#   G(t) = max(0, 1 - t / tau_max)   [censoring survival]
compute_A_tau <- function(tau_eval, tau_max,
                          pi0, shape0, scale0,
                          pi1, shape1, scale1,
                          n_points = 500) {
  if (tau_eval <= 0) return(0)
  t_grid <- seq(0, tau_eval, length.out = n_points + 1)
  dt <- t_grid[2] - t_grid[1]

  A_val <- 0
  for (i in 1:n_points) {
    t <- t_grid[i]
    if (t < 1e-10) next

    S0 <- S_cure(t, pi0, shape0, scale0)
    S1 <- S_cure(t, pi1, shape1, scale1)
    h0 <- h_cure(t, pi0, shape0, scale0)
    h1 <- h_cure(t, pi1, shape1, scale1)

    G_t <- max(0, 1 - t / tau_max)
    w_t <- G_t * 2 * S0 * S1 / (S0 + S1 + 1e-15)

    A_val <- A_val + w_t * (h1 - h0) * dt
  }
  return(A_val)
}

logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))

################################################################################
# Scenario parameters
################################################################################

# Control uncured: Weibull(shape = 2, scale = 1)
k <- 2
gamma0 <- 1

# Treatment uncured: Weibull(shape = 2, scale = gamma1) with HR_u = 0.5
# Under PH Weibull: scale_1 = scale_0 / HR^(1/shape)
HR <- 0.5
gamma1 <- gamma0 / HR^(1 / k)  # = sqrt(2)

# Cure fractions for the three scenarios:
# Scenario A (solid):  pi0 = 0,   pi1 = 0
# Scenario B (dashed): pi0 = 0,   pi1 = 0.2
# Scenario C (dotted): pi0 = 0.2, OR = 1.5 => pi1 = expit(logit(0.2) + log(1.5))
pi0_C <- 0.2
OR_C <- 1.5
pi1_C <- expit(logit(pi0_C) + log(OR_C))

cat("Scenario parameters:\n")
cat("  Control uncured: Weibull(shape =", k, ", scale =", gamma0, ")\n")
cat("  Treatment uncured: Weibull(shape =", k, ", scale =", round(gamma1, 4), ")\n")
cat("  HR_u =", HR, "\n\n")
cat("  Scenario A (solid):  pi0 = 0, pi1 = 0\n")
cat("  Scenario B (dashed): pi0 = 0, pi1 = 0.2\n")
cat("  Scenario C (dotted): pi0 =", pi0_C, ", pi1 =", round(pi1_C, 4),
    " (OR =", OR_C, ")\n\n")

################################################################################
# Compute A(tau) on a dense percentile grid
################################################################################

# Follow-up percentiles: 1% to 99.5%
pct_grid <- seq(1, 99.5, by = 0.5)

# Map percentiles to times on the control uncured Weibull(2,1) scale
time_grid <- qweibull(pct_grid / 100, shape = k, scale = gamma0)

# tau_max for the censoring survival G(t)
tau_max <- max(time_grid) * 1.05

cat("Computing A(tau) for", length(pct_grid), "follow-up percentiles...\n")

# Scenario A: pi0 = 0, pi1 = 0 (no cure)
A_scA <- sapply(time_grid, function(tau) {
  compute_A_tau(tau, tau_max,
                pi0 = 0, shape0 = k, scale0 = gamma0,
                pi1 = 0, shape1 = k, scale1 = gamma1)
})

# Scenario B: pi0 = 0, pi1 = 0.2 (cure in treatment only)
A_scB <- sapply(time_grid, function(tau) {
  compute_A_tau(tau, tau_max,
                pi0 = 0, shape0 = k, scale0 = gamma0,
                pi1 = 0.2, shape1 = k, scale1 = gamma1)
})

# Scenario C: pi0 = 0.2, pi1 from OR = 1.5 (cure in both)
A_scC <- sapply(time_grid, function(tau) {
  compute_A_tau(tau, tau_max,
                pi0 = pi0_C, shape0 = k, scale0 = gamma0,
                pi1 = pi1_C, shape1 = k, scale1 = gamma1)
})

cat("Done.\n\n")

################################################################################
# Find tau* (minimum of the dotted/Scenario C curve)
################################################################################

idx_min_C <- which.min(A_scC)
tau_star_pct <- pct_grid[idx_min_C]

cat("tau* (minimum of dotted curve):\n")
cat("  Percentile:", tau_star_pct, "\n")
cat("  A(tau*):", round(A_scC[idx_min_C], 6), "\n\n")

################################################################################
# Plot Figure 5
################################################################################

plot_df <- data.frame(
  pct = rep(pct_grid, 3),
  A_tau = c(A_scA, A_scB, A_scC),
  scenario = rep(c("A", "B", "C"), each = length(pct_grid))
)

fig5 <- ggplot(plot_df, aes(x = pct, y = A_tau)) +
  geom_line(data = plot_df[plot_df$scenario == "A", ],
            linetype = "solid", linewidth = 1.0, color = "black") +
  geom_line(data = plot_df[plot_df$scenario == "B", ],
            linetype = "dashed", linewidth = 1.0, color = "black") +
  geom_line(data = plot_df[plot_df$scenario == "C", ],
            linetype = "dotted", linewidth = 1.0, color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.4) +
  labs(
    x = "Follow-up (control arm uncured quantile)",
    y = expression(A(tau))
  ) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.minor = element_blank()
  )

ggsave("figure_weighted_diff_threecurves_marked.pdf", fig5, width = 10, height = 6)
ggsave("figure_weighted_diff_threecurves_marked.png", fig5, width = 10, height = 6,
       dpi = 300)

cat("Saved: figure_weighted_diff_threecurves_marked.pdf\n")
cat("Saved: figure_weighted_diff_threecurves_marked.png\n")

################################################################################
# Verification summary
################################################################################

cat("\n=== Verification Summary ===\n")
cat("Scenario A (solid, no cure):\n")
cat("  A(tau) at 50th pct:", round(A_scA[pct_grid == 50], 6), "\n")
cat("  A(tau) at 90th pct:", round(A_scA[pct_grid == 90], 6), "\n")
cat("  Monotone decreasing?", all(diff(A_scA[pct_grid >= 5]) <= 1e-10), "\n\n")

cat("Scenario B (dashed, cure in treatment only):\n")
cat("  A(tau) at 50th pct:", round(A_scB[pct_grid == 50], 6), "\n")
cat("  A(tau) at 90th pct:", round(A_scB[pct_grid == 90], 6), "\n")
cat("  Monotone decreasing?", all(diff(A_scB[pct_grid >= 5]) <= 1e-10), "\n\n")

cat("Scenario C (dotted, cure in both, OR = 1.5):\n")
cat("  A(tau) at 50th pct:", round(A_scC[pct_grid == 50], 6), "\n")
cat("  A(tau) at 90th pct:", round(A_scC[pct_grid == 90], 6), "\n")
cat("  tau* percentile:", tau_star_pct, "\n")
cat("  Non-monotone?", !all(diff(A_scC[pct_grid >= 5]) <= 1e-10), "\n")
