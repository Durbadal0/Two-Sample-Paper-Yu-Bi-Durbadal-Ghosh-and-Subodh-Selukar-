################################################################################
# 02_supplement_scenarios.R
#
# Reproduces: Supplementary Material Section B — Simulation Scenarios
#
# Generates 3-panel figures for five theoretical scenarios showing:
#   Panel A: Survival functions S(t) for control and treatment arms
#   Panel B: Hazard functions h(t) for both arms
#   Panel C: Weighted average hazard difference A(tau) and cumulative hazard
#            difference Delta H(t) = H_1(t) - H_0(t)
#
# Control uncured distribution: Weibull(shape = 2, scale = 1)
#
# Scenarios:
#   1. No cure (pi0 = pi1 = 0), HR_u = 0.5
#   2. Cure in treatment only (pi0 = 0, pi1 = 0.2), HR_u = 0.5
#   3. Cure in both, identical uncured (pi0 = 0.2, OR = 1.5), HR_u = 1
#   4. Cure in both, different uncured (pi0 = 0.2, OR = 1.5), HR_u = 0.5
#   5. Equal cure fractions, different uncured (pi0 = pi1 = 0.2), HR_u = 0.5
#
# Outputs:
#   figure_S_scenario1.pdf through figure_S_scenario5.pdf
#
# Required packages: ggplot2, gridExtra, grid
################################################################################

library(ggplot2)
library(gridExtra)
library(grid)

################################################################################
# Helper functions
################################################################################

logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))

# Weibull survival
S_weibull <- function(t, shape, scale) {
  exp(-(t / scale)^shape)
}

# Weibull hazard
h_weibull <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1)
}

# Mixture cure survival: S(t) = pi + (1 - pi) * S_u(t)
S_cure <- function(t, pi_g, shape, scale) {
  pi_g + (1 - pi_g) * S_weibull(t, shape, scale)
}

# Mixture cure hazard: h(t) = (1 - pi) * f_u(t) / S(t)
h_cure <- function(t, pi_g, shape, scale) {
  S_uc <- S_weibull(t, shape, scale)
  f_uc <- h_weibull(t, shape, scale) * S_uc  # f_u(t) = h_u(t) * S_u(t)
  S_mix <- pi_g + (1 - pi_g) * S_uc
  (1 - pi_g) * f_uc / S_mix
}

# Compute A(tau) via left Riemann sum
compute_A_tau <- function(tau_eval, tau_max,
                          pi0, shape0, scale0,
                          pi1, shape1, scale1) {
  if (tau_eval <= 0.1) return(0)
  n_points <- 200
  t_grid <- seq(0.01, tau_eval, length.out = n_points)
  dt <- tau_eval / n_points

  A_val <- 0
  for (t in t_grid) {
    S0 <- S_cure(t, pi0, shape0, scale0)
    S1 <- S_cure(t, pi1, shape1, scale1)
    h0 <- h_cure(t, pi0, shape0, scale0)
    h1 <- h_cure(t, pi1, shape1, scale1)

    G_t <- max(0, 1 - t / tau_max)
    w_t <- G_t * 2 * S0 * S1 / (S0 + S1 + 1e-10)

    A_val <- A_val + w_t * (h1 - h0) * dt
  }
  return(A_val)
}

################################################################################
# Function to create 3-panel figure for a scenario
#
# Line styles:
#   Panel A & B: Treatment = solid, Control = dashed
#   Panel C: A(tau) = solid, Delta H(t) = dashed
################################################################################

create_scenario_figure <- function(pi0, pi1,
                                   shape0, scale0,
                                   shape1, scale1,
                                   tau_max = 5) {
  t_grid <- seq(0.01, tau_max, length.out = 200)
  tau_grid <- seq(0.5, tau_max, length.out = 50)

  # Compute survival and hazard
  S0 <- sapply(t_grid, function(t) S_cure(t, pi0, shape0, scale0))
  S1 <- sapply(t_grid, function(t) S_cure(t, pi1, shape1, scale1))
  h0 <- sapply(t_grid, function(t) h_cure(t, pi0, shape0, scale0))
  h1 <- sapply(t_grid, function(t) h_cure(t, pi1, shape1, scale1))

  # Cumulative hazard difference: Delta H(t) = H_1(t) - H_0(t)
  H0 <- -log(pmax(S0, 1e-10))
  H1 <- -log(pmax(S1, 1e-10))
  delta_H <- H1 - H0

  # Weighted average hazard difference A(tau)
  A_tau <- sapply(tau_grid, function(tau) {
    compute_A_tau(tau, tau_max, pi0, shape0, scale0, pi1, shape1, scale1)
  })

  # Panel A: Survival functions
  pA <- ggplot() +
    geom_line(aes(x = t_grid, y = S1), linetype = "solid",
              linewidth = 0.8, color = "black") +
    geom_line(aes(x = t_grid, y = S0), linetype = "dashed",
              linewidth = 0.8, color = "black") +
    labs(x = "t", y = "S(t)", title = "(A)") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16))

  # Panel B: Hazard functions
  pB <- ggplot() +
    geom_line(aes(x = t_grid, y = h1), linetype = "solid",
              linewidth = 0.8, color = "black") +
    geom_line(aes(x = t_grid, y = h0), linetype = "dashed",
              linewidth = 0.8, color = "black") +
    labs(x = "t", y = "h(t)", title = "(B)") +
    theme_bw() +
    theme(plot.title = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16))

  # Panel C: A(tau) and Delta H(t)
  pC <- ggplot() +
    geom_line(aes(x = tau_grid, y = A_tau), linetype = "solid",
              linewidth = 0.8, color = "black") +
    geom_line(aes(x = t_grid, y = delta_H), linetype = "dashed",
              linewidth = 0.8, color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    labs(x = "t", y = "Value", title = "(C)") +
    theme_bw() +
    theme(plot.title = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16))

  grid.arrange(pA, pB, pC, ncol = 3)
}

################################################################################
# Generate scenario figures
################################################################################

# Control parameters: Weibull(shape = 2, scale = 1)
shape0 <- 2
scale0 <- 1

# Treatment scale for HR_u = 0.5: scale_1 = scale_0 / HR^(1/shape)
scale1_HR05 <- scale0 / (0.5^(1 / shape0))  # = sqrt(2)

# pi1 when pi0 = 0.2 and OR = 1.5
pi1_OR15 <- expit(logit(0.2) + log(1.5))  # approx 0.2727

cat("Generating Supplementary Section B figures...\n\n")

# --- Scenario 1: No cure, HR_u = 0.5 ---
cat("  Scenario 1: No cure (pi0 = pi1 = 0), HR_u = 0.5\n")
fig1 <- create_scenario_figure(
  pi0 = 0, pi1 = 0,
  shape0 = shape0, scale0 = scale0,
  shape1 = shape0, scale1 = scale1_HR05,
  tau_max = 3
)
ggsave("figure_S_scenario1.pdf", fig1, width = 14, height = 5)

# --- Scenario 2: Cure in treatment only, HR_u = 0.5 ---
cat("  Scenario 2: Cure in treatment only (pi0 = 0, pi1 = 0.2), HR_u = 0.5\n")
fig2 <- create_scenario_figure(
  pi0 = 0, pi1 = 0.2,
  shape0 = shape0, scale0 = scale0,
  shape1 = shape0, scale1 = scale1_HR05,
  tau_max = 4
)
ggsave("figure_S_scenario2.pdf", fig2, width = 14, height = 5)

# --- Scenario 3: Cure in both, identical uncured (OR = 1.5, HR_u = 1) ---
cat("  Scenario 3: Cure in both (pi0 = 0.2, OR = 1.5), identical uncured\n")
fig3 <- create_scenario_figure(
  pi0 = 0.2, pi1 = pi1_OR15,
  shape0 = shape0, scale0 = scale0,
  shape1 = shape0, scale1 = scale0,  # HR_u = 1 => same scale
  tau_max = 5
)
ggsave("figure_S_scenario3.pdf", fig3, width = 14, height = 5)

# --- Scenario 4: Cure in both, different uncured (OR = 1.5, HR_u = 0.5) ---
cat("  Scenario 4: Cure in both (pi0 = 0.2, OR = 1.5), HR_u = 0.5\n")
fig4 <- create_scenario_figure(
  pi0 = 0.2, pi1 = pi1_OR15,
  shape0 = shape0, scale0 = scale0,
  shape1 = shape0, scale1 = scale1_HR05,
  tau_max = 5
)
ggsave("figure_S_scenario4.pdf", fig4, width = 14, height = 5)

# --- Scenario 5: Equal cure fractions, different uncured (HR_u = 0.5) ---
cat("  Scenario 5: Same cure (pi0 = pi1 = 0.2), HR_u = 0.5\n")
fig5 <- create_scenario_figure(
  pi0 = 0.2, pi1 = 0.2,
  shape0 = shape0, scale0 = scale0,
  shape1 = shape0, scale1 = scale1_HR05,
  tau_max = 5
)
ggsave("figure_S_scenario5.pdf", fig5, width = 14, height = 5)

cat("\nAll Section B figures saved.\n")
