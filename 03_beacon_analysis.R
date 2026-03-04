################################################################################
# 03_beacon_analysis.R
#
# Reproduces: Main text Section 4 — Data Example (Figure 6)
#             Supplementary Material Section C — BEACON-Immuno Trial Analysis
#               C1: Hypothesis test p-values (Table)
#               C2: AIC for mixture cure models by arm (Table)
#               C3: Parameter estimates for best-fitting models (Table)
#               Fitted survival function figure
#
# Uses reconstructed individual patient-level data from the BEACON-Immuno
# trial (NCT02308527), comparing dinutuximab beta (dB) plus chemotherapy
# versus chemotherapy alone for relapsed high-risk neuroblastoma.
#
# Data file: ipd_pfs_jco.csv (reconstructed from published KM curves)
#
# Outputs:
#   figure_beacon_final.pdf  (Main text Figure 6: KM + fitted curves + A(t))
#   figure_S1_survival.pdf   (Supplement: fitted survival curves only)
#   Console output: hypothesis test table, AIC table, parameter estimates
#
# Required packages: survival, flexsurv, flexsurvcure, ggplot2, gridExtra
################################################################################

library(survival)
library(flexsurv)
library(flexsurvcure)
library(ggplot2)
library(gridExtra)

# Save all console output to a text file (in addition to printing)
sink("beacon_results.txt", split = TRUE)

################################################################################
# Load data
################################################################################

dat <- read.csv("ipd_pfs_jco.csv", stringsAsFactors = FALSE)
dat$arm_factor <- factor(dat$arm, levels = c("No dB", "dB"))

dat_dB <- dat[dat$arm == "dB", ]
dat_nodB <- dat[dat$arm == "No dB", ]

cat("=== BEACON-Immuno Trial Data ===\n\n")
cat("Total patients:", nrow(dat), "\n")
cat("dB arm: N =", nrow(dat_dB), ", Events =", sum(dat_dB$event),
    "(", round(100 * mean(dat_dB$event), 1), "%)\n")
cat("No dB arm: N =", nrow(dat_nodB), ", Events =", sum(dat_nodB$event),
    "(", round(100 * mean(dat_nodB$event), 1), "%)\n\n")

################################################################################
# SECTION C1: Hypothesis tests
#
# We implement the Fleming-Harrington weighted log-rank test manually to
# support arbitrary (rho, gamma) combinations.
################################################################################

# General weighted log-rank test using Fleming-Harrington G^{rho,gamma} weights
weighted_logrank <- function(time, status, group, rho = 0, gamma = 0) {
  event_times <- sort(unique(time[status == 1]))
  if (length(event_times) == 0) return(list(chisq = 0, pval = 1))

  U <- 0  # test statistic numerator
  V <- 0  # variance
  S_km <- 1.0  # pooled KM estimate (S(t-))

  for (tj in event_times) {
    Y  <- sum(time >= tj)                     # total at risk
    Y1 <- sum(time >= tj & group == 1)        # group 1 at risk
    Y0 <- Y - Y1
    d  <- sum(time == tj & status == 1)       # total events
    d1 <- sum(time == tj & status == 1 & group == 1)

    # FH weight: w = S(t-)^rho * (1 - S(t-))^gamma
    w <- (S_km^rho) * ((1 - S_km)^gamma)

    if (Y > 1) {
      U <- U + w * (d1 - Y1 * d / Y)
      V <- V + w^2 * Y1 * Y0 * d * (Y - d) / (Y^2 * (Y - 1))
    }
    # Update pooled KM
    if (Y > 0) S_km <- S_km * (1 - d / Y)
  }

  if (V > 0) {
    chisq <- U^2 / V
    pval <- 1 - pchisq(chisq, df = 1)
  } else {
    chisq <- 0; pval <- 1
  }
  list(chisq = chisq, pval = pval)
}

cat("--- Hypothesis Tests (Supplement C1) ---\n\n")

group_indicator <- as.integer(dat$arm_factor == "dB")

# 1. Standard log-rank (rho = 0, gamma = 0)
lr <- weighted_logrank(dat$time, dat$event, group_indicator, rho = 0, gamma = 0)

# 2. FH(1,0) early-weighted
fh10 <- weighted_logrank(dat$time, dat$event, group_indicator, rho = 1, gamma = 0)

# 3. FH(0,1) late-weighted
fh01 <- weighted_logrank(dat$time, dat$event, group_indicator, rho = 0, gamma = 1)

# 4. FH(1,1)
fh11 <- weighted_logrank(dat$time, dat$event, group_indicator, rho = 1, gamma = 1)

cat(sprintf("  %-30s  p = %.3f\n", "Standard log-rank", lr$pval))
cat(sprintf("  %-30s  p = %.3f\n", "FH(1,0) early-weighted", fh10$pval))
cat(sprintf("  %-30s  p = %.3f\n", "FH(0,1) late-weighted", fh01$pval))
cat(sprintf("  %-30s  p = %.3f\n", "FH(1,1) optimal", fh11$pval))

################################################################################
# SECTION C2: AIC for mixture cure models by arm
#
# Fit 5 latency distributions per arm using flexsurvcure with identity link
# (so the theta parameter directly equals the cure fraction).
################################################################################

cat("\n--- Mixture Cure Model Fitting (Supplement C2) ---\n")

distributions <- c("exp", "weibull", "gamma", "lnorm", "gengamma")
dist_names <- c("Exponential", "Weibull", "Gamma", "Log-normal", "Gen. Gamma")

fit_cure_models <- function(data, arm_name) {
  cat("\n", arm_name, ":\n")
  tau <- max(data$time)

  results <- data.frame(
    Distribution = dist_names,
    AIC = NA_real_,
    pi_n = NA_real_,
    r_n = NA_real_,
    stringsAsFactors = FALSE
  )
  models <- list()

  for (i in seq_along(distributions)) {
    tryCatch({
      fit <- flexsurvcure(Surv(time, event) ~ 1, data = data,
                          dist = distributions[i], link = "identity")

      pi_n <- fit$res["theta", "est"]
      S_tau <- summary(fit, t = tau, type = "survival")[[1]]$est
      S_uc_tau <- max(0, (S_tau - pi_n) / (1 - pi_n))
      r_n <- ifelse(S_tau > 0, S_uc_tau / S_tau, 0)

      results$AIC[i] <- AIC(fit)
      results$pi_n[i] <- pi_n
      results$r_n[i] <- r_n
      models[[distributions[i]]] <- fit
    }, error = function(e) {
      cat("    ", dist_names[i], ": FAILED -", e$message, "\n")
    })
  }

  best_idx <- which.min(results$AIC)
  cat("  Best model:", dist_names[best_idx],
      "(AIC =", round(results$AIC[best_idx], 2), ")\n")
  cat("  pi_n =", round(results$pi_n[best_idx], 3),
      ", r_n =", round(results$r_n[best_idx], 3), "\n")

  # RECeUS-AIC decision
  cure_ok <- (results$pi_n[best_idx] > 0.025) & (results$r_n[best_idx] < 0.05)
  cat("  RECeUS criteria met:", cure_ok, "\n")

  list(results = results, models = models, best_idx = best_idx,
       pi_n = results$pi_n[best_idx], r_n = results$r_n[best_idx],
       cure_ok = cure_ok)
}

dB_fits <- fit_cure_models(dat_dB, "Dinutuximab beta (dB)")
nodB_fits <- fit_cure_models(dat_nodB, "Chemotherapy alone (No dB)")

# Print AIC table
cat("\n--- AIC Table (Supplement C2) ---\n")
aic_table <- data.frame(
  Distribution = dist_names,
  AIC_dB = round(dB_fits$results$AIC, 2),
  AIC_nodB = round(nodB_fits$results$AIC, 2)
)
print(aic_table, row.names = FALSE)

################################################################################
# SECTION C3: Parameter estimates for best-fitting models
################################################################################

cat("\n--- Parameter Estimates (Supplement C3) ---\n")

best_dB_dist <- distributions[dB_fits$best_idx]
best_nodB_dist <- distributions[nodB_fits$best_idx]
best_model_dB <- dB_fits$models[[best_dB_dist]]
best_model_nodB <- nodB_fits$models[[best_nodB_dist]]

cat("\ndB arm (best distribution:", dist_names[dB_fits$best_idx], "):\n")
print(best_model_dB$res[, c("est", "L95%", "U95%", "se")])

cat("\nNo dB arm (best distribution:", dist_names[nodB_fits$best_idx], "):\n")
print(best_model_nodB$res[, c("est", "L95%", "U95%", "se")])

################################################################################
# Mixture Cure Model LRT (arm effects on cure fraction AND latency)
#
# Tests H0: S_0(t) = S_1(t) for all t
# vs    H1: S_0(t) != S_1(t) for some t
#
# For each distribution: fit null (no arm effect) and alt (arm on cure + latency),
# then compute the LRT statistic.
################################################################################

cat("\n--- Mixture Cure Model LRT ---\n")

get_anc_for_dist <- function(dist) {
  switch(dist,
    "exp"      = list(rate = ~arm_factor),
    "weibull"  = list(scale = ~arm_factor),
    "gamma"    = list(rate = ~arm_factor),
    "lnorm"    = list(meanlog = ~arm_factor),
    "gengamma" = list(mu = ~arm_factor),
    NULL
  )
}

lrt_results <- data.frame(
  Distribution = dist_names,
  AIC_null = NA_real_,
  AIC_alt = NA_real_,
  LRT_stat = NA_real_,
  df = NA_integer_,
  p_value = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(distributions)) {
  d <- distributions[i]
  anc_spec <- get_anc_for_dist(d)
  tryCatch({
    fit_null <- flexsurvcure(Surv(time, event) ~ 1, data = dat,
                              dist = d, mixture = TRUE)
    fit_alt <- flexsurvcure(Surv(time, event) ~ arm_factor, anc = anc_spec,
                             data = dat, dist = d, mixture = TRUE)

    lr_stat <- max(2 * (fit_alt$loglik - fit_null$loglik), 0)
    df_diff <- length(coef(fit_alt)) - length(coef(fit_null))
    p_val <- 1 - pchisq(lr_stat, df = df_diff)

    lrt_results$AIC_null[i] <- AIC(fit_null)
    lrt_results$AIC_alt[i] <- AIC(fit_alt)
    lrt_results$LRT_stat[i] <- lr_stat
    lrt_results$df[i] <- df_diff
    lrt_results$p_value[i] <- p_val

    cat(sprintf("  %-12s: LRT = %6.3f, df = %d, p = %.4f (AIC_alt = %.1f)\n",
                dist_names[i], lr_stat, df_diff, p_val, AIC(fit_alt)))
  }, error = function(e) {
    cat(sprintf("  %-12s: FAILED - %s\n", dist_names[i], e$message))
  })
}

# Select best MCM LRT by AIC of the alternative model (best-fitting joint model)
valid_lrt <- lrt_results[!is.na(lrt_results$p_value), ]
if (nrow(valid_lrt) > 0) {
  best_lrt_idx <- which.min(valid_lrt$AIC_alt)
  mcm_lrt_pval <- valid_lrt$p_value[best_lrt_idx]
  cat("\n  Best MCM LRT (by AIC of alt model):",
      valid_lrt$Distribution[best_lrt_idx],
      ", p =", round(mcm_lrt_pval, 3), "\n")
}

# Complete hypothesis test table (Supplement C1)
cat("\n\n=== Complete Hypothesis Test Table (Supplement C1) ===\n")
cat(sprintf("  %-30s  p-value\n", "Method"))
cat("  ", paste(rep("-", 45), collapse = ""), "\n")
cat(sprintf("  %-30s  %.3f\n", "Standard log-rank", lr$pval))
cat(sprintf("  %-30s  %.3f\n", "FH(1,0) early-weighted", fh10$pval))
cat(sprintf("  %-30s  %.3f\n", "FH(0,1) late-weighted", fh01$pval))
cat(sprintf("  %-30s  %.3f\n", "FH(1,1) optimal", fh11$pval))
if (nrow(valid_lrt) > 0) {
  cat(sprintf("  %-30s  %.3f\n", "Mixture cure model LRT", mcm_lrt_pval))
}

################################################################################
# MAIN TEXT FIGURE 6: KM plot with overlaid cure model fits + A(t)
################################################################################

cat("\n--- Generating Figure 6 ---\n")

# Time grid for fitted curves
t_max <- max(dat$time)
t_grid <- seq(0.01, t_max, length.out = 200)

# Fitted survival from best per-arm models
S_dB <- sapply(t_grid, function(t)
  summary(best_model_dB, t = t, type = "survival")[[1]]$est)
S_nodB <- sapply(t_grid, function(t)
  summary(best_model_nodB, t = t, type = "survival")[[1]]$est)

# Fitted hazard for A(tau) computation
t_fine <- seq(0.01, t_max, length.out = 500)
h_dB_fine <- sapply(t_fine, function(t)
  summary(best_model_dB, t = t, type = "hazard")[[1]]$est)
h_nodB_fine <- sapply(t_fine, function(t)
  summary(best_model_nodB, t = t, type = "hazard")[[1]]$est)
S_dB_fine <- sapply(t_fine, function(t)
  summary(best_model_dB, t = t, type = "survival")[[1]]$est)
S_nodB_fine <- sapply(t_fine, function(t)
  summary(best_model_nodB, t = t, type = "survival")[[1]]$est)

# Compute A(tau) from fitted models
# w(t) = G(t) * 2*S_0(t)*S_1(t) / (S_0(t) + S_1(t))
# G(t) = max(0, 1 - t/tau_max)
# tau_max set to maximum observed time
tau_max_data <- max(dat$time)
A_tau_values <- numeric(length(t_grid))

for (j in seq_along(t_grid)) {
  tau_j <- t_grid[j]
  idx <- which(t_fine <= tau_j)
  if (length(idx) < 2) { A_tau_values[j] <- 0; next }

  t_sub <- t_fine[idx]
  dt <- diff(t_sub)
  A_val <- 0
  for (ii in seq_along(dt)) {
    S0_i <- S_nodB_fine[idx[ii]]   # Control = chemotherapy alone
    S1_i <- S_dB_fine[idx[ii]]     # Treatment = dB
    h0_i <- h_nodB_fine[idx[ii]]
    h1_i <- h_dB_fine[idx[ii]]
    G_i <- max(0, 1 - t_sub[ii] / tau_max_data)
    w_i <- G_i * 2 * S0_i * S1_i / (S0_i + S1_i + 1e-15)
    A_val <- A_val + w_i * (h1_i - h0_i) * dt[ii]
  }
  A_tau_values[j] <- A_val
}

cat("  A(t) at t = 9 months:",
    round(A_tau_values[which.min(abs(t_grid - 9))], 4), "\n")
cat("  A(t) at t = 24 months:",
    round(A_tau_values[which.min(abs(t_grid - 24))], 4), "\n")
cat("  Minimum A(t):", round(min(A_tau_values), 4),
    "at t =", round(t_grid[which.min(A_tau_values)], 1), "months\n")

# KM curves
km_fit <- survfit(Surv(time, event) ~ arm_factor, data = dat)

create_step_data <- function(km_fit, strata_idx) {
  if (strata_idx == 1) {
    idx <- 1:km_fit$strata[1]
  } else {
    idx <- (km_fit$strata[1] + 1):length(km_fit$time)
  }
  data.frame(time = c(0, km_fit$time[idx]), surv = c(1, km_fit$surv[idx]))
}

km_nodB <- create_step_data(km_fit, 1)
km_dB <- create_step_data(km_fit, 2)

# Fitted survival data for overlay
fitted_surv_df <- data.frame(
  time = c(t_grid, t_grid),
  surv = c(S_dB, S_nodB),
  arm = rep(c("dB", "nodB"), each = length(t_grid))
)

# Panel A: KM + fitted mixture cure survival
# KM: thick lines; Fitted: thin gray lines
# dB = solid, No dB = dashed
p1 <- ggplot() +
  geom_step(data = km_dB, aes(x = time, y = surv),
            linetype = "solid", linewidth = 0.9, color = "black") +
  geom_step(data = km_nodB, aes(x = time, y = surv),
            linetype = "dashed", linewidth = 0.9, color = "black") +
  geom_line(data = fitted_surv_df[fitted_surv_df$arm == "dB", ],
            aes(x = time, y = surv),
            linetype = "solid", linewidth = 0.6, color = "gray40") +
  geom_line(data = fitted_surv_df[fitted_surv_df$arm == "nodB", ],
            aes(x = time, y = surv),
            linetype = "dashed", linewidth = 0.6, color = "gray40") +
  labs(x = "Time (months)", y = "Progression-Free Survival", title = "(A)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

# Panel B: Weighted average hazard difference A(t)
plot_data <- data.frame(time = t_grid, A_tau = A_tau_values)
p2 <- ggplot(plot_data, aes(x = time, y = A_tau)) +
  geom_line(linewidth = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(x = "Time (months)", y = expression(A(t)), title = "(B)") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave("figure_beacon_final.pdf", combined_plot, width = 14, height = 6)
cat("Saved: figure_beacon_final.pdf\n")

################################################################################
# SUPPLEMENT FIGURE: Fitted survival curves only
################################################################################

p_surv <- ggplot() +
  geom_line(data = fitted_surv_df[fitted_surv_df$arm == "dB", ],
            aes(x = time, y = surv),
            linetype = "solid", linewidth = 0.9, color = "black") +
  geom_line(data = fitted_surv_df[fitted_surv_df$arm == "nodB", ],
            aes(x = time, y = surv),
            linetype = "dashed", linewidth = 0.9, color = "black") +
  labs(x = "Time (months)", y = "S(t)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

ggsave("figure_S1_survival.pdf", p_surv, width = 7, height = 5)
cat("Saved: figure_S1_survival.pdf\n")

