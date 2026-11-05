# ================================================================
# Figure 2 (Section 5.1): comparison on hyperbolic space H^3.
#
# Data:  X = exp_p(u) on the hyperboloid model, tangent radius U(3,8);
#        Y built in the tangent space from the direction of u.
#        Dependence strength is the explained-variance ratio
#        Dep = Var(signal) / (Var(signal) + noise^2), taken from a fixed
#        reference sample of 5000 directions.
#
# Methods compared, 199 permutations and 100 replications:
#   (a) DistCov : manifoldtest::manifold.dcov.test, base point p = (1,0,0,0)
#   (b) BallCov : hyperbolic distance matrices        -> Ball::bcov.test
#   (c) Raw dCov: ambient R^4 hyperboloid coordinates -> energy::dcov.test
#   (d) Geo dCov: hyperbolic distance matrices        -> energy::dcov.test
#   (e) HHG     : hyperbolic distance matrices        -> HHG (Python)
# (d) and (e) are computed from the samples written by
# hyperbolic_shared_data("distcov"), which are those of (a)-(c).
#
# Scans:  n = 30, 50, 80, 120, 180, 250   (Dep = 0.31)
#         Dep = 0.97 ... 0.04             (n = 150)
#
# Output (in the working directory):
#   Fig_Hyperbolic_Power_vs_n.png, Fig_Hyperbolic_Power_vs_Dep.png
#   hyperbolic_power_vs_n.csv, hyperbolic_power_vs_dep.csv
# ================================================================

suppressMessages({library(manifoldtest.simulations); library(Ball)})

# ================================================================
# 1. Hyperbolic geometry (hyperboloid model H^3), from manifoldtest
# ================================================================

P0 <- c(1, 0, 0, 0)                                  # base point of the design
exp_hyperbolic <- function(u) exp_map_hyperbolic(u, P0)
log_hyperbolic <- function(x) log_map_hyperbolic(x, P0)

# ================================================================
# 2. Data generation (directional dependence in tangent space)
# ================================================================

# Fixed signal weights (generated once, used for all replications)
set.seed(777)
signal_weights <- matrix(runif(9, 0.5, 1.5), 3, 3)   # 3 signals x 3 coords
cross_alpha <- 0.15
cross_beta  <- 0.10

#' Deterministic signal as a function of the unit direction d.
#' The dependence is DIRECTIONAL: it does not depend on the radius.
#' @param d unit direction matrix (n x 3)
#' @return signal matrix (n x 3)
signal_fun <- function(d) {
  v <- matrix(0, nrow(d), 3)
  for (j in 1:3) {
    v[, j] <- d %*% signal_weights[j, ] +
              cross_alpha * d[, 1] * d[, 2] +
              cross_beta  * d[, 2] * d[, 3]
  }
  v
}

#' Generate X on H^3 with LARGE hyperbolic radius (near the boundary).
#' @param n    sample size
#' @param seed random seed
#' @return hyperboloid points (n x 4)
gen_X <- function(n, seed) {
  set.seed(seed)
  r <- runif(n, 3, 8)                       # large radius -> saturation
  d <- matrix(rnorm(3 * n), n, 3)
  d <- d / sqrt(rowSums(d^2))               # uniform on S^2
  u <- r * d
  exp_hyperbolic(u)
}

#' Generate dependent Y on H^3 from X.
#' v = signal(direction of u) + noise,  Y = exp_p(v).
#' @param X           hyperboloid points (n x 4)
#' @param noise_level sd of the additive noise on v
#' @param seed        random seed
#' @return hyperboloid points (n x 4)
gen_dep_Y <- function(X, noise_level, seed) {
  set.seed(seed)
  u <- log_hyperbolic(X)
  d <- u / sqrt(rowSums(u^2))
  signal <- signal_fun(d)
  v <- signal + matrix(rnorm(nrow(X) * 3, 0, noise_level), nrow(X), 3)
  exp_hyperbolic(v)
}

# ================================================================
# 3. Reference explained-variance ratio (Dep)
# ================================================================

# Fixed large reference sample (independent of the test sample size)
set.seed(42)
n_ref <- 5000
u_ref <- log_hyperbolic(gen_X(n_ref, seed = 42))
d_ref <- u_ref / sqrt(rowSums(u_ref^2))
signal_ref <- signal_fun(d_ref)
Var_signal <- mean(apply(signal_ref, 2, var))
cat(sprintf("Reference Var(signal) = %.4f\n", Var_signal))

#' Map a noise level to the explained-variance ratio Dep.
dep_from_noise <- function(noise) Var_signal / (Var_signal + noise^2)

# ================================================================
# 4. Test functions
# ================================================================

#' Three independence tests on hyperbolic data.
#' @return p-values: (DistCov, BallCov, Raw dCov)
run_tests <- function(X, Y, n_perm = 199, seed0) {
  p <- c(NA, NA, NA)
  # (a) DistCov: log map + dCov, through manifoldtest
  tryCatch({
    p[1] <- paper_test(X, Y, "hyperbolic", px_base = P0, R = n_perm)
  }, error = function(e) message("test failed: ", conditionMessage(e)))
  # (b) BallCov: hyperbolic distance + bcov
  tryCatch({
    Dx <- hyperbolic_dist(X); Dy <- hyperbolic_dist(Y)
    rb <- bcov.test(Dx, Dy, num.permutations = n_perm); p[2] <- rb$p.value
  }, error = function(e) message("test failed: ", conditionMessage(e)))
  # (c) Raw dCov: ambient R^4 hyperboloid coordinates + dCov.
  # Ball::bcov.test leaves the global random stream in a non-reproducible
  # state, so the stream is re-seeded before the permutations of Raw dCov
  set.seed(seed0 + 901)
  tryCatch({
    rr <- dcov.test(X, Y, R = n_perm); p[3] <- rr$p.value
  }, error = function(e) message("test failed: ", conditionMessage(e)))
  p
}

#' Power analysis for one scenario.
#' @return vector of power: (DistCov, BallCov, Raw dCov) in percent
run_power <- function(n_replicates, n, noise_level, n_perm = 199) {
  s <- c(0, 0, 0)
  for (r in 1:n_replicates) {
    seed0 <- 10000 + r * 7
    X <- gen_X(n, seed0)
    Y <- gen_dep_Y(X, noise_level, seed0 + 100)
    pv <- run_tests(X, Y, n_perm, seed0)
    if (!is.na(pv[1]) && pv[1] < 0.05) s[1] <- s[1] + 1
    if (!is.na(pv[2]) && pv[2] < 0.05) s[2] <- s[2] + 1
    if (!is.na(pv[3]) && pv[3] < 0.05) s[3] <- s[3] + 1
  }
  s / n_replicates * 100
}

# ================================================================
# 5. Main experiments
# ================================================================

cat("\n========================================\n")
cat("Hyperbolic H^3 comparison: DistCov vs BallCov vs Raw dCov\n")
cat("========================================\n")

n_reps <- 100
n_perm <- 199

# ---- Experiment 1: Power vs Sample Size (fixed dependence) ----
# Fixed at a moderate dependence where the methods differ clearly
# (smoke test: at noise = 1.6, Dep ~ 0.31, DistCov ~ 100% vs BallCov ~ 23%)
cat("\nExperiment 1: Power vs Sample Size\n")
noise_fixed <- 1.6
dep_fixed  <- dep_from_noise(noise_fixed)
n_values   <- c(30, 50, 80, 120, 180, 250)
cat(sprintf("Fixed noise = %.2f, Dep = %.3f\n", noise_fixed, dep_fixed))

power_n <- data.frame()
for (n in n_values) {
  cat(sprintf("  n = %d ...\n", n))
  pw <- run_power(n_reps, n, noise_fixed, n_perm)
  power_n <- rbind(power_n, data.frame(
    n = n, DistCov = pw[1], BallCov = pw[2], Raw_dCov = pw[3]))
  cat(sprintf("    DistCov=%.1f%% BallCov=%.1f%% Raw=%.1f%%\n", pw[1], pw[2], pw[3]))
}

# ---- Experiment 2: Power vs Dependence (fixed n, varying Dep) ----
# Noise levels: Dep roughly evenly spread on a linear axis
cat("\nExperiment 2: Power vs Dependence\n")
n_fixed <- 150
noise_values <- c(0.20, 0.40, 0.70, 1.00, 1.50, 2.20, 3.50, 5.00)
dep_values   <- sapply(noise_values, dep_from_noise)
cat("Noise -> Dep mapping:\n")
print(data.frame(noise = noise_values, Dep = round(dep_values, 3)), row.names = FALSE)

power_dep <- data.frame()
for (i in seq_along(noise_values)) {
  nl <- noise_values[i]
  cat(sprintf("  noise = %.2f (Dep = %.3f) ...\n", nl, dep_values[i]))
  pw <- run_power(n_reps, n_fixed, nl, n_perm)
  power_dep <- rbind(power_dep, data.frame(
    noise = nl, Dep = dep_values[i],
    DistCov = pw[1], BallCov = pw[2], Raw_dCov = pw[3]))
  cat(sprintf("    DistCov=%.1f%% BallCov=%.1f%% Raw=%.1f%%\n", pw[1], pw[2], pw[3]))
}

# ================================================================
# 6. Geo dCov and HHG, on the shared samples
# ================================================================

cat("\nGeo dCov and HHG\n")
shared <- hyperbolic_shared_data("distcov")
gi <- geo_dcov_curves(shared, n_perm)
hhg <- run_hhg(attr(shared, "file"), "hyperbolic", n_perm)

rate <- function(p) 100 * mean(p < 0.05)
hhg_pow <- aggregate(cbind(HHG = p_value) ~ scan + xval, data = hhg, FUN = rate)
gi_pow  <- aggregate(Geo_dCov ~ scan + xval, data = gi, FUN = rate)

power_n <- merge(power_n, hhg_pow[hhg_pow$scan == "n", c("xval", "HHG")],
                 by.x = "n", by.y = "xval")
power_n <- merge(power_n, gi_pow[gi_pow$scan == "n", c("xval", "Geo_dCov")],
                 by.x = "n", by.y = "xval")
power_dep <- merge(power_dep, hhg_pow[hhg_pow$scan == "dep", c("xval", "HHG")],
                   by.x = "noise", by.y = "xval")
power_dep <- merge(power_dep, gi_pow[gi_pow$scan == "dep", c("xval", "Geo_dCov")],
                   by.x = "noise", by.y = "xval")
power_dep <- power_dep[order(power_dep$Dep, decreasing = TRUE), ]

write.csv(power_n, "hyperbolic_power_vs_n.csv", row.names = FALSE)
write.csv(power_dep, "hyperbolic_power_vs_dep.csv", row.names = FALSE)

# ================================================================
# 7. Figure 2 (two panels sharing a page -> legend on the first only)
# ================================================================

METHODS <- c("DistCov", "BallCov", "HHG", "Raw_dCov", "Geo_dCov")
p1 <- power_plot(power_n, "n", "n", METHODS,
                 "Power vs Sample Size on Hyperbolic Space",
                 "H^3,  fixed Dep = 0.31",
                 "Sample size N", TRUE, xbreaks = power_n$n)
save_power_plot(p1, "Fig_Hyperbolic_Power_vs_n.png")

p2 <- power_plot(power_dep, c("noise", "Dep"), "Dep", METHODS,
                 "Power vs Dependence on Hyperbolic Space",
                 "H^3,  fixed N = 150",
                 expression("Dependence strength " * R^2), FALSE)
save_power_plot(p2, "Fig_Hyperbolic_Power_vs_Dep.png")

cat("\n===== Summary =====\n")
cat("--- Power vs Sample Size (Dep =", round(dep_fixed, 3), ") ---\n")
print(power_n, row.names = FALSE)
cat("\n--- Power vs Dependence ---\n")
print(power_dep, row.names = FALSE)
cat("\nDone!\n")
