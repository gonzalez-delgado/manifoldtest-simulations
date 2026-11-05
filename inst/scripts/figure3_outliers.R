# ================================================================
# Figure 3 (Section 5.2): outliers in the X-margin, on H^3.
#
# Data:  as in figure2_hyperbolic.R, with tangent radius U(0.5,2) and the
#        cross terms switched off (a = b = 0), so that the signal is the
#        linear term W d alone.
#
# Outlier model: a fraction eta of the X-margin is replaced in the tangent
#        space, u_i <- A * w_i with w_i uniform on S^2; the Y-margin v is
#        left as generated.
#
# Methods compared, 199 permutations and 80 replications:
#   DistCov  : manifoldtest::manifold.dcov.test(center_outward = FALSE)
#   RankCov  : manifoldtest::manifold.dcov.test(center_outward = TRUE), i.e.
#              distance covariance on the center-outward ranks
#   (both with the base point p = (1,0,0,0) of the design)
#   BallCov  : bcov.test( d_g(X,X'), d_g(Y,Y') )
#   Raw dCov : dcov.test( X, Y ), ambient R^4 coordinates
#   Geo dCov : dcov.test( d_g(X,X'), d_g(Y,Y') )
#   HHG      : HHG test on d_g (Python)
# The last four are computed from the samples written by
# hyperbolic_shared_data("outlier"), which are those of DistCov and RankCov.
#
# Sweeps:
#   eta = 0.00 ... 0.30   (n = 150, Dep = 0.25, A = 50)
#   n   = 50 ... 300      (eta = 0.10, Dep = 0.25, A = 50)
#   Dep = 0.15 ... 0.50   (n = 150, eta = 0.10, A = 50)
#   A   = 10 ... 70       (n = 150, eta = 0.10, Dep = 0.25)
#
# Output (in the working directory):
#   Fig_Outlier_Power_vs_eps.png, Fig_Outlier_Power_vs_n.png,
#   Fig_Outlier_Power_vs_Dep.png, Fig_Outlier_Power_vs_M.png,
#   outlier_power_vs_{eps,n,dep,M}.csv
# ================================================================

suppressMessages({library(manifoldtest.simulations); library(Ball)})

# ---- Hyperbolic geometry on H^3 (base p = (1,0,0,0)), from manifoldtest ----
P0 <- c(1, 0, 0, 0)
exp_hyperbolic <- function(u) exp_map_hyperbolic(u, P0)
log_hyperbolic <- function(x) log_map_hyperbolic(x, P0)

# ---- Signal weights & data generation (moderate radius bulk) ----
set.seed(777)
W <- matrix(runif(9, 0.5, 1.5), 3, 3)
signal_fun <- function(dir) dir %*% t(W)

gen_clean <- function(n, seed) {
  set.seed(seed)
  # bulk: moderate radius r ~ U(0.5, 2)
  r <- runif(n, 0.5, 2)
  d <- matrix(rnorm(3 * n), n, 3); d <- d / sqrt(rowSums(d^2))
  u <- r * d
  X <- exp_hyperbolic(u)
  list(X = X, u = u, d = d)
}

gen_dep_v <- function(u, noise_level, seed) {
  set.seed(seed)
  d <- u / sqrt(rowSums(u^2))
  signal <- signal_fun(d)
  v <- signal + matrix(rnorm(nrow(u) * 3, 0, noise_level), nrow(u), 3)
  list(v = v, Y = exp_hyperbolic(v))
}

# ---- Reference explained variance ----
set.seed(42)
n_ref <- 5000
d_ref <- matrix(rnorm(3 * n_ref), n_ref, 3); d_ref <- d_ref / sqrt(rowSums(d_ref^2))
Var_signal <- mean(apply(signal_fun(d_ref), 2, var))
cat(sprintf("Reference Var(signal) = %.4f\n", Var_signal))
noise_for_dep <- function(dep) sqrt(Var_signal * (1 / dep - 1))
dep_for_noise <- function(noise) Var_signal / (Var_signal + noise^2)

# ---- Contaminate the X-margin (u) with outliers ----
contaminate_u <- function(u, frac, M, seed) {
  set.seed(seed)
  n <- nrow(u); d <- ncol(u)
  n_out <- round(frac * n)
  idx <- if (n_out > 0) sample(1:n, n_out) else integer(0)
  if (n_out > 0) {
    w <- matrix(rnorm(n_out * d), n_out, d)
    w <- w / sqrt(rowSums(w^2))
    u[idx, ] <- M * w
  }
  u
}

# ---- Tests: DistCov vs RankCov, through manifoldtest ----
# u, v are the tangent coordinates at p; the data on H^3 are exp_p(u), exp_p(v)
run_tests <- function(u, v, n_perm = 199) {
  X <- exp_hyperbolic(u); Y <- exp_hyperbolic(v)
  p <- c(NA, NA)
  tryCatch({ p[1] <- paper_test(X, Y, "hyperbolic", px_base = P0, R = n_perm) },
           error = function(e) message("test failed: ", conditionMessage(e)))
  tryCatch({ p[2] <- paper_test(X, Y, "hyperbolic", px_base = P0, R = n_perm,
                                center_outward = TRUE) },
           error = function(e) message("test failed: ", conditionMessage(e)))
  p
}

#' Power for one scenario with a given contamination.
run_power <- function(n, noise_level, frac, M, n_replicates, n_perm = 199) {
  s <- c(0, 0)
  for (r in 1:n_replicates) {
    seed0 <- 10000 + r * 7
    cl <- gen_clean(n, seed0)
    dv <- gen_dep_v(cl$u, noise_level, seed0 + 100)
    u_cont <- contaminate_u(cl$u, frac, M, seed0 + 200)
    pv <- run_tests(u_cont, dv$v, n_perm)
    if (!is.na(pv[1]) && pv[1] < 0.05) s[1] <- s[1] + 1
    if (!is.na(pv[2]) && pv[2] < 0.05) s[2] <- s[2] + 1
  }
  s / n_replicates * 100
}

# ================================================================
# Sweeps
# ================================================================
n_reps <- 80
n_perm <- 199

cat("\n========================================\n")
cat("Rank vs plain dCov under outliers, H^3\n")
cat("========================================\n")

# ---- (1) Power vs outlier fraction (main) ----
# Weaker dependence + extreme outliers expose the plain-dCov failure.
cat("\n[1] Power vs outlier fraction eps\n")
n_main  <- 150
dep_main <- 0.25
noise_main <- noise_for_dep(dep_main)
M_main  <- 50
eps_values <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.30)
cat(sprintf("n=%d, Dep=%.2f, M=%d\n", n_main, dep_main, M_main))
pow_eps <- data.frame()
for (ep in eps_values) {
  cat(sprintf("  eps = %.2f ...\n", ep))
  pw <- run_power(n_main, noise_main, ep, M_main, n_reps, n_perm)
  pow_eps <- rbind(pow_eps, data.frame(eps = ep, DistCov = pw[1], RankCov = pw[2]))
  cat(sprintf("    DistCov=%.1f%% RankCov=%.1f%%\n", pw[1], pw[2]))
}

# ---- (2) Power vs sample size ----
cat("\n[2] Power vs sample size n\n")
eps_n  <- 0.10
n_values <- c(50, 100, 150, 200, 300)
pow_n <- data.frame()
for (nn in n_values) {
  cat(sprintf("  n = %d ...\n", nn))
  pw <- run_power(nn, noise_main, eps_n, M_main, n_reps, n_perm)
  pow_n <- rbind(pow_n, data.frame(n = nn, DistCov = pw[1], RankCov = pw[2]))
  cat(sprintf("    DistCov=%.1f%% RankCov=%.1f%%\n", pw[1], pw[2]))
}

# ---- (3) Power vs dependence ----
cat("\n[3] Power vs dependence Dep\n")
n_dep <- 150
eps_dep <- 0.10
dep_values <- c(0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
pow_dep <- data.frame()
for (dp in dep_values) {
  cat(sprintf("  Dep = %.2f ...\n", dp))
  nl <- noise_for_dep(dp)
  pw <- run_power(n_dep, nl, eps_dep, M_main, n_reps, n_perm)
  pow_dep <- rbind(pow_dep, data.frame(Dep = dp, DistCov = pw[1], RankCov = pw[2]))
  cat(sprintf("    DistCov=%.1f%% RankCov=%.1f%%\n", pw[1], pw[2]))
}

# ---- (4) Power vs outlier magnitude ----
cat("\n[4] Power vs outlier magnitude M\n")
n_M <- 150
eps_M <- 0.10
M_values <- c(10, 20, 30, 40, 50, 70)
pow_M <- data.frame()
for (Mm in M_values) {
  cat(sprintf("  M = %d ...\n", Mm))
  pw <- run_power(n_M, noise_main, eps_M, Mm, n_reps, n_perm)
  pow_M <- rbind(pow_M, data.frame(M = Mm, DistCov = pw[1], RankCov = pw[2]))
  cat(sprintf("    DistCov=%.1f%% RankCov=%.1f%%\n", pw[1], pw[2]))
}

# ================================================================
# BallCov, Raw dCov, Geo dCov and HHG, on the shared samples
# ================================================================
cat("\nBallCov, Raw dCov, Geo dCov and HHG\n")
shared <- hyperbolic_shared_data("outlier")

# ---- BallCov and Raw dCov ----
parts <- split(shared, list(shared$scan, shared$xval, shared$rep), drop = TRUE)
res <- vector("list", length(parts))
t0 <- Sys.time()
for (k in seq_along(parts)) {
  s <- parts[[k]]
  rp <- s$rep[1]
  X <- exp_hyperbolic(as.matrix(s[, c("u1", "u2", "u3")]))
  Y <- exp_hyperbolic(as.matrix(s[, c("v1", "v2", "v3")]))
  seed0 <- 10000 + rp * 7
  pb <- NA; pr <- NA
  set.seed(seed0 + 900)
  tryCatch({ pb <- bcov.test(hyperbolic_dist(X), hyperbolic_dist(Y),
                             num.permutations = n_perm)$p.value },
           error = function(e) message("test failed: ", conditionMessage(e)))
  set.seed(seed0 + 901)
  tryCatch({ pr <- dcov.test(X, Y, R = n_perm)$p.value },
           error = function(e) message("test failed: ", conditionMessage(e)))
  res[[k]] <- data.frame(scan = s$scan[1], xval = s$xval[1], rep = rp,
                         BallCov = pb, Raw_dCov = pr)
  if (k %% 100 == 0)
    cat(sprintf("  BallCov / Raw dCov %d/%d  (%.1f min)\n", k, length(parts),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
ex <- do.call(rbind, res)

# ---- Geo dCov and HHG ----
gi <- geo_dcov_curves(shared, n_perm)
hhg <- run_hhg(attr(shared, "file"), "hyperbolic", n_perm)

rate <- function(p) 100 * mean(p < 0.05)
ex_pow  <- aggregate(cbind(BallCov, Raw_dCov) ~ scan + xval, data = ex, FUN = rate)
gi_pow  <- aggregate(Geo_dCov ~ scan + xval, data = gi, FUN = rate)
hhg_pow <- aggregate(cbind(HHG = p_value) ~ scan + xval, data = hhg, FUN = rate)

o <- list(eps = pow_eps, n = pow_n, dep = pow_dep, M = pow_M)
KEY <- c(eps = "eps", n = "n", dep = "Dep", M = "M")
for (sc in names(o)) {
  k <- KEY[[sc]]
  df <- o[[sc]]
  df <- merge(df, ex_pow[ex_pow$scan == sc, c("xval", "BallCov", "Raw_dCov")],
              by.x = k, by.y = "xval")
  df <- merge(df, hhg_pow[hhg_pow$scan == sc, c("xval", "HHG")],
              by.x = k, by.y = "xval")
  df <- merge(df, gi_pow[gi_pow$scan == sc, c("xval", "Geo_dCov")],
              by.x = k, by.y = "xval")
  o[[sc]] <- df[order(df[[k]]), ]
  write.csv(o[[sc]], sprintf("outlier_power_vs_%s.csv", sc), row.names = FALSE)
}

# ================================================================
# Figure 3 (four panels sharing a page -> legend on the first only)
# ================================================================
METHODS <- c("DistCov", "BallCov", "HHG", "Raw_dCov", "Geo_dCov", "RankCov")

p3 <- power_plot(o$eps, "eps", "eps", METHODS,
                 "Power vs Outlier Fraction under Outliers (H^3)",
                 expression(paste("N = 150,  Dep = 0.25,  A = 50")),
                 expression("Outlier fraction " * eta), TRUE, xbreaks = o$eps$eps)
save_power_plot(p3, "Fig_Outlier_Power_vs_eps.png")

p4 <- power_plot(o$n, "n", "n", METHODS,
                 "Power vs Sample Size under Outliers (H^3)",
                 expression(paste(eta, " = 0.10,  Dep = 0.25,  A = 50")),
                 "Sample size N", FALSE, xbreaks = o$n$n)
save_power_plot(p4, "Fig_Outlier_Power_vs_n.png")

p5 <- power_plot(o$dep, "Dep", "Dep", METHODS,
                 "Power vs Dependence under Outliers (H^3)",
                 expression(paste("N = 150,  ", eta, " = 0.10,  A = 50")),
                 expression("Dependence strength " * R^2), FALSE, xbreaks = o$dep$Dep)
save_power_plot(p5, "Fig_Outlier_Power_vs_Dep.png")

p6 <- power_plot(o$M, "M", "M", METHODS,
                 "Power vs Outlier Magnitude under Outliers (H^3)",
                 expression(paste("N = 150,  ", eta, " = 0.10,  Dep = 0.25")),
                 "Outlier magnitude A", FALSE, xbreaks = o$M$M)
save_power_plot(p6, "Fig_Outlier_Power_vs_M.png")

cat("\n===== Summary =====\n")
for (sc in names(o)) {
  cat(sprintf("--- Power vs %s ---\n", KEY[[sc]]))
  print(o[[sc]], row.names = FALSE)
}
cat("\nDone!\n")
