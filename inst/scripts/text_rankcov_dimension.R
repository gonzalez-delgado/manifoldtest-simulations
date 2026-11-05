# ================================================================
# RankCov experiments on hyperbolic space and on the sphere.
#
#   DistCov and RankCov (= dCov of the center-outward ranks), both run with
#   manifoldtest::manifold.dcov.test at the base point of the design.
#
# Two manifolds, run separately: hyperbolic space H^d (hyperboloid model)
# and the sphere S^d (base point = north pole).  For each:
#   1. Power vs sample size   (fixed Dep, varying n)
#   2. Power vs dependence    (fixed n, varying Dep = R^2)
#   3. Power vs dimension     (fixed n and Dep, varying d)
#
# Dependence strength is the explained-variance ratio
#   Dep = Var(signal) / (Var(signal) + noise^2), estimated from a fixed
#   reference sample.
#
# Output (in the working directory):
#   rank_{Hyperbolic,Sphere}_power_vs_{n,dep,dim}.csv
#   Rank_{Hyperbolic,Sphere}_Power_vs_{n,Dep,Dim}.png
# ================================================================

suppressMessages(library(manifoldtest.simulations))
suppressMessages({library(ggplot2); library(reshape2)})

# ================================================================
# 1. Manifold geometry (dimension-parameterised)
# ================================================================

#' Hyperbolic geometry on H^d (hyperboloid model, base p = (1,0,...,0)).
make_hyperbolic <- function() {
  list(
    exp = function(u) exp_map_hyperbolic(u),            # at (1, 0, ..., 0)
    log = function(x) log_map_hyperbolic(x, c(1, rep(0, ncol(x) - 1))),
    base = function(x) c(1, rep(0, ncol(x) - 1)),
    manifold = "hyperbolic",
    r_low = 3, r_high = 8,   # large radius -> saturation effect
    label = "Hyperbolic"
  )
}

#' Sphere geometry on S^d (base point = north pole (0,...,0,1)).
#' max_radius: the sphere exponential map is injective only for ||u|| < pi;
#' large noise can push ||v|| beyond pi and break the logarithmic map, so
#' we clamp v to the ball of radius max_radius before exponentiating.
make_sphere <- function() {
  north <- function(d) c(rep(0, d), 1)
  list(
    exp = function(u) exp_map_sphere(u, north(ncol(u))),
    log = function(x) log_map_sphere(x, north(ncol(x) - 1)),
    base = function(x) north(ncol(x) - 1),
    manifold = manifold_sphere(),
    r_low = 0.2, r_high = 1.2,   # X within injectivity radius
    max_radius = pi - 0.1,       # clamp ||v|| to keep log map well-defined
    label = "Sphere"
  )
}

# ================================================================
# 3. Data generation (directional dependence in tangent space)
# ================================================================

# Fixed signal weights (one set per manifold, seeded)
make_signal_weights <- function(d, seed) {
  set.seed(seed)
  matrix(runif(d * d, 0.5, 1.5), d, d)   # d signals x d coords
}

#' Deterministic signal as a function of the unit direction d (dim p).
signal_fun <- function(dir, W) {
  dir %*% t(W)                            # n x p, each col = dir . w_j
}

#' Generate X on the manifold at prescribed radius.
gen_X <- function(geo, n, d, seed) {
  set.seed(seed)
  r <- runif(n, geo$r_low, geo$r_high)
  dir <- matrix(rnorm(d * n), n, d); dir <- dir / sqrt(rowSums(dir^2))
  u <- r * dir
  geo$exp(u)
}

#' Generate dependent Y on the manifold from X (directional signal).
#' If geo$max_radius is set (sphere), v is clamped to the injectivity ball.
gen_dep_Y <- function(geo, X, W, noise_level, seed) {
  set.seed(seed)
  u <- geo$log(X)
  dir <- u / sqrt(rowSums(u^2))
  signal <- signal_fun(dir, W)
  v <- signal + matrix(rnorm(nrow(X) * ncol(W), 0, noise_level), nrow(X), ncol(W))
  if (!is.null(geo$max_radius)) {
    vn <- sqrt(rowSums(v^2))
    too_large <- vn > geo$max_radius
    if (any(too_large)) v[too_large, ] <- v[too_large, ] *
      (geo$max_radius / vn[too_large])
  }
  geo$exp(v)
}

# ================================================================
# 4. Reference explained-variance ratio (Dep)
# ================================================================

ref_dep <- function(geo, d, seed = 42, n_ref = 5000) {
  set.seed(seed)
  dir <- matrix(rnorm(d * n_ref), n_ref, d); dir <- dir / sqrt(rowSums(dir^2))
  W <- make_signal_weights(d, seed = 777)
  signal <- signal_fun(dir, W)
  Var_signal <- mean(apply(signal, 2, var))
  list(Var_signal = Var_signal, W = W,
       dep = function(noise) Var_signal / (Var_signal + noise^2))
}

# ================================================================
# 5. Test functions
# ================================================================

#' Two independence tests on manifold data, through manifoldtest:
#'   (1) DistCov : manifold.dcov.test(center_outward = FALSE)
#'   (2) RankCov : manifold.dcov.test(center_outward = TRUE)
#' @return p-values (DistCov, RankCov)
run_tests <- function(geo, X, Y, n_perm = 199) {
  p <- c(NA, NA)
  tryCatch({ p[1] <- paper_test(X, Y, geo$manifold, px_base = geo$base(X),
                                py_base = geo$base(Y), R = n_perm) },
           error = function(e) message("test failed: ", conditionMessage(e)))
  tryCatch({ p[2] <- paper_test(X, Y, geo$manifold, px_base = geo$base(X),
                                py_base = geo$base(Y), R = n_perm,
                                center_outward = TRUE) },
           error = function(e) message("test failed: ", conditionMessage(e)))
  p
}

#' Power analysis for one scenario.
run_power <- function(geo, d, W, noise_level, n_replicates, n, n_perm = 199) {
  s <- c(0, 0)
  for (r in 1:n_replicates) {
    seed0 <- 10000 + r * 7
    X <- gen_X(geo, n, d, seed0)
    Y <- gen_dep_Y(geo, X, W, noise_level, seed0 + 100)
    pv <- run_tests(geo, X, Y, n_perm)
    if (!is.na(pv[1]) && pv[1] < 0.05) s[1] <- s[1] + 1
    if (!is.na(pv[2]) && pv[2] < 0.05) s[2] <- s[2] + 1
  }
  s / n_replicates * 100
}

# ================================================================
# 6. Full experiment driver for one manifold
# ================================================================

run_manifold <- function(geo, d0 = 3, n_reps = 80, n_perm = 199,
                         out_tag) {
  cat("\n========================================\n")
  cat(sprintf("%s manifold, dimension d = %d\n", geo$label, d0))
  cat("========================================\n")

  # Reference R^2 in the base dimension
  ref <- ref_dep(geo, d0)
  W <- ref$W
  cat(sprintf("Var(signal) = %.4f\n", ref$Var_signal))

  # ---- Experiment A: Power vs Sample Size ----
  cat("\n[A] Power vs Sample Size\n")
  noise_fixed <- 1.5
  dep_fixed <- ref$dep(noise_fixed)
  n_values <- c(30, 50, 80, 120, 180, 250)
  cat(sprintf("Fixed noise = %.2f, Dep = %.3f\n", noise_fixed, dep_fixed))
  power_n <- data.frame()
  for (n in n_values) {
    cat(sprintf("  n = %d ...\n", n))
    pw <- run_power(geo, d0, W, noise_fixed, n_reps, n, n_perm)
    power_n <- rbind(power_n, data.frame(n = n, DistCov = pw[1], RankCov = pw[2]))
    cat(sprintf("    DistCov=%.1f%% RankCov=%.1f%%\n", pw[1], pw[2]))
  }

  # ---- Experiment B: Power vs Dependence ----
  cat("\n[B] Power vs Dependence\n")
  n_fixed <- 150
  noise_values <- c(0.20, 0.40, 0.70, 1.00, 1.50, 2.20, 3.50, 5.00)
  dep_values <- sapply(noise_values, ref$dep)
  power_dep <- data.frame()
  for (i in seq_along(noise_values)) {
    cat(sprintf("  noise = %.2f (Dep = %.3f) ...\n", noise_values[i], dep_values[i]))
    pw <- run_power(geo, d0, W, noise_values[i], n_reps, n_fixed, n_perm)
    power_dep <- rbind(power_dep, data.frame(
      noise = noise_values[i], Dep = dep_values[i],
      DistCov = pw[1], RankCov = pw[2]))
    cat(sprintf("    DistCov=%.1f%% RankCov=%.1f%%\n", pw[1], pw[2]))
  }

  # ---- Experiment C: Power vs Dimension ----
  cat("\n[C] Power vs Dimension\n")
  d_values <- c(2, 3, 4, 5)
  n_dim <- 150
  noise_dim <- 2.5     # Dep ~ 0.15; weak enough that higher d dilutes power
  power_dim <- data.frame()
  for (dd in d_values) {
    cat(sprintf("  d = %d ...\n", dd))
    ref_dd <- ref_dep(geo, dd)
    W_dd <- ref_dd$W
    pw <- run_power(geo, dd, W_dd, noise_dim, n_reps, n_dim, n_perm)
    power_dim <- rbind(power_dim, data.frame(
      d = dd, Dep = ref_dd$dep(noise_dim),
      DistCov = pw[1], RankCov = pw[2]))
    cat(sprintf("    DistCov=%.1f%% RankCov=%.1f%% (Dep=%.3f)\n", pw[1], pw[2],
                ref_dd$dep(noise_dim)))
  }

  # ---- Save data ----
  write.csv(power_n, paste0("rank_", out_tag, "_power_vs_n.csv"), row.names = FALSE)
  write.csv(power_dep, paste0("rank_", out_tag, "_power_vs_dep.csv"), row.names = FALSE)
  write.csv(power_dim, paste0("rank_", out_tag, "_power_vs_dim.csv"), row.names = FALSE)
  cat("Data saved.\n")

  list(power_n = power_n, power_dep = power_dep, power_dim = power_dim,
       dep_fixed = dep_fixed, noise_fixed = noise_fixed, n_fixed = n_fixed)
}

# ================================================================
# 7. Plots
# ================================================================

colors <- c("DistCov (Log + dCov)" = "#E41A1C",
            "RankCov (dCov of F)"  = "#984EA3")
lty <- c("DistCov (Log + dCov)" = "solid",
         "RankCov (dCov of F)"  = "dashed")
shp <- c("DistCov (Log + dCov)" = 16,
         "RankCov (dCov of F)"  = 18)

fmt2 <- function(df, idvars) {
  m <- melt(df, id.vars = idvars, variable.name = "Method", value.name = "Power")
  m$Method <- factor(m$Method, levels = c("DistCov", "RankCov"),
                     labels = c("DistCov (Log + dCov)", "RankCov (dCov of F)"))
  m
}

theme_power <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.box = "horizontal",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
        panel.grid.minor = element_blank())

plot_manifold <- function(res, geo, out_tag, d0) {
  # A: power vs n
  df_n <- fmt2(res$power_n, "n")
  p1 <- ggplot(df_n, aes(x = n, y = Power, color = Method, linetype = Method, shape = Method)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
    scale_color_manual(values = colors) + scale_linetype_manual(values = lty) +
    scale_shape_manual(values = shp) +
    scale_x_continuous(breaks = res$power_n$n) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
    geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.4, alpha = 0.5) +
    labs(title = paste0("Power vs Sample Size on ", geo$label),
         subtitle = sprintf("Fixed Dep = %.2f, d = %d", res$dep_fixed, d0),
         x = "Sample Size (n)", y = "Power (%)") + theme_power
  ggsave(paste0("Rank_", out_tag, "_Power_vs_n.png"), p1, width = 8, height = 6.5, dpi = 150)

  # B: power vs dep
  df_d <- fmt2(res$power_dep, c("noise", "Dep"))
  p2 <- ggplot(df_d, aes(x = Dep, y = Power, color = Method, linetype = Method, shape = Method)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
    scale_color_manual(values = colors) + scale_linetype_manual(values = lty) +
    scale_shape_manual(values = shp) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
    geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.4, alpha = 0.5) +
    labs(title = paste0("Power vs Dependence on ", geo$label),
         subtitle = sprintf("Fixed n = %d, d = %d", res$n_fixed, d0),
         x = expression("Dependence Strength " * (R^2)), y = "Power (%)") + theme_power
  ggsave(paste0("Rank_", out_tag, "_Power_vs_Dep.png"), p2, width = 8, height = 6.5, dpi = 150)

  # C: power vs dimension
  df_dim <- fmt2(res$power_dim, c("d", "Dep"))
  p3 <- ggplot(df_dim, aes(x = d, y = Power, color = Method, linetype = Method, shape = Method)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
    scale_color_manual(values = colors) + scale_linetype_manual(values = lty) +
    scale_shape_manual(values = shp) +
    scale_x_continuous(breaks = res$power_dim$d) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
    geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.4, alpha = 0.5) +
    labs(title = paste0("Power vs Dimension on ", geo$label),
         subtitle = sprintf("Fixed n = %d, Dep ~ %.2f", res$n_fixed, res$power_dim$Dep[1]),
         x = "Manifold Dimension (d)", y = "Power (%)") + theme_power
  ggsave(paste0("Rank_", out_tag, "_Power_vs_Dim.png"), p3, width = 8, height = 6.5, dpi = 150)

  cat("Plots saved.\n")
}

# ================================================================
# 8. Main
# ================================================================
n_reps  <- 80
n_perm  <- 199

# ---- Hyperbolic ----
res_hyp <- run_manifold(make_hyperbolic(), d0 = 3, n_reps, n_perm, "Hyperbolic")
plot_manifold(res_hyp, make_hyperbolic(), "Hyperbolic", 3)

# ---- Sphere ----
res_sph <- run_manifold(make_sphere(), d0 = 3, n_reps, n_perm, "Sphere")
plot_manifold(res_sph, make_sphere(), "Sphere", 3)

cat("\nAll done!\n")
