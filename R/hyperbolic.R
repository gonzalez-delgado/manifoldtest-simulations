# ================================================================
# Shared data of the hyperbolic experiments (Sections 5.1 and 5.2).
#
# The benchmark curves (Geo-dCov, BallCov and Raw-dCov under outliers, HHG)
# and the base-point scan of Section 4.3 are computed from the tangent-space
# data written by hyperbolic_shared_data(), so that every procedure sees
# identical samples, in R and in Python.
# ================================================================

# Tangent-space data of one experiment, written to `file` and read back.
#
# The generators and the replication loops reproduce, seed for seed, those of
#   figure2_hyperbolic.R   (experiment "distcov")
#   figure3_outliers.R     (experiment "outlier")
#
# Exported per replication:
#   u : tangent vector of X (contaminated, for the "outlier" experiment)
#   v : tangent vector of Y
# The hyperbolic distance follows from these exactly, via
#   cosh d_H(x,y) = cosh r_x cosh r_y - sinh r_x sinh r_y cos(theta).
#
# `scan` and `xval` optionally restrict the output to one scenario. The data
# are returned as read back from the file (15 significant digits), which is
# what the Python HHG reads too.
hyperbolic_shared_data <- function(experiment = c("distcov", "outlier"),
                                   file = tempfile(fileext = ".csv"),
                                   scan = NULL, xval = NULL) {
  experiment <- match.arg(experiment)
  n_ref <- 5000

  # ---- Hyperbolic geometry (H^3, base p = (1,0,0,0)), from manifoldtest ----
  P0 <- c(1, 0, 0, 0)
  exp_hyperbolic <- function(u) exp_map_hyperbolic(u, P0)
  log_hyperbolic <- function(x) log_map_hyperbolic(x, P0)

  # ---- Experiment 1 "distcov" -- as in figure2_hyperbolic.R ----
  set.seed(777)
  signal_weights <- matrix(runif(9, 0.5, 1.5), 3, 3)
  cross_alpha <- 0.15
  cross_beta  <- 0.10
  signal_fun1 <- function(d) {
    v <- matrix(0, nrow(d), 3)
    for (j in 1:3) {
      v[, j] <- d %*% signal_weights[j, ] +
                cross_alpha * d[, 1] * d[, 2] +
                cross_beta  * d[, 2] * d[, 3]
    }
    v
  }
  gen_u1 <- function(n, seed) {            # tangent vector of X (r ~ U(3,8))
    set.seed(seed)
    r <- runif(n, 3, 8)
    d <- matrix(rnorm(3 * n), n, 3)
    d <- d / sqrt(rowSums(d^2))
    r * d
  }
  gen_v1 <- function(u, noise_level, seed) {   # tangent vector of Y
    set.seed(seed)
    sig <- signal_fun1(u / sqrt(rowSums(u^2)))
    sig + matrix(rnorm(nrow(u) * 3, 0, noise_level), nrow(u), 3)
  }

  # ---- Experiment 2 "outlier" -- as in figure3_outliers.R ----
  set.seed(777)
  W <- matrix(runif(9, 0.5, 1.5), 3, 3)
  signal_fun2 <- function(dir) dir %*% t(W)

  set.seed(42)
  d_ref2 <- matrix(rnorm(3 * n_ref), n_ref, 3)
  d_ref2 <- d_ref2 / sqrt(rowSums(d_ref2^2))
  Var_signal2 <- mean(apply(signal_fun2(d_ref2), 2, stats::var))
  noise_for_dep2 <- function(dep) sqrt(Var_signal2 * (1 / dep - 1))

  gen_u2 <- function(n, seed) {            # tangent vector of X (r ~ U(0.5,2))
    set.seed(seed)
    r <- runif(n, 0.5, 2)
    d <- matrix(rnorm(3 * n), n, 3)
    d <- d / sqrt(rowSums(d^2))
    r * d
  }
  gen_v2 <- function(u, noise_level, seed) {
    set.seed(seed)
    sig <- signal_fun2(u / sqrt(rowSums(u^2)))
    sig + matrix(rnorm(nrow(u) * 3, 0, noise_level), nrow(u), 3)
  }
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

  # ---- Replication list (same scan grids, same seed convention) ----
  n_reps1 <- 100; n_reps2 <- 80
  noise_fixed <- 1.6
  noise_main  <- noise_for_dep2(0.25)
  M_main      <- 50

  reps <- list()
  push <- function(exp, scan, xval, rep, n, noise, M, eps)
    list(exp = exp, scan = scan, xval = xval, rep = rep,
         n = n, noise = noise, M = M, eps = eps)

  if (experiment == "distcov") {
    for (n in c(30, 50, 80, 120, 180, 250))
      for (r in 1:n_reps1) reps[[length(reps) + 1]] <-
        push("distcov", "n", n, r, n, noise_fixed, NA, NA)
    for (nl in c(0.20, 0.40, 0.70, 1.00, 1.50, 2.20, 3.50, 5.00))
      for (r in 1:n_reps1) reps[[length(reps) + 1]] <-
        push("distcov", "dep", nl, r, 150, nl, NA, NA)
  } else {
    grid2 <- list(
      eps = list(vals = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.30)),
      n   = list(vals = c(50, 100, 150, 200, 300)),
      dep = list(vals = c(0.15, 0.20, 0.25, 0.30, 0.40, 0.50)),
      M   = list(vals = c(10, 20, 30, 40, 50, 70)))
    for (sc in names(grid2)) {
      for (xv in grid2[[sc]]$vals) for (r in 1:n_reps2) {
        cfg <- switch(sc,
          eps = list(n = 150, noise = noise_main, M = M_main, eps = xv),
          n   = list(n = xv,  noise = noise_main, M = M_main, eps = 0.10),
          dep = list(n = 150, noise = noise_for_dep2(xv), M = M_main, eps = 0.10),
          M   = list(n = 150, noise = noise_main, M = xv, eps = 0.10))
        reps[[length(reps) + 1]] <-
          push("outlier", sc, xv, r, cfg$n, cfg$noise, cfg$M, cfg$eps)
      }
    }
  }
  keep <- vapply(reps, function(rp)
    (is.null(scan) || rp$scan == scan) &&
    (is.null(xval) || abs(rp$xval - xval) < 1e-9), logical(1))
  reps <- reps[keep]
  if (length(reps) == 0) stop("no replication matches `scan` and `xval`")

  # ---- Generate and stream out ----
  con <- file(file, "w")
  cat("exp,scan,xval,rep,n,i,u1,u2,u3,v1,v2,v3\n", file = con)
  for (rp in reps) {
    if (rp$exp == "distcov") {
      # as in figure2_hyperbolic.R: X = exp(u), then Y from log(X)
      u <- log_hyperbolic(exp_hyperbolic(gen_u1(rp$n, 10000 + rp$rep * 7)))
      v <- gen_v1(u, rp$noise, 10000 + rp$rep * 7 + 100)
    } else {
      u0 <- gen_u2(rp$n, 10000 + rp$rep * 7)
      v  <- gen_v2(u0, rp$noise, 10000 + rp$rep * 7 + 100)
      u  <- contaminate_u(u0, rp$eps, rp$M, 10000 + rp$rep * 7 + 200)
    }
    n <- nrow(u)
    write.table(data.frame(exp = rp$exp, scan = rp$scan, xval = rp$xval,
                           rep = rp$rep, n = n, i = seq_len(n),
                           u1 = u[, 1], u2 = u[, 2], u3 = u[, 3],
                           v1 = v[, 1], v2 = v[, 2], v3 = v[, 3]),
                con, sep = ",", row.names = FALSE, col.names = FALSE)
  }
  close(con)
  d <- read.csv(file)
  attr(d, "file") <- file
  d
}

# Geo dCov p-values: distance covariance fed the intrinsic geodesic distance,
#   dcov.test( d_g(X,X'), d_g(Y,Y'), R = n_perm ),
# for every replication of data written by hyperbolic_shared_data().
geo_dcov_curves <- function(d, n_perm = 199) {
  P0 <- c(1, 0, 0, 0)
  parts <- split(d, list(d$exp, d$scan, d$xval, d$rep), drop = TRUE)
  res <- vector("list", length(parts))
  t0 <- Sys.time()
  for (k in seq_along(parts)) {
    s <- parts[[k]]
    X <- exp_map_hyperbolic(as.matrix(s[, c("u1", "u2", "u3")]), P0)
    Y <- exp_map_hyperbolic(as.matrix(s[, c("v1", "v2", "v3")]), P0)
    set.seed(10000 + s$rep[1] * 7 + 902)
    p <- NA
    tryCatch({ p <- dcov.test(hyperbolic_dist(X), hyperbolic_dist(Y),
                              R = n_perm)$p.value },
             error = function(e) message("test failed: ", conditionMessage(e)))
    res[[k]] <- data.frame(exp = s$exp[1], scan = s$scan[1], xval = s$xval[1],
                           rep = s$rep[1], n = nrow(X), Geo_dCov = p)
    if (k %% 250 == 0)
      cat(sprintf("  Geo dCov %d/%d  (%.1f min)\n", k, length(parts),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  do.call(rbind, res)
}
