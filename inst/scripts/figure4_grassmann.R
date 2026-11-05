# ================================================================
# Figure 4 (Section 5.3): power on the Grassmann manifold Gr(4,2).
#
# Geometry: points are n x k Stiefel matrices U with U'U = I_k; tangent
# vectors are horizontal lifts (U0'Delta = 0); the geodesic distance is
# sqrt(sum_i theta_i^2), where cos(theta_i) are the singular values of U'V.
#
#   base point    U0, grassmann_basepoint() (inst/extdata)
#   response curve gamma(theta) = U0 cos(theta) + Q sin(theta),
#                 Q an n x k matrix of orthonormal columns with U0'Q = 0;
#                 arc length L = pi*sqrt(k), every point in the cut domain
#                 of U0 except the one at arc length L/2
#   X             exp_{U0}(r d),  r ~ U(0.2, 1.4),  d uniform on the unit
#                 sphere of the tangent space (dimension k(n-k) = 4)
#   signal        a(X) = pi * (rank{<Delta, w>} - 1/2) / N,  Delta = log_{U0} X,
#                 w a fixed direction
#   response      theta = (a(X) + (pi/2) B) mod pi,  B ~ Bernoulli(1/2),
#                 Y = gamma(theta)
#
# On the curve the intrinsic distance is the circle distance, so
#     d(z1,z2) + d(z1, z2 + L/2) = L/2   for all z1, z2 on the curve,
# and averaging over the four combinations of (B, B') gives
#     E[d_g(Y,Y') | X, X'] = L/4    for every pair X, X'.
#
# Scans:  N = 100, 150, 200, 300, 500   (dependence held fixed)
#         null                          (Y independent of X, N = 200)
#
# Tests, 199 permutations and 100 replications:
#   DistCov  : manifoldtest::manifold.dcov.test on manifold_grassmann(4, 2),
#              base point U0
#   Raw dCov : dcov.test( U, V ), U and V the Stiefel matrices themselves
#   BallCov  : bcov.test( d_g(X,X'), d_g(Y,Y') )
#   Geo dCov : dcov.test( d_g(X,X'), d_g(Y,Y') )
#   HHG      : HHG test on d_g (Python)
#
# The points are generated with manifoldtest::exp_map_grassmann and written
# to a temporary file as Stiefel representatives (x1..x8 = as.vector(U) of
# X, y1..y8 of Y), which both the R tests and the Python HHG read. The R
# tests are spread over all but two of the available cores.
#
# Output (in the working directory):
#   Fig_Grassmann_Power_vs_N.png
#   grassmann_power.csv   the numbers behind the figure
#   grassmann_size.csv    rejection rates of the null scan
# ================================================================
suppressMessages({
  library(manifoldtest.simulations)
  library(Ball); library(data.table)
})

N_DIM <- 4; K_DIM <- 2; DIM <- N_DIM * K_DIM
N_REPS <- 100
R_LO <- 0.2; R_HI <- 1.4

U0 <- grassmann_basepoint()
horiz <- function(M) M - U0 %*% (t(U0) %*% M)
unit_norm <- function(M) { s <- sqrt(sum(M^2)); if (s < 1e-14) matrix(0, nrow(M), ncol(M)) else M / s }

# Tangent coordinates at U0 of a horizontal lift Delta (n x k, U0' Delta = 0),
# in the basis of the orthogonal complement of U0 used by manifoldtest
U0_perp <- qr.Q(qr(U0), complete = TRUE)[, -seq_len(K_DIM), drop = FALSE]
lift_to_coords <- function(D) t(apply(D, 1, function(row)
  as.vector(crossprod(U0_perp, matrix(row, N_DIM, K_DIM)))))

# ---- the closed geodesic ----------------------------------------------
set.seed(999)
Q <- qr.Q(qr(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))))
set.seed(4242)
W_SIG <- rnorm(DIM)          # fixed direction that defines the signal
gamma <- function(theta) {
  if (length(theta) == 1) theta <- rep(theta, K_DIM)
  U0 %*% diag(cos(theta), K_DIM, K_DIM) + Q %*% diag(sin(theta), K_DIM, K_DIM)
}

# ---- data of one replication --------------------------------------------
gen_block <- function(N, seed, dependent = TRUE) {
  set.seed(seed)
  Dx <- matrix(0, N, DIM)
  for (i in 1:N)
    Dx[i, ] <- as.vector(runif(1, R_LO, R_HI) * unit_norm(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))))
  if (dependent) {
    a <- pi * (rank(as.vector(Dx %*% W_SIG)) - 0.5) / N
    B <- rbinom(N, 1, 0.5)
    theta <- (a + (pi / 2) * B) %% pi
  } else {
    theta <- runif(N, 0, pi)
  }
  X <- exp_map_grassmann(lift_to_coords(Dx), U0)         # points of Gr(4,2)
  # points on the curve, stored as the Stiefel representative of
  # exp_U0(log_U0(gamma(theta))), as in the pipeline of the paper (Raw-dCov,
  # which works on the Stiefel entries, depends on the representative)
  G <- t(sapply(theta, function(t) as.vector(gamma(t))))
  Y <- exp_map_grassmann(log_map_grassmann(G, U0, N_DIM, K_DIM), U0)
  list(X = X, Y = Y, theta = theta)
}

# ---- replicate list ----------------------------------------------------
reps <- list()
n_vals <- c(100, 150, 200, 300, 500)
for (N in n_vals) for (r in 1:N_REPS) reps[[length(reps) + 1]] <-
  list(scan = "n", xval = N, N = N, r = r, dep_flag = TRUE)
for (r in 1:N_REPS) reps[[length(reps) + 1]] <-
  list(scan = "null", xval = 0, N = 200, r = r, dep_flag = FALSE)
cat("replications:", length(reps), "\n")

# ---- generate and stream out ------------------------------------------
INPUT <- tempfile("grassmann_input", fileext = ".csv")
con <- file(INPUT, "w")
cols <- c(paste0("x", 1:DIM), paste0("y", 1:DIM))
cat(paste(c("scan", "xval", "rep", "n", "i", cols), collapse = ","), "\n",
    sep = "", file = con)
t0 <- Sys.time()
for (k in seq_along(reps)) {
  rp <- reps[[k]]
  blk <- gen_block(rp$N, 70000 + rp$r * 7, rp$dep_flag)
  for (i in 1:rp$N)
    cat(paste(c(rp$scan, rp$xval, rp$r, rp$N, i,
                as.vector(blk$X[i, ]), as.vector(blk$Y[i, ])), collapse = ","),
        "\n", sep = "", file = con)
  if (k %% 50 == 0)
    cat(sprintf("  %d/%d  (%.1f min)\n", k, length(reps),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
close(con)
cat(sprintf("data generated in %.1f min.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- the R and Python pipelines build the same distances --------------
NR <- 40
d40 <- fread(INPUT, nrows = NR)
D_R <- grassmann_dist(as.matrix(d40[, paste0("x", 1:DIM), with = FALSE]), N_DIM, K_DIM)
D_py <- python_grassmann_dist(INPUT, NR)
cat(sprintf("R vs Python distance matrices, %d rows: max difference %.3e\n",
            NR, max(abs(D_R - D_py))))
stopifnot(max(abs(D_R - D_py)) < 1e-9)

# ---- the four R tests on one replication -------------------------------
run_rep <- function(X, Y, rep_id) {
  lab <- c("DistCov", "Raw_dCov", "BallCov", "Geo_dCov")
  p <- setNames(rep(NA_real_, 4), lab)
  Dx <- as.dist(grassmann_dist(X, N_DIM, K_DIM))
  Dy <- as.dist(grassmann_dist(Y, N_DIM, K_DIM))
  seed0 <- 10000 + rep_id * 7
  set.seed(seed0 + 900)
  tryCatch(p["DistCov"]  <- paper_test(X, Y, GR, px_base = as.vector(U0), R = N_PERM),
           error = function(e) message("test failed: ", conditionMessage(e)))
  set.seed(seed0 + 901)
  tryCatch(p["Raw_dCov"] <- dcov.test(X, Y, R = N_PERM)$p.value,
           error = function(e) message("test failed: ", conditionMessage(e)))
  # bcov.test re-seeds itself and ignores the global stream, so its seed is
  # passed explicitly for the permutations to differ across replications
  tryCatch(p["BallCov"]  <- bcov.test(Dx, Dy, num.permutations = N_PERM,
                                      seed = seed0 + 903)$p.value,
           error = function(e) message("test failed: ", conditionMessage(e)))
  set.seed(seed0 + 904)
  tryCatch(p["Geo_dCov"] <- dcov.test(Dx, Dy, R = N_PERM)$p.value,
           error = function(e) message("test failed: ", conditionMessage(e)))
  p
}

# ---- DistCov, Raw dCov, BallCov, Geo dCov -------------------------------
N_PERM <- 199
N_CORES <- if (.Platform$OS.type == "windows") 1L else
  max(1L, parallel::detectCores() - 2L)
GR <- manifold_grassmann(N_DIM, K_DIM)

d <- fread(INPUT)
XCOL <- paste0("x", 1:DIM); YCOL <- paste0("y", 1:DIM)
idx <- split(seq_len(nrow(d)), paste(d$scan, d$xval, d$rep, sep = "|"))
cat(sprintf("R tests: %d replications on %d cores\n", length(idx), N_CORES))
t0 <- Sys.time()
res <- parallel::mclapply(names(idx), function(key) {
  s <- d[idx[[key]]]
  p <- run_rep(as.matrix(s[, ..XCOL]), as.matrix(s[, ..YCOL]), s$rep[1])
  data.frame(scan = s$scan[1], xval = s$xval[1], rep = s$rep[1],
             n = nrow(s), as.list(p))
}, mc.cores = N_CORES)
cur <- rbindlist(res)
cat(sprintf("R tests done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- HHG ---------------------------------------------------------------
hhg <- as.data.table(run_hhg(INPUT, "grassmann", N_PERM))

# ---- power and figure --------------------------------------------------

METHODS <- c("DistCov", "BallCov", "HHG", "Raw_dCov", "Geo_dCov")

d <- merge(cur, hhg[, .(scan, xval, rep, HHG)], by = c("scan", "xval", "rep"))

pow <- d[, lapply(.SD, function(p) round(100 * mean(p < 0.05), 1)),
         .SDcols = METHODS, by = .(scan, xval)]
fwrite(pow, "grassmann_power.csv")

cat("replications per scenario:\n"); print(d[, .N, by = .(scan, xval)][order(scan, xval)])
cat("\nnull scan (Y independent of X), rejection rate at the 5% level:\n")
print(pow[scan == "null"], row.names = FALSE)
fwrite(pow[scan == "null"], "grassmann_size.csv")

# ---- figure ------------------------------------------------------------
pn <- pow[scan == "n"][order(xval)]
# printed at 0.8\textwidth, with the legend on one line
p <- power_plot(as.data.frame(pn), "xval", "xval", METHODS,
                "Power on the Grassmann manifold", NULL, "Sample size N",
                xbreaks = c(100, 200, 300, 500), legend_rows = 1)
save_power_plot(p, "Fig_Grassmann_Power_vs_N.png", w = 5.6, h = 3.5)
print(pn, row.names = FALSE)
