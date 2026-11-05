# ================================================================
# Shared helpers, available to every script after
#     library(manifoldtest.simulations)
#
# The tests of the paper (DistCov and RankCov) and the geometry they use
# (logarithmic and exponential maps, center-outward ranks) come from the
# manifoldtest package. This file only adds
#   - the pairwise geodesic distance matrices needed by the benchmark tests
#     (Geo-dCov, BallCov), which the package does not provide;
#   - the unit sphere, used in Section 4.3 and built with
#     manifoldtest::new_manifold();
#   - the fixed base point of the Grassmann experiment.
# ================================================================

# Seed of the quasi-uniform grid directions of the center-outward ranks
# (passed as `seed` to manifold.dcov.test), as in the simulations of the paper
GRID_SEED <- 20260723

# ---- Pairwise geodesic distance matrices (benchmarks only) --------------

# Flat torus [0,1)^d
torus_dist <- function(X) {
  X <- as.matrix(X)
  D <- matrix(0, nrow(X), nrow(X))
  for (j in seq_len(ncol(X))) {
    da <- abs(outer(X[, j], X[, j], "-"))
    D <- D + pmin(da, 1 - da)^2
  }
  as.dist(sqrt(D))
}

# Hyperboloid model: cosh d(x, y) = -<x, y>_L
hyperbolic_dist <- function(X) {
  X <- as.matrix(X)
  M <- X[, -1, drop = FALSE] %*% t(X[, -1, drop = FALSE]) - outer(X[, 1], X[, 1])
  D <- acosh(pmax(-M, 1 + 1e-12))
  diag(D) <- 0
  as.dist(D)
}

# Grassmann manifold Gr(n, k), rows as.vector(U) of n x k Stiefel matrices:
# square root of the sum of the squared principal angles. The eigenvalues of
# the k x k matrices (U_i^T U_j)^T (U_i^T U_j) are computed in closed form,
# for all pairs at once (k = 2 or 3).
grassmann_dist <- function(X, n, k) {
  X <- as.matrix(X)
  N <- nrow(X)
  U <- array(X, c(N, n, k))
  M <- array(0, c(N, N, k, k))
  for (r in 1:k) for (s in 1:k) M[, , r, s] <- U[, , r] %*% t(U[, , s])
  S <- array(0, c(N, N, k, k))
  for (a in 1:k) for (r in 1:k) for (s in 1:k)
    S[, , r, s] <- S[, , r, s] + M[, , a, r] * M[, , a, s]
  th <- acos(pmin(sqrt(pmax(sym_eig(S), 0)), 1))
  D <- sqrt(rowSums(th^2, dims = 2))
  diag(D) <- 0
  D
}

# Eigenvalues, in descending order, of the symmetric k x k matrices stacked in
# the last two dimensions of S (closed forms for k = 2 and k = 3)
sym_eig <- function(S) {
  k <- dim(S)[3]
  if (k == 2) {
    a <- S[, , 1, 1]; b <- S[, , 1, 2]; c <- S[, , 2, 2]
    q <- (a + c) / 2
    h <- sqrt(((a - c) / 2)^2 + b^2)
    out <- array(0, c(dim(a), 2)); out[, , 1] <- q + h; out[, , 2] <- q - h
    return(out)
  }
  if (k == 3) {
    a11 <- S[, , 1, 1]; a22 <- S[, , 2, 2]; a33 <- S[, , 3, 3]
    a12 <- S[, , 1, 2]; a13 <- S[, , 1, 3]; a23 <- S[, , 2, 3]
    q  <- (a11 + a22 + a33) / 3
    p2 <- (a11 - q)^2 + (a22 - q)^2 + (a33 - q)^2 + 2 * (a12^2 + a13^2 + a23^2)
    p  <- sqrt(p2 / 6)
    ps <- ifelse(p > 1e-300, p, 1)
    b11 <- (a11 - q) / ps; b22 <- (a22 - q) / ps; b33 <- (a33 - q) / ps
    b12 <- a12 / ps; b13 <- a13 / ps; b23 <- a23 / ps
    r <- (b11 * b22 * b33 + 2 * b12 * b13 * b23 -
          b11 * b23^2 - b22 * b13^2 - b33 * b12^2) / 2
    phi <- acos(pmin(pmax(r, -1), 1)) / 3
    e1 <- q + 2 * p * cos(phi); e3 <- q + 2 * p * cos(phi + 2 * pi / 3)
    e2 <- 3 * q - e1 - e3
    out <- array(0, c(dim(a11), 3))
    out[, , 1] <- e1; out[, , 2] <- e2; out[, , 3] <- e3
    return(out)
  }
  stop("sym_eig is implemented for k = 2 and k = 3")
}

# Unit sphere S^d in R^(d+1)
sphere_dist <- function(X) {
  X <- as.matrix(X)
  as.dist(acos(pmin(pmax(X %*% t(X), -1), 1)))
}

# ---- The unit sphere as a manifoldtest manifold -------------------------

# Orthonormal basis of the tangent space at p (columns), by Gram-Schmidt on
# the canonical basis projected onto the complement of p. At the north pole
# (0, ..., 0, 1) it is e_1, ..., e_d, so that the tangent coordinates are the
# first d coordinates of the tangent vector.
sphere_tangent_basis <- function(p) {
  m <- length(p)
  B <- NULL
  for (j in seq_len(m)) {
    w <- diag(m)[, j] - p[j] * p
    if (!is.null(B)) w <- w - B %*% crossprod(B, w)
    if (sqrt(sum(w^2)) > 1e-8) B <- cbind(B, w / sqrt(sum(w^2)))
    if (!is.null(B) && ncol(B) == m - 1) break
  }
  B
}

log_map_sphere <- function(X, p) {
  X <- as.matrix(X)
  cos_t <- pmin(pmax(as.vector(X %*% p), -1), 1)
  V <- X - outer(cos_t, p)
  nv <- sqrt(rowSums(V^2))
  V <- V * (acos(cos_t) / ifelse(nv > 1e-12, nv, 1))
  V %*% sphere_tangent_basis(p)
}

exp_map_sphere <- function(V, p) {
  V <- as.matrix(V)
  Vt <- V %*% t(sphere_tangent_basis(p))
  r <- sqrt(rowSums(Vt^2))
  outer(cos(r), p) + Vt * (sin(r) / ifelse(r > 1e-12, r, 1))
}

manifold_sphere <- function() {
  new_manifold(
    name = "sphere",
    representation = "unit vectors in R^(d+1)",
    contains = function(X) all(abs(rowSums(as.matrix(X)^2) - 1) < 1e-8),
    log_map = log_map_sphere,
    dist = function(X, p) acos(pmin(pmax(as.vector(as.matrix(X) %*% p), -1), 1)),
    is_valid = function(X, p) {
      if (length(p) != ncol(X) || abs(sum(p^2) - 1) > 1e-8) stop("Invalid base point.")
      all(as.vector(as.matrix(X) %*% p) > -1 + 1e-12)
    },
    base_point = list(
      auto = function(X, ...) {
        m <- colMeans(as.matrix(X))
        list(p = m / sqrt(sum(m^2)), convergence = NA)
      }
    )
  )
}

# ---- Base point of the Grassmann experiment -----------------------------

# U0 (4 x 2), the Q factor of a Gaussian matrix drawn after set.seed(2026),
# stored in inst/extdata so that the R and the Python pipelines use the same
# matrix
grassmann_basepoint <- function() {
  as.matrix(read.csv(system.file("extdata", "grassmann_basepoint.csv",
                                 package = "manifoldtest.simulations",
                                 mustWork = TRUE)))
}

# ---- Tests of the paper through manifold.dcov.test ----------------------

# DistCov (center_outward = FALSE) or RankCov (TRUE) p-value, with the base
# points fixed as in the paper. The warning that DistCov raises on hyperbolic
# space is expected here: comparing DistCov with RankCov is the point.
paper_test <- function(X, Y, manifold, px_base, py_base = px_base,
                       center_outward = FALSE, R = 199) {
  res <- suppressWarnings(
    manifold.dcov.test(X, Y, manifold = manifold, px_base = px_base,
                       py_base = py_base, R = R, center_outward = center_outward,
                       seed = GRID_SEED, verbose = FALSE)
  )
  # manifold.dcov.test replaces a base point lying on the cut locus of the
  # data: make sure the base points of the design were kept
  if (!isTRUE(all.equal(as.vector(res$px_base), as.vector(px_base))) ||
      !isTRUE(all.equal(as.vector(res$py_base), as.vector(py_base)))) {
    stop("The prescribed base point was replaced (cut locus of the data).")
  }
  res$dcov_test$p.value
}
