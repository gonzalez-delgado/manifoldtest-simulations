# ================================================================
# DistCov rejection rate on the Section 1.2 torus design, for the base
# points p = (0,0,0), (0.25,0.25,0.25), (0.5,0.5,0.5) and (0.1,0.4,0.7).
#
# The design, the seeds and the number of replications are those of
# figure1_torus.R with N = 300; only the base point of the
# inverse exponential map changes. DistCov is run with
# manifoldtest::manifold.dcov.test.
#
# Output: console table only.
# ================================================================
suppressMessages(library(manifoldtest.simulations))

NP <- 199; NREPS <- 100; N <- 300
N_DIM <- 3; DIM <- N_DIM * N_DIM

gen <- function(n, seed) {
  set.seed(seed)
  C <- as.matrix(expand.grid(rep(list(c(0.002, 0.998)), N_DIM)))
  X <- matrix(0, n, N_DIM)
  for (i in 1:n) X[i, ] <- rnorm(N_DIM, C[sample(1:8, 1), ], 0.05)
  X <- X %% 1
  cd <- pmin(abs(X - 0.5), 1 - abs(X - 0.5))
  S <- as.vector(cd %*% c(1.2, 0.8, 1.5)) +
       0.15 * (cd[, 1] * cd[, 2] + cd[, 1] * cd[, 3] + cd[, 2] * cd[, 3])
  a <- (rank(S) - 0.5) / n
  list(X = X, Y = matrix(((a + 0.5 * rbinom(n, 1, 0.5)) %% 1), n, N_DIM))
}

cat(sprintf("DistCov rejection rate at N = %d, %d replications\n", N, NREPS))
for (p in list(c(0, 0, 0), c(0.25, 0.25, 0.25), c(0.5, 0.5, 0.5), c(0.1, 0.4, 0.7))) {
  s <- 0
  for (r in 1:NREPS) {
    g <- gen(N, 4000 + r)
    if (paper_test(g$X, g$Y, "torus", px_base = p, R = NP) < .05) s <- s + 1
  }
  cat(sprintf("  p = (%.2f, %.2f, %.2f) : %d%%\n", p[1], p[2], p[3], s))
}
