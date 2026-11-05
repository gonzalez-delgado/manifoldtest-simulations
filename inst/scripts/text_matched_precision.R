# ================================================================
# Matched-precision experiment (Section 5.4: 100 / 100 / 79 on R^2, H^2 and
# T^2 at N = 100, and 100 on all three at N = 300).
#
# One and the same tangent-space problem is realised on three spaces -- the
# Euclidean plane R^2, hyperbolic space H^2 and the flat torus T^2 -- and
# each is tested with DistCov (manifoldtest::manifold.dcov.test: inverse
# exponential map, then distance covariance). The first two spaces are of
# strong negative type, the third is not of negative type at all. DistCov
# reaches full power on all three at N = 300, while the deficit that appears
# on the torus at smaller N is a cut-locus effect that vanishes with N.
#
# Output: console table only.
# ================================================================
suppressMessages(library(manifoldtest.simulations))

N_VALUES <- c(100, 300); REP <- 100; NP <- 199
W <- matrix(c(1.0, -0.3, 0.4, 1.2), 2, 2)   # rows (1.0, 0.4) and (-0.3, 1.2)
P_H2 <- c(1, 0, 0)                          # base point on H^2
P_T2 <- c(0.5, 0.5)                         # base point on T^2

for (N in N_VALUES) {
  res <- c(R2 = 0, H2 = 0, T2 = 0)
  for (b in seq_len(REP)) {
    set.seed(88000 + b)
    U <- matrix(runif(2 * N, -0.45, 0.45), N, 2)           # tangent-space data
    V <- (U / sqrt(rowSums(U^2))) %*% t(W) + matrix(rnorm(2 * N, sd = 0.1), N, 2)

    # (a) Euclidean R^2: the logarithmic map is the identity, so DistCov is
    #     the distance covariance test of the data themselves
    res["R2"] <- res["R2"] + (dcov.test(U, V, R = NP)$p.value < .05)

    # (b) hyperbolic H^2 (strong negative type)
    X <- exp_map_hyperbolic(U, P_H2); Y <- exp_map_hyperbolic(V, P_H2)
    res["H2"] <- res["H2"] + (paper_test(X, Y, "hyperbolic", px_base = P_H2, R = NP) < .05)

    # (c) flat torus T^2 (not even of negative type). The data are wrapped
    #     around 0 and mapped at the base point (0.5, 0.5), so that they
    #     straddle its cut locus: this is the cut-locus effect of the text
    Xt <- U %% 1; Yt <- V %% 1
    res["T2"] <- res["T2"] + (paper_test(Xt, Yt, "torus", px_base = P_T2, R = NP) < .05)
  }

  cat(sprintf("\nDistCov, same tangent-space dependence on three spaces (N=%d, %d replications, nominal 5%%)\n",
              N, REP))
  cat(sprintf("  Euclidean R^2  : %5.1f%%\n", 100 * res["R2"] / REP))
  cat(sprintf("  Hyperbolic H^2 : %5.1f%%\n", 100 * res["H2"] / REP))
  cat(sprintf("  Flat torus T^2 : %5.1f%%\n", 100 * res["T2"] / REP))
}
