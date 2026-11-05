# ================================================================
# Matched-precision experiment (quoted in Section 4.3).
#
# One and the same tangent-space problem is realised on three spaces -- the
# Euclidean plane R^2, hyperbolic space H^2 and the flat torus T^2 -- and each
# is tested twice:
#
#     DistCov         : manifoldtest::manifold.dcov.test (inverse exponential
#                       map, then distance covariance)
#     dCov(intrinsic) : distance covariance fed the space's own metric
#
# The first two spaces are of strong negative type, the third is not of
# negative type at all. If the mapping were losing information, DistCov would
# degrade on the torus; it does not: it reaches full power on all three at
# N = 300, while the deficit that appears on the torus at smaller N is a
# cut-locus effect that vanishes with N, not a loss of validity.
#
# Output: console table only.
# ================================================================
suppressMessages(library(manifoldtest.simulations))

N_VALUES <- c(100, 300); REP <- 100; NP <- 199
W <- matrix(c(1.0, -0.3, 0.4, 1.2), 2, 2)   # rows (1.0, 0.4) and (-0.3, 1.2)
P_H2 <- c(1, 0, 0)                          # base point on H^2
P_T2 <- c(0.5, 0.5)                         # base point on T^2

for (N in N_VALUES) {
  res <- c(R2_Dis = 0, R2_int = 0, H2_Dis = 0, H2_int = 0, T2_Dis = 0, T2_int = 0)
  for (b in seq_len(REP)) {
    set.seed(88000 + b)
    U <- matrix(runif(2 * N, -0.45, 0.45), N, 2)           # tangent-space data
    V <- (U / sqrt(rowSums(U^2))) %*% t(W) + matrix(rnorm(2 * N, sd = 0.1), N, 2)

    # (a) Euclidean R^2: the intrinsic metric IS the Euclidean one, so the two
    #     tests are the same computation and give the same p-value
    r2 <- dcov.test(U, V, R = NP)$p.value < .05
    res["R2_Dis"] <- res["R2_Dis"] + r2; res["R2_int"] <- res["R2_int"] + r2

    # (b) hyperbolic H^2 (strong negative type)
    X <- exp_map_hyperbolic(U, P_H2); Y <- exp_map_hyperbolic(V, P_H2)
    res["H2_Dis"] <- res["H2_Dis"] + (paper_test(X, Y, "hyperbolic", px_base = P_H2, R = NP) < .05)
    res["H2_int"] <- res["H2_int"] +
      (dcov.test(hyperbolic_dist(X), hyperbolic_dist(Y), R = NP)$p.value < .05)

    # (c) flat torus T^2 (not even of negative type). The data are wrapped
    #     around 0 and mapped at the base point (0.5, 0.5), so that they
    #     straddle its cut locus: this is the cut-locus effect of the text
    Xt <- U %% 1; Yt <- V %% 1
    res["T2_Dis"] <- res["T2_Dis"] + (paper_test(Xt, Yt, "torus", px_base = P_T2, R = NP) < .05)
    res["T2_int"] <- res["T2_int"] +
      (dcov.test(torus_dist(Xt), torus_dist(Yt), R = NP)$p.value < .05)
  }

  cat(sprintf("\nSame tangent-space dependence, realised on three spaces (N=%d, %d replications, nominal 5%%)\n",
              N, REP))
  cat(sprintf("%16s %22s %9s %16s\n", "space", "type", "DistCov", "dCov(intrinsic)"))
  for (row in list(c("Euclidean R^2", "strong negative type", "R2"),
                   c("Hyperbolic H^2", "strong negative type", "H2"),
                   c("Flat torus T^2", "not of negative type", "T2"))) {
    cat(sprintf("%16s %22s %8.1f%% %15.1f%%\n", row[1], row[2],
                100 * res[paste0(row[3], "_Dis")] / REP,
                100 * res[paste0(row[3], "_int")] / REP))
  }
}
