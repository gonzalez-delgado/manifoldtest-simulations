# ================================================================
# Antipodal design on S^2 (Section 5.4: 6.7% / 13.3% / 100%).
#
# X uniform on S^2 and, given X, Y uniform on the antipodal pair {X, -X}.
#
# Tests (199 permutations, 60 replications, N = 200 and 400):
#     Geo-dCov : distance covariance of the geodesic distance
#     DistCov  : manifoldtest::manifold.dcov.test on the sphere (built with
#                manifoldtest::new_manifold, see R/helpers.R), base point
#                the north pole
#
# Output: console table only.
# ================================================================
suppressMessages(library(manifoldtest.simulations))

NP <- 199; REPS <- 60; SEED0 <- 4000
NORTH <- c(0, 0, 1)
S2 <- manifold_sphere()

rand_s2 <- function(n) {
  z <- runif(n, -1, 1); ph <- runif(n, 0, 2 * pi); s <- sqrt(1 - z^2)
  cbind(s * cos(ph), s * sin(ph), z)
}

gen_antip <- function(n, seed) {
  set.seed(seed); X <- rand_s2(n); B <- rbinom(n, 1, 0.5)
  list(X = X, Y = X * ifelse(B == 0, 1, -1))
}

power <- function(gen, n) {
  out <- c(Geo_dCov = 0, DistCov = 0)
  for (k in 0:(REPS - 1)) {
    g <- gen(n, SEED0 + k)
    # the two tests share the same permutations
    set.seed(SEED0 + k + 70000)
    out["Geo_dCov"] <- out["Geo_dCov"] +
      (dcov.test(sphere_dist(g$X), sphere_dist(g$Y), R = NP)$p.value < .05)
    set.seed(SEED0 + k + 70000)
    out["DistCov"] <- out["DistCov"] +
      (paper_test(g$X, g$Y, S2, px_base = NORTH, R = NP) < .05)
  }
  100 * out / REPS
}

cat("S^2 antipodal design     N   Geo-dCov   DistCov\n")
for (n in c(200, 400)) {
  r <- power(gen_antip, n)
  cat(sprintf("                      %4d %10.1f %9.1f\n", n, r[1], r[2]))
}
