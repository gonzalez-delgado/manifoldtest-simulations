# ================================================================
# Rejection rates on S^2 under three designs and three tests
# (Section 5.4: 6.7 / 13.3 / 100).
#
# Designs:
#     antipodal : Y | X uniform on {X, -X}
#     rotate    : Y = R_z(pi * X_3 + N(0, 0.3)) X
#     chordflip : Y uniform on S^2, with the sign of its first coordinate
#                 tied to that of X
#
# Tests (199 permutations, 60 replications):
#     intrinsic : distance covariance of the geodesic distance
#     chord     : distance covariance of the chord (ambient R^3) distance
#     DistCov   : manifoldtest::manifold.dcov.test on the sphere (built with
#                 manifoldtest::new_manifold, see R/helpers.R), base
#                 point the north pole
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
rot_z <- function(X, a) {
  cbind(cos(a) * X[, 1] - sin(a) * X[, 2], sin(a) * X[, 1] + cos(a) * X[, 2], X[, 3])
}

gen_antip <- function(n, seed) {
  set.seed(seed); X <- rand_s2(n); B <- rbinom(n, 1, 0.5)
  list(X = X, Y = X * ifelse(B == 0, 1, -1))
}
gen_rotate <- function(n, seed) {          # rotation angle driven by latitude + noise
  set.seed(seed); X <- rand_s2(n)
  list(X = X, Y = rot_z(X, pi * X[, 3] + rnorm(n, 0, 0.3)))
}
gen_chord <- function(n, seed) {           # dependence simple in the ambient R^3 coords
  set.seed(seed); X <- rand_s2(n); Y <- rand_s2(n)
  Y[, 1] <- sign(X[, 1]) * abs(Y[, 1])
  list(X = X, Y = Y / sqrt(rowSums(Y^2)))
}

power <- function(gen, n) {
  out <- c(intrinsic = 0, chord = 0, DistCov = 0)
  for (k in 0:(REPS - 1)) {
    g <- gen(n, SEED0 + k)
    # the three tests share the same permutations
    set.seed(SEED0 + k + 70000)
    out["intrinsic"] <- out["intrinsic"] +
      (dcov.test(sphere_dist(g$X), sphere_dist(g$Y), R = NP)$p.value < .05)
    set.seed(SEED0 + k + 70000)
    out["chord"] <- out["chord"] + (dcov.test(g$X, g$Y, R = NP)$p.value < .05)
    set.seed(SEED0 + k + 70000)
    out["DistCov"] <- out["DistCov"] +
      (paper_test(g$X, g$Y, S2, px_base = NORTH, R = NP) < .05)
  }
  100 * out / REPS
}

cat("S^2  design      n   intrinsic  raw(chord)  DistCov\n")
for (d in list(list("antipodal", gen_antip), list("rotate", gen_rotate),
               list("chordflip", gen_chord))) {
  for (n in c(100, 200, 400)) {
    r <- power(d[[2]], n)
    cat(sprintf("     %-10s %4d %10.1f %11.1f %8.1f\n", d[[1]], n, r[1], r[2], r[3]))
  }
}
