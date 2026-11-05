# ================================================================
# DistCov power as a function of the base point (Section 4.3: 83% -> 75%,
# 71%, 32%, 6%), on the hyperbolic design of Figure 2 (Dep scan at
# noise = 5.0, Dep = 0.044, N = 150).
#
# The manifold points are re-mapped at base points p' at hyperbolic distance
# t = 0.5, 1, 2, 4 from the data-generating base point p0 = (1,0,0,0), in
# three directions; X and Y are always mapped at the SAME p'. Taking p' = p0
# reproduces the published curve and serves as a control. DistCov is run
# with manifoldtest::manifold.dcov.test; for a given replication, all the
# base points share the same permutations.
#
# Data  : hyperbolic_shared_data(), the samples of Figure 2
# Output: basepoint_power.csv
# ================================================================
suppressMessages(library(manifoldtest.simulations))

NP <- 199
P0 <- c(1, 0, 0, 0)

# candidate base points: hyperbolic distance t from p0, three directions
W3 <- diag(3)
cands <- list(list(t = 0, w = W3[1, ]))
for (t in c(0.5, 1, 2, 4)) for (k in 1:3) cands[[length(cands) + 1]] <- list(t = t, w = W3[k, ])

rows <- list()
for (st in list(list(tag = "Dep=0.044, N=150", exp = "distcov", scan = "dep", xval = 5.0))) {
  s0 <- hyperbolic_shared_data(st$exp, scan = st$scan, xval = st$xval)
  reps <- sort(unique(s0$rep))
  rej <- numeric(length(cands))
  for (r in seq_along(reps)) {
    s <- s0[s0$rep == reps[r], ]
    X <- exp_map_hyperbolic(as.matrix(s[, c("u1", "u2", "u3")]), P0)
    Y <- exp_map_hyperbolic(as.matrix(s[, c("v1", "v2", "v3")]), P0)
    for (j in seq_along(cands)) {
      p <- as.vector(exp_map_hyperbolic(matrix(cands[[j]]$t * cands[[j]]$w, 1), P0))
      set.seed(4242 + r)                      # same permutations for every p'
      if (paper_test(X, Y, "hyperbolic", px_base = p, R = NP) < .05) rej[j] <- rej[j] + 1
    }
  }
  n <- length(reps)
  cat(sprintf("\n=== %s   reps=%d ===\n", st$tag, n))
  cat(sprintf("%26s | %9s %9s %9s\n", "distance t of p' from p0", "dir1", "dir2", "dir3"))
  cat(sprintf("%26.1f | %9.1f   (p' = p0, control)\n", 0, 100 * rej[1] / n))
  for (g in seq(2, length(cands), by = 3)) {
    vals <- 100 * rej[g:(g + 2)] / n
    cat(sprintf("%26.1f | %9.1f %9.1f %9.1f\n", cands[[g]]$t, vals[1], vals[2], vals[3]))
    rows[[length(rows) + 1]] <- data.frame(setting = st$tag, t = cands[[g]]$t,
                                           dir1 = vals[1], dir2 = vals[2], dir3 = vals[3])
  }
}
write.csv(do.call(rbind, rows), "basepoint_power.csv", row.names = FALSE)
