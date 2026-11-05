# ================================================================
# Figure 1 (Section 1.2): rejection rates on the torus T^3.
#
# X-margin: eight cluster centres at the vertices of the cube
# {0.002, 0.998}^3; an observation chooses a centre uniformly, is perturbed
# by Gaussian noise of standard deviation 0.05, and has every coordinate
# reduced modulo 1.
#
# Signal and response, with cd(x)_j = min(|x_j - 0.5|, 1 - |x_j - 0.5|):
#     S(X) = sum_j cd(X)_j w_j + 0.15 sum_{i<i'} cd(X)_i cd(X)_i'
#     a(X) = (rank{S(X)} - 1/2) / N
#     Y    = ((a + 0.5 B) mod 1) * 1_3,    B ~ Bernoulli(1/2)
# so that Y lies on the diagonal circle {(u,u,u)} of T^3, a closed geodesic
# of arc length sqrt(3).
#
# Tests: DistCov (manifoldtest::manifold.dcov.test on the torus, base point
# 0), Raw dCov (chart [0,1)^3), Geo dCov (torus distance), 199 permutations,
# 100 replications.
#
# Scans:  N = 100, 200, 300, 500, 800   (dependence held fixed)
#         null                          (Y independent of X, N = 300)
#
# Output (in the working directory):
#   Fig_Torus_Motivation.png
#   torus_power.csv, torus_size.csv   the numbers behind the figure
# ================================================================
suppressMessages({
  library(manifoldtest.simulations)
  library(ggplot2); library(reshape2)
})

N_PERM <- 199
N_REPS <- 100
NC <- 8                 # corner clusters
SD_NOISE <- 0.05        # spread of a cluster, on the torus
P_CORNER <- 0.002       # the corner coordinates are 0.002 and 0.998
W <- c(1.2, 0.8, 1.5)   # weights of the three coordinates in the signal

# ---- the design --------------------------------------------------------
gen <- function(n, d = 3, seed, null = FALSE) {
  set.seed(seed)
  C <- as.matrix(expand.grid(rep(list(c(P_CORNER, 1 - P_CORNER)), d)))
  X <- matrix(0, n, d)
  for (i in 1:n) X[i, ] <- rnorm(d, C[sample(1:NC, 1), ], SD_NOISE)
  X <- X %% 1
  cd <- pmin(abs(X - 0.5), 1 - abs(X - 0.5))          # distance to the centre
  S <- as.vector(cd %*% W) + 0.15 * (cd[, 1] * cd[, 2] + cd[, 1] * cd[, 3] +
                                     cd[, 2] * cd[, 3])
  a <- (rank(S) - 0.5) / n
  B <- rbinom(n, 1, 0.5)
  if (null) a <- runif(n)                             # Y independent of X
  list(X = X, Y = matrix(((a + 0.5 * B) %% 1), n, d))
}

# ---- the three procedures ----------------------------------------------
run_rep <- function(n, seed, null = FALSE) {
  g <- gen(n, seed = seed, null = null); p0 <- c(0, 0, 0)
  c(DistCov  = paper_test(g$X, g$Y, "torus", px_base = p0, R = N_PERM),
    Raw_dCov = dcov.test(dist(g$X), dist(g$Y), R = N_PERM)$p.value,
    Geo_dCov = dcov.test(torus_dist(g$X), torus_dist(g$Y), R = N_PERM)$p.value)
}

cat(sprintf("%5s %9s %9s %9s\n", "N", "DistCov", "Raw dCov", "Geo dCov"))
out <- NULL; t0 <- Sys.time()
for (n in c(100, 200, 300, 500, 800)) {
  P <- t(sapply(1:N_REPS, function(r) run_rep(n, 4000 + r)))
  pw <- round(100 * colMeans(P < 0.05), 1)
  cat(sprintf("%5d %9.1f %9.1f %9.1f\n", n, pw[1], pw[2], pw[3])); flush.console()
  out <- rbind(out, data.frame(N = n, DistCov = pw[1], Raw_dCov = pw[2], Geo_dCov = pw[3]))
}
write.csv(out, "torus_power.csv", row.names = FALSE)

P0 <- t(sapply(1:N_REPS, function(r) run_rep(300, 4000 + r, null = TRUE)))
sz <- round(100 * colMeans(P0 < 0.05), 1)
cat(sprintf("\nnull (Y independent of X), N = 300 : DistCov %.1f  Raw %.1f  Geo %.1f\n",
            sz[1], sz[2], sz[3]))
write.csv(data.frame(N = 300, DistCov = sz[1], Raw_dCov = sz[2], Geo_dCov = sz[3]),
          "torus_size.csv", row.names = FALSE)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- figure ------------------------------------------------------------

res <- out
df <- melt(res, id.vars = "N", variable.name = "Method", value.name = "Power")
df$Method <- factor(df$Method, levels = c("DistCov", "Raw_dCov", "Geo_dCov"),
                    labels = c("DistCov", "Raw dCov", "Geo dCov"))
PAL <- c(DistCov = "#E41A1C", "Raw dCov" = "#984EA3", "Geo dCov" = "#00BFC4")
LTY <- c(DistCov = "solid",   "Raw dCov" = "dotdash", "Geo dCov" = "twodash")
SHP <- c(DistCov = 16,        "Raw dCov" = 3,         "Geo dCov" = 4)

p <- ggplot(df, aes(x = N, y = Power, color = Method, linetype = Method,
                    shape = Method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_color_manual(values = PAL) + scale_linetype_manual(values = LTY) +
  scale_shape_manual(values = SHP) +
  scale_x_continuous(breaks = res$N) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  geom_hline(yintercept = 5, linetype = "dotted", linewidth = 0.4, alpha = 0.6) +
  annotate("text", x = min(res$N), y = 1.2, label = expression(alpha == 0.05),
           hjust = 0, size = 3.0, color = "gray40") +
  labs(title = "Rejection rates on the flat torus",
       x = "Sample size N", y = "Power (%)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.width = unit(1.0, "cm"), legend.key.height = unit(0.65, "cm"),
        legend.box = "horizontal", legend.margin = margin(t = 4, b = 0),
        legend.box.margin = margin(t = -2),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        panel.grid.minor = element_blank(),
        plot.margin = margin(6, 8, 2, 6))
ggsave("Fig_Torus_Motivation.png", p, width = 7.2, height = 5.4, dpi = 150)
cat("saved Fig_Torus_Motivation.png\n")
print(res, row.names = FALSE)
