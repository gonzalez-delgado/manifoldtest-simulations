# ================================================================
# Distance covariance of the Grassmann design of figure4_grassmann.R on a
# large sample (Section 5.3: 6.3e-4 / 1.3e-2 / 1.0e-3, N = 1500).
#
# On the response curve gamma the intrinsic distance is the circle distance,
# so E[d_g(Y,Y') | X, X'] = L/4 for every pair X, X': the distance
# covariance of the design is essentially zero for the geodesic metric and
# positive for the tangent chart. The same statistic against an independent
# copy of the response gives the finite-sample bias of the estimator.
#
# Output: console only.
# ================================================================
suppressMessages(library(manifoldtest.simulations))

N_DIM <- 4; K_DIM <- 2; DIM <- N_DIM * K_DIM
U0 <- grassmann_basepoint()
horiz <- function(M) M - U0 %*% (t(U0) %*% M)
unit_norm <- function(M) { s <- sqrt(sum(M^2)); if (s < 1e-14) matrix(0, nrow(M), ncol(M)) else M / s }

# The maps of manifoldtest, written for horizontal lifts Delta (n x k,
# U0' Delta = 0): manifoldtest works with the coordinates U0_perp' Delta
U0_perp <- qr.Q(qr(U0), complete = TRUE)[, -seq_len(K_DIM), drop = FALSE]
exp_grassmann <- function(U0, Delta)
  matrix(exp_map_grassmann(matrix(crossprod(U0_perp, Delta), 1), U0), N_DIM, K_DIM)
log_grassmann <- function(V)
  U0_perp %*% matrix(log_map_grassmann(matrix(V, 1), U0, N_DIM, K_DIM), N_DIM - K_DIM, K_DIM)
# pairwise geodesic distance matrix of the package, for (N, n, k) arrays
gr_dist <- function(U) grassmann_dist(matrix(U, dim(U)[1]), N_DIM, K_DIM)

# ---- the response curve, as in figure4_grassmann.R ----------------------
set.seed(999)
Q <- qr.Q(qr(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))))
set.seed(4242); W_SIG <- rnorm(DIM)   # fixed direction defining the signal
gamma <- function(theta) { if (length(theta) == 1) theta <- rep(theta, K_DIM)
  U0 %*% diag(cos(theta), K_DIM, K_DIM) + Q %*% diag(sin(theta), K_DIM, K_DIM) }

# ---- distance covariance on N = 1500 observations -----------------------
V2 <- function(D1, D2) { n <- nrow(D1)
  A <- D1 - rowMeans(D1) - rep(colMeans(D1), each = n) + mean(D1)
  B <- D2 - rowMeans(D2) - rep(colMeans(D2), each = n) + mean(D2)
  sum(A * B) / n^2 }
N <- 1500; set.seed(777)
Xt <- matrix(0, N, DIM)
for (i in 1:N) Xt[i,] <- as.vector(runif(1, 0.2, 1.4) * unit_norm(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))))
a <- pi * (rank(as.vector(Xt %*% W_SIG)) - 0.5) / N
B <- rbinom(N, 1, 0.5); theta <- (a + (pi/2)*B) %% pi
Yt <- t(sapply(theta, function(t) as.vector(log_grassmann(gamma(t)))))
Xs <- array(0, c(N, N_DIM, K_DIM)); Ys <- array(0, c(N, N_DIM, K_DIM))
for (i in 1:N) { Xs[i,,] <- exp_grassmann(U0, matrix(Xt[i,], N_DIM, K_DIM)); Ys[i,,] <- gamma(theta[i]) }
cat(sprintf("distance covariance, geodesic metric   : %.2e\n", V2(gr_dist(Xs), gr_dist(Ys))))
cat(sprintf("distance covariance, tangent chart     : %.2e\n", V2(as.matrix(dist(Xt)), as.matrix(dist(Yt)))))
# the same statistic against an independent copy of Y, which is the
# finite-sample bias of V^2_N and the level at which the estimate sits
Yb <- Ys[sample(N), , , drop = FALSE]
cat(sprintf("against an independent copy of Y       : %.2e\n",
            V2(as.matrix(dist(Xt)), gr_dist(Yb))))
