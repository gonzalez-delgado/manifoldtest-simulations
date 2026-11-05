# ================================================================
# Self-checks of the designs; they are not figures of the paper.
#
# 1. Torus design of figure1_torus.R:
#    - the diagonal of T^3 is a closed geodesic of arc length sqrt(3) and
#      the torus distance restricted to it is the circle distance;
#    - the antipodal identity holds on it, so the conditional distance
#      profile of the response is constant;
#    - the population distance covariance of the geodesic metric is
#      annihilated, while the chart and the raw coordinates keep theirs.
#
# 2. Grassmann design of figure4_grassmann.R:
#    - the exponential and logarithmic maps of Gr(n,k) are inverse to each
#      other and the exponential map is a radial isometry;
#    - the fast geodesic distance agrees with a direct loop over the
#      singular values of U^T V;
#    - the response curve gamma is a closed geodesic of arc length
#      pi*sqrt(k), lies in the cut domain of U0, and carries the circle
#      distance and the antipodal identity the design rests on;
#    - the population distance covariance of the design is essentially zero
#      for the geodesic metric and positive for the tangent chart.
#
# 3. Negative type of the metrics used in the paper. (X,d) is of negative
#    type iff for every finite set x_1..x_N and every vector a with
#    sum(a) = 0 one has a'Da <= 0, where D = (d(x_i,x_j)); equivalently the
#    doubly centred matrix B = -P D P, P = I - 11'/N, is positive
#    semidefinite. A single negative eigenvalue of B is therefore a witness
#    that d is not of negative type. Distance covariance needs the strictly
#    stronger condition of strong negative type (Lyons 2013), so a space
#    that fails here is outside the consistency guarantee of the
#    intrinsic-metric test of Section 1.2.
#
# Each block runs in its own environment and sets its own seeds.
#
# Output: console, and negative_type.csv (space, N, min_eig, min_eig/max_eig)
# ================================================================
suppressMessages(library(manifoldtest.simulations))

# ================================================================
# 1. Torus design
# ================================================================
cat("==== 1. Torus design ====\n")
local({
  W <- c(1.2, 0.8, 1.5); NC <- 8; P_CORNER <- 0.002; SD_NOISE <- 0.05
  gen <- function(n, d = 3, seed, null = FALSE) {
    set.seed(seed)
    C <- as.matrix(expand.grid(rep(list(c(P_CORNER, 1 - P_CORNER)), d)))
    X <- matrix(0, n, d)
    for (i in 1:n) X[i,] <- rnorm(d, C[sample(1:NC,1),], SD_NOISE)
    X <- X %% 1
    cd <- pmin(abs(X-0.5), 1-abs(X-0.5))
    S <- as.vector(cd %*% W) + 0.15*(cd[,1]*cd[,2]+cd[,1]*cd[,3]+cd[,2]*cd[,3])
    a <- (rank(S) - 0.5)/n
    B <- rbinom(n, 1, 0.5)
    if (null) a <- runif(n)
    list(X = X, Y = matrix(((a + 0.5*B) %% 1), n, d), theta = (a + 0.5*B) %% 1)
  }
  V2 <- function(D1, D2) { n <- nrow(D1)
    A <- D1 - rowMeans(D1) - rep(colMeans(D1), each = n) + mean(D1)
    B <- D2 - rowMeans(D2) - rep(colMeans(D2), each = n) + mean(D2)
    sum(A*B)/n^2 }

  # ---- 1-2. the diagonal circle and the antipodal identity ---------------
  u <- seq(1e-4, 1 - 1e-4, length.out = 200)
  Dg <- as.matrix(torus_dist(cbind(u, u, u)))
  dth <- abs(outer(u, u, "-"))
  cat(sprintf("diameter of the diagonal circle      : %.6f   (sqrt(3)/2 = %.6f)\n",
              max(Dg), sqrt(3)/2))
  cat(sprintf("d_T on the diagonal = sqrt(3)*min(|du|,1-|du|) : %.2e\n",
              max(abs(Dg - sqrt(3) * pmin(dth, 1 - dth)))))
  set.seed(3); i <- sample(200, 60); j <- sample(200, 60)
  sum_id <- sapply(seq_along(i), function(k)
    Dg[i[k], j[k]] + sqrt(3) * min(abs(u[i[k]] - ((u[j[k]] + 0.5) %% 1)),
                                   1 - abs(u[i[k]] - ((u[j[k]] + 0.5) %% 1))))
  cat(sprintf("antipodal identity, sd of the sum    : %.2e   (sqrt(3)/2 = %.6f)\n",
              sd(sum_id), sqrt(3)/2))

  # ---- 3. population distance covariance ---------------------------------
  N <- 3000; g <- gen(N, seed = 4001)
  gn <- gen(N, seed = 4001, null = TRUE)
  TX <- as.matrix(log_map_torus(g$X, c(0,0,0))); TY <- as.matrix(log_map_torus(g$Y, c(0,0,0)))
  cat(sprintf("\npopulation V^2, tangent chart (DistCov) : %.6f\n",
              V2(as.matrix(dist(TX)), as.matrix(dist(TY)))))
  cat(sprintf("population V^2, raw chart     (Raw dCov): %.6f\n",
              V2(as.matrix(dist(g$X)), as.matrix(dist(g$Y)))))
  cat(sprintf("population V^2, geodesic      (Geo dCov): %.8f   <- the obstruction\n",
              V2(as.matrix(torus_dist(g$X)), as.matrix(torus_dist(g$Y)))))
  cat("\n            the same three statistics against an independent response\n")
  cat("            (each is the finite-sample bias of the estimator, at N = 3000):\n")
  cat(sprintf("              chart %.6f   raw %.6f   geodesic %.8f\n",
              V2(as.matrix(dist(TX)), as.matrix(dist(log_map_torus(gn$Y, c(0,0,0))))),
              V2(as.matrix(dist(g$X)), as.matrix(dist(gn$Y))),
              V2(as.matrix(torus_dist(g$X)), as.matrix(torus_dist(gn$Y)))))
  cat(sprintf("mean d_T(Y,Y') : %.6f   (sqrt(3)/4 = %.6f)\n",
              mean(torus_dist(g$Y)), sqrt(3)/4))
})

# ================================================================
# 2. Grassmann design
# ================================================================
cat("\n==== 2. Grassmann design ====\n")
local({
  N_DIM <- 4; K_DIM <- 2; DIM <- N_DIM * K_DIM
  U0 <- grassmann_basepoint()
  horiz <- function(M) M - U0 %*% (t(U0) %*% M)
  unit  <- function(M) { s <- sqrt(sum(M^2)); if (s < 1e-14) matrix(0, nrow(M), ncol(M)) else M / s }

  # The maps of manifoldtest, written for horizontal lifts Delta (n x k,
  # U0' Delta = 0): manifoldtest works with the coordinates U0_perp' Delta
  U0_perp <- qr.Q(qr(U0), complete = TRUE)[, -seq_len(K_DIM), drop = FALSE]
  exp_grassmann <- function(U0, Delta)
    matrix(exp_map_grassmann(matrix(crossprod(U0_perp, Delta), 1), U0), N_DIM, K_DIM)
  log_grassmann <- function(V)
    U0_perp %*% matrix(log_map_grassmann(matrix(V, 1), U0, N_DIM, K_DIM), N_DIM - K_DIM, K_DIM)
  # geodesic distance between two points, from the principal angles
  gdist2 <- function(U, V) {
    s <- pmin(pmax(svd(t(U) %*% V)$d, -1), 1); sqrt(sum(acos(s)^2))
  }
  # fast pairwise distance matrix of the package (R/helpers.R), for (N, n, k) arrays
  gr_dist <- function(U, V) grassmann_dist(matrix(U, dim(U)[1]), N_DIM, K_DIM)

  # ---- 1. the two maps ---------------------------------------------------
  set.seed(1); e1 <- 0; e2 <- 0; e3 <- 0
  for (t in 1:300) {
    D <- unit(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))) * runif(1, 0.05, 1.5)
    X <- exp_grassmann(U0, D)
    e1 <- max(e1, max(abs(t(X) %*% X - diag(K_DIM))))
    e2 <- max(e2, max(abs(log_grassmann(X) - D)))
    e3 <- max(e3, abs(gdist2(U0, X) - sqrt(sum(D^2))))
  }
  cat(sprintf("Stiefel defect                : %.2e\n", e1))
  cat(sprintf("log(exp(Delta)) - Delta       : %.2e\n", e2))
  cat(sprintf("d_g(U0, exp(Delta)) - |Delta| : %.2e\n", e3))

  # ---- 2. fast geodesic distance vs a direct loop ------------------------
  gr_dist_loop <- function(U, V) {
    N <- dim(U)[1]; D <- matrix(0, N, N)
    for (i in 1:N) for (j in 1:N) if (i != j) D[i, j] <- gdist2(U[i,,], V[j,,])
    D
  }
  set.seed(5); M <- 8
  Dl <- lapply(1:M, function(i) unit(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))) * runif(1, 0.1, 1.5))
  U <- array(0, c(M, N_DIM, K_DIM)); for (i in 1:M) U[i,,] <- exp_grassmann(U0, Dl[[i]])
  cat(sprintf("gr_dist vs svd loop           : %.2e\n", max(abs(gr_dist(U, U) - gr_dist_loop(U, U)))))

  # ---- 3. the response curve ---------------------------------------------
  set.seed(999)
  Q <- qr.Q(qr(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))))
  set.seed(4242); W_SIG <- rnorm(DIM)   # fixed direction defining the signal
  gamma <- function(theta) { if (length(theta) == 1) theta <- rep(theta, K_DIM)
    U0 %*% diag(cos(theta), K_DIM, K_DIM) + Q %*% diag(sin(theta), K_DIM, K_DIM) }
  L <- pi * sqrt(K_DIM)
  th <- seq(1e-4, pi - 1e-4, length.out = 200)
  rad <- sapply(th, function(t) sqrt(sum(log_grassmann(gamma(t))^2)))
  cat(sprintf("\nchart radius along gamma, max : %.4f   (cut domain reached at L/2 = %.4f)\n",
              max(rad), L/2))
  cat(sprintf("  the ball {|Delta| < pi/2} covers only radii up to %.4f\n", pi/2))
  set.seed(11); i <- sample(200, 40); j <- sample(200, 40)
  dcur <- sapply(seq_along(i), function(k) gdist2(gamma(th[i[k]]), gamma(th[j[k]])))
  dth <- abs(th[i] - th[j])
  cat(sprintf("circle-distance identity      : %.2e\n",
              max(abs(dcur - sqrt(K_DIM) * pmin(dth, pi - dth)))))
  zz <- sapply(seq_along(i), function(k)
    gdist2(gamma(th[i[k]]), gamma(th[j[k]])) +
    gdist2(gamma(th[i[k]]), gamma((th[j[k]] + pi/2) %% pi)))
  cat(sprintf("antipodal identity, sd of sum : %.2e   (should be 0, L/2 = %.4f)\n", sd(zz), L/2))

  # ---- 4. population distance covariance ---------------------------------
  V2 <- function(D1, D2) { n <- nrow(D1)
    A <- D1 - rowMeans(D1) - rep(colMeans(D1), each = n) + mean(D1)
    B <- D2 - rowMeans(D2) - rep(colMeans(D2), each = n) + mean(D2)
    sum(A * B) / n^2 }
  N <- 1500; set.seed(777)
  Xt <- matrix(0, N, DIM)
  for (i in 1:N) Xt[i,] <- as.vector(runif(1, 0.2, 1.4) * unit(horiz(matrix(rnorm(DIM), N_DIM, K_DIM))))
  a <- pi * (rank(as.vector(Xt %*% W_SIG)) - 0.5) / N
  B <- rbinom(N, 1, 0.5); theta <- (a + (pi/2)*B) %% pi
  Yt <- t(sapply(theta, function(t) as.vector(log_grassmann(gamma(t)))))
  Xs <- array(0, c(N, N_DIM, K_DIM)); Ys <- array(0, c(N, N_DIM, K_DIM))
  for (i in 1:N) { Xs[i,,] <- exp_grassmann(U0, matrix(Xt[i,], N_DIM, K_DIM)); Ys[i,,] <- gamma(theta[i]) }
  cat(sprintf("\npopulation V^2, tangent chart : %.6f\n", V2(as.matrix(dist(Xt)), as.matrix(dist(Yt)))))
  cat(sprintf("population V^2, geodesic      : %.6f   <- the obstruction\n", V2(gr_dist(Xs,Xs), gr_dist(Ys,Ys))))
  cat(sprintf("mean d_g(Y,Y')                : %.6f   (L/4 = %.6f)\n",
              mean(gr_dist(Ys, Ys)), L/4))
  # the same statistic against an independent copy of Y, which is the
  # finite-sample bias of V^2_N and the level at which the estimate sits
  Yb <- Ys[sample(N), , , drop = FALSE]
  cat(sprintf("V^2 against an independent Y  : %.6f   <- finite-sample bias\n",
              V2(as.matrix(dist(Xt)), gr_dist(Yb, Yb))))
})

# ================================================================
# 3. Negative type
# ================================================================
cat("\n==== 3. Negative type ====\n")
local({
  negtype_eigs <- function(D) {
    n <- nrow(D); P <- diag(n) - matrix(1, n, n) / n
    B <- -P %*% D %*% P; B <- (B + t(B)) / 2
    ev <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
    c(min_eig = min(ev), ratio = min(ev) / max(ev))
  }

  # geodesic distance on Gr(n,k), from arbitrary representatives
  gr_dist_loop <- function(U, V) {
    N <- dim(U)[1]; D <- matrix(0, N, N)
    for (i in 1:N) for (j in 1:N) if (i != j) {
      s <- svd(t(U[i, , ]) %*% V[j, , ])$d
      D[i, j] <- sqrt(sum(acos(pmin(pmax(s, -1), 1))^2))
    }
    D
  }
  rand_stiefel <- function(n, k, N) {
    out <- array(0, c(N, n, k))
    for (i in 1:N) out[i, , ] <- qr.Q(qr(matrix(rnorm(n * k), n, k)))
    out
  }
  torus_dist <- function(X) {                      # T^n, as in Section 1.2
    n <- nrow(X); D <- matrix(0, n, n)
    for (i in 1:n) for (j in 1:n) {
      dd <- abs(X[i, ] - X[j, ])
      D[i, j] <- sqrt(sum(pmin(dd, 1 - dd)^2))
    }
    D
  }

  rows <- list()
  add <- function(space, N, D) {
    e <- negtype_eigs(D)
    rows[[length(rows) + 1]] <<- data.frame(space = space, N = N,
                                            min_eig = e["min_eig"], ratio = e["ratio"])
  }

  N <- 400
  set.seed(3)
  E <- matrix(rnorm(N * 21), N, 21)                       # the tangent chart (Euclidean)
  add("R^21 (tangent chart)", N, as.matrix(dist(E)))

  S <- matrix(rnorm(N * 3), N, 3); S <- S / sqrt(rowSums(S^2))
  add("S^2, geodesic", N, acos(pmin(pmax(S %*% t(S), -1), 1)))

  r <- runif(N, 0, 3); th <- runif(N, 0, 2 * pi)
  H <- cbind(cosh(r), sinh(r) * cos(th), sinh(r) * sin(th))
  add("H^2, geodesic", N, acosh(pmax(-(H %*% diag(c(-1, 1, 1)) %*% t(H)), 1 + 1e-12)))

  add("T^2, geodesic", N, torus_dist(matrix(runif(N * 2), N, 2)))
  add("T^3, geodesic", N, torus_dist(matrix(runif(N * 3), N, 3)))

  for (nk in list(c(4, 2), c(10, 3))) for (N in c(150, 400)) {
    n <- nk[1]; k <- nk[2]
    set.seed(7); G <- rand_stiefel(n, k, N)
    add(sprintf("Gr(%d,%d), geodesic", n, k), N, gr_dist_loop(G, G))
    Pr <- t(apply(G, 1, function(A) as.vector(A %*% t(A))))
    add(sprintf("Gr(%d,%d), projection distance", n, k), N, as.matrix(dist(Pr)))
  }

  out <- do.call(rbind, rows)
  write.csv(out, "negative_type.csv", row.names = FALSE)
  cat(sprintf("%-32s %6s %14s %14s\n", "space", "N", "min eigenvalue", "min/max"))
  for (i in seq_len(nrow(out)))
    cat(sprintf("%-32s %6d %14.2e %14.2e   %s\n", out$space[i], out$N[i],
                out$min_eig[i], out$ratio[i],
                ifelse(out$ratio[i] > -1e-8, "negative type", "NOT of negative type")))
})
