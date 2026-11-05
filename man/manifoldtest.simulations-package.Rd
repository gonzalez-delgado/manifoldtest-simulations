\name{manifoldtest.simulations-package}
\alias{manifoldtest.simulations-package}
\alias{manifoldtest.simulations}
\alias{GRID_SEED}
\alias{paper_test}
\alias{torus_dist}
\alias{hyperbolic_dist}
\alias{grassmann_dist}
\alias{sphere_dist}
\alias{log_map_sphere}
\alias{exp_map_sphere}
\alias{manifold_sphere}
\alias{grassmann_basepoint}
\alias{hyperbolic_shared_data}
\alias{geo_dcov_curves}
\alias{run_hhg}
\alias{python_grassmann_dist}
\alias{power_plot}
\alias{save_power_plot}
\docType{package}
\title{Simulations for "Independence Test on Riemannian Manifolds"}
\description{
  Code that produces the numerical experiments of Liu, Li, González-Delgado
  and Hallin (2026), "Independence test on Riemannian manifolds". The tests
  proposed in the paper, DistCov and RankCov, are run with
  \code{manifoldtest::manifold.dcov.test}; this package adds the benchmark
  procedures and one R script per figure and per number quoted in the text.

  The scripts are installed in
  \code{system.file("scripts", package = "manifoldtest.simulations")} and
  write their results to the working directory. See the README of the
  repository for the map between the scripts and the paper.

  The exported functions are helpers shared by the scripts:
  \describe{
    \item{\code{paper_test}, \code{GRID_SEED}}{DistCov or RankCov p-value
      at the base points of a design.}
    \item{\code{torus_dist}, \code{hyperbolic_dist}, \code{grassmann_dist},
      \code{sphere_dist}}{pairwise geodesic distance matrices for the
      benchmark tests.}
    \item{\code{manifold_sphere}, \code{log_map_sphere},
      \code{exp_map_sphere}}{the unit sphere as a \pkg{manifoldtest}
      manifold.}
    \item{\code{grassmann_basepoint}}{the fixed base point of the
      Grassmann experiment.}
    \item{\code{hyperbolic_shared_data}, \code{geo_dcov_curves}}{the
      samples of the hyperbolic experiments shared by all procedures, and
      the Geo-dCov test on them.}
    \item{\code{run_hhg}, \code{python_grassmann_dist}}{the HHG benchmark,
      run in Python (\code{inst/python/hhg.py}).}
    \item{\code{power_plot}, \code{save_power_plot}}{power curves of
      Figures 2 and 3.}
  }
}
\keyword{internal}
