# ================================================================
# HHG benchmark (Heller, Heller & Gorfine, 2013), run in Python.
#
# The statistic and its permutation test are in inst/python/hhg.py; the
# figure scripts call it through run_hhg(), so that each figure is produced
# by a single R script. The Python interpreter is `python3`, or the one given
# by the environment variable MANIFOLDTEST_PYTHON.
# ================================================================

python_bin <- function() Sys.getenv("MANIFOLDTEST_PYTHON", "python3")

hhg_script <- function() {
  system.file("python", "hhg.py", package = "manifoldtest.simulations",
              mustWork = TRUE)
}

# HHG p-values for every replication of a data file:
#   space "hyperbolic": the file written by hyperbolic_shared_data(),
#                       result columns exp, scan, xval, rep, n, p_value
#   space "grassmann" : the file written by figure4_grassmann.R,
#                       result columns scan, xval, rep, n, HHG
# The permutations of a replication are seeded by its scenario key, so the
# p-values do not depend on which other replications are in the file. The
# replications are spread over all but two of the available cores.
run_hhg <- function(input, space = c("hyperbolic", "grassmann"), n_perm = 199,
                    output = tempfile(fileext = ".csv")) {
  space <- match.arg(space)
  status <- system2(python_bin(),
                    c(shQuote(hhg_script()), space, shQuote(input),
                      shQuote(output), "--n-perm", n_perm))
  if (!identical(status, 0L)) stop("the Python HHG run failed (status ", status, ")")
  read.csv(output, na.strings = character(0))
}

# Geodesic distance matrix of the first `nrows` rows of a Grassmann data
# file, computed by the Python pipeline (used to check it against the R one)
python_grassmann_dist <- function(input, nrows, output = tempfile(fileext = ".txt")) {
  status <- system2(python_bin(),
                    c(shQuote(hhg_script()), "grassmann", shQuote(input),
                      shQuote(output), "--distances", nrows))
  if (!identical(status, 0L)) stop("the Python distance run failed (status ", status, ")")
  as.matrix(read.csv(output, header = FALSE))
}
