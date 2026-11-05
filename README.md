# manifoldtest-simulations

Code that produces the numerical experiments presented in:

> Liu, H., Li, M., González-Delgado, J. and Hallin, M. (2026). Independence
> test on Riemannian manifolds.

The two tests proposed in the paper, **DistCov** and **RankCov**, are run
with the R package [`manifoldtest`](https://github.com/gonzalez-delgado/manifoldtest)
(function `manifold.dcov.test()`). The benchmark procedures (Raw-dCov, Geo-dCov, BallCov,
HHG) are implemented here.

## Running the simulations

### Install

To reproduce the numerical experiments, first install `manifoldtest-simulations` from GitHub. This also installs `manifoldtest` and the
required R dependencies:

```r
install.packages("remotes")
remotes::install_github("gonzalez-delgado/manifoldtest-simulations")
```

or, from a local clone of the repository:

```bash
git clone https://github.com/gonzalez-delgado/manifoldtest-simulations
cd manifoldtest-simulations
R CMD INSTALL .
```

Figures 2, 3 and 4 also need Python 3 with `numpy` and `pandas` (section 4).

### Run everything

To reproduce all the numerical experiments presented in the paper, the user
can run the script `run_all.sh` from the root of the clone:

```bash
bash run_all.sh            # results in ./results
bash run_all.sh my_dir     # results in ./my_dir
```

The figures, the tables of rejection rates behind them, and one log per
script (`<script>.log`) are written to the output directory. Running
everything takes a few hours on a multi-core machine.

### Run a single figure or number

To reproduce a single figure or numerical result, the user can
run its corresponding R script, that writes its results to the working directory. From the root
of the clone:

```bash
mkdir -p results && cd results
Rscript ../inst/scripts/figure2_hyperbolic.R
```

The scripts are also installed with the package, so they can be run
without a clone, for instance from R:

```r
setwd("my_results")
source(system.file("scripts", "figure2_hyperbolic.R",
                   package = "manifoldtest.simulations"))
```

The individual R scripts are available at `inst/scripts/`:

| Item in the paper | Script |
|---|---|
| Figure 1 (Section 1.2) | `figure1_torus.R` |
| Figure 2 (Section 5.1) | `figure2_hyperbolic.R` |
| Figure 3 (Section 5.2) | `figure3_outliers.R` |
| Figure 4 (Section 5.3) | `figure4_grassmann.R` |
| Section 4.3, rejection rate 83 → 75, 71, 32, 6 | `text_basepoint_hyperbolic.R` |
| Section 4.3, torus rejection rate 54 → 5 | `text_basepoint_torus.R` |
| Section 5.3, distance covariance 6.3×10⁻⁴ / 1.3×10⁻² / 1.0×10⁻³ (N = 1500) | `text_grassmann_population_dcov.R` |
| Section 5.4, sphere antipodal 6.7 / 13.3 / 100 | `text_sphere_antipodal.R` |
| Section 5.4, 100 / 100 / 79 on R², H², T² (N = 100), 100 on all three (N = 300) | `text_matched_precision.R` |

The scripts are independent of each other and can be run in any order or in
parallel. `figure4_grassmann.R` and the HHG step of Figures 2–4 use all but
two of the available cores.

The real data analyses presented in Section 6 can be reproduced with the code in the
[`examples/`](https://github.com/gonzalez-delgado/manifoldtest/tree/main/examples) folder of
the `manifoldtest` repository.

### Requirements

- R (>= 4.0) with `manifoldtest` (>= 1.0.0), `energy`, `Ball`, `data.table`,
  `ggplot2` and `reshape2`, all installed with the package.
- Python 3 with `numpy` and `pandas`, for the HHG benchmark. The scripts
  call `python3`; set the environment variable `MANIFOLDTEST_PYTHON` to use
  another interpreter, e.g. `MANIFOLDTEST_PYTHON=/path/to/venv/bin/python`.
