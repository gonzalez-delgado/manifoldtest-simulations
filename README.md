# Simulations for "Independence test on Riemannian manifolds"

Code that produces the numerical experiments presented in:

> Liu, H., Li, M., González-Delgado, J. and Hallin, M. (2026). Independence
> test on Riemannian manifolds.

The two tests proposed in the paper, **DistCov** and **RankCov**, are run
with the R package [`manifoldtest`](https://github.com/gonzalez-delgado/manifoldtest)
(function `manifold.dcov.test()`). The benchmark procedures (Raw-dCov, Geo-dCov, BallCov,
HHG) are implemented here.

The repository is an R package, `manifoldtest.simulations`, with one R script
per figure and per number quoted in the text.

---

## 1. Figure and number map

Each item of the paper is produced by a single R script, in `inst/scripts/`:

| Item in the paper | Script |
|---|---|
| Figure 1 (Section 1.2) | `figure1_torus.R` |
| Figure 2 (Section 5.1) | `figure2_hyperbolic.R` |
| Figure 3 (Section 5.2) | `figure3_outliers.R` |
| Figure 4 (Section 5.3) | `figure4_grassmann.R` |
| Section 4.3, rejection rate 83 → 75, 71, 32, 6 | `text_basepoint_hyperbolic.R` |
| Section 4.3, torus rejection rate 54 → 5 | `text_basepoint_torus.R` |
| Section 5.4, 100 / 100 / 79 on R², H², T² (N = 100) | `text_matched_precision.R` |
| Section 5.4, sphere antipodal 6.7 / 13.3 / 100 | `text_sphere_antipodal.R` |
| RankCov power versus dimension | `text_rankcov_dimension.R` |
| Checks of the designs (geometry, population V², negative type) | `design_checks.R` |

The real data analyses presented in Section 6 can be reproduced with the code in the `examples/` folder of
the `manifoldtest` repository.

---

## 2. Layout

```
manifoldtest-simulations/
├── DESCRIPTION, NAMESPACE      package metadata
├── R/                          functions shared by the scripts
│   ├── helpers.R               benchmark distance matrices; the sphere (built
│   │                           with manifoldtest::new_manifold); paper_test(),
│   │                           a thin wrapper of manifold.dcov.test
│   ├── hyperbolic.R            samples of Figures 2 and 3 shared by all the
│   │                           procedures; Geo-dCov on them
│   ├── hhg.R                   run_hhg(), which calls inst/python/hhg.py
│   └── plot.R                  power curves of Figures 2 and 3
├── man/                        help page of the package
├── inst/
│   ├── scripts/                one R script per figure or number (section 1)
│   ├── python/hhg.py           HHG test (Heller, Heller and Gorfine, 2013)
│   └── extdata/grassmann_basepoint.csv   fixed base point U0 (4 x 2)
└── run_all.sh                  runs every script, in order
```

The HHG benchmark is written in Python; the R scripts that need it call it
themselves, so each figure is still produced by running a single R script.

---

## 3. Running the simulations

### Install

Install the package from GitHub. This also installs `manifoldtest` and the
other R dependencies:

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

From the root of the clone:

```bash
bash run_all.sh            # results in ./results
bash run_all.sh my_dir     # results in ./my_dir
```

The figures, the tables of rejection rates behind them, and one log per
script (`<script>.log`) are written to the output directory. Running
everything takes a few hours on a multi-core machine.

### Run a single figure or number

Each script writes its results to the **working directory**. From the root
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

The scripts are independent of each other and can be run in any order or in
parallel. `figure4_grassmann.R` and the HHG step of Figures 2–4 use all but
two of the available cores.

---

## 4. Requirements

- R (>= 4.0) with `manifoldtest` (>= 1.0.0), `energy`, `Ball`, `data.table`,
  `ggplot2`, `reshape2` and `gtable`, all installed with the package.
- Python 3 with `numpy` and `pandas`, for the HHG benchmark. The scripts
  call `python3`; set the environment variable `MANIFOLDTEST_PYTHON` to use
  another interpreter, e.g. `MANIFOLDTEST_PYTHON=/path/to/venv/bin/python`.
