#!/usr/bin/env bash
# Runs every simulation of the paper, one R script per figure or number.
# Install the package first (see README.md), then, from the root of the
# repository:
#     bash run_all.sh [output directory]        (default: results)
# Figures, result tables and one log per script are written to the output
# directory. The scripts are independent and could also be run in parallel.
# Running everything takes a few hours on a multi-core machine.
set -euo pipefail
SCRIPTS="$(cd "$(dirname "$0")" && pwd)/inst/scripts"
OUT="${1:-results}"
mkdir -p "$OUT"
cd "$OUT"

for s in \
    figure1_torus \
    figure2_hyperbolic \
    figure3_outliers \
    figure4_grassmann \
    text_basepoint_hyperbolic \
    text_basepoint_torus \
    text_matched_precision \
    text_sphere_antipodal \
    text_grassmann_population_dcov
do
  echo "==== $s"
  Rscript "$SCRIPTS/$s.R" 2>&1 | tee "$s.log"
done
