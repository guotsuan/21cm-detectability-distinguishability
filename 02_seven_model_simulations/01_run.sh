#!/usr/bin/env bash
set -euo pipefail
ROOT="$HOME/eor_model_test"
PKG="$ROOT/eor_seven130_v1"
BATCH="$ROOT/seven_models130_three_seeds_v1"
conda run --no-capture-output -n eor-models python "$PKG/self_test.py"
# Geometry preflight does not generate physical cubes.
conda run --no-capture-output -n eor-models python "$PKG/batch130.py" preflight \
 --out "$BATCH/geometry_preflight.json"
if [[ ! -f "$BATCH/batch.json" ]]; then
 conda run --no-capture-output -n eor-models python "$PKG/batch130.py" plan --batch "$BATCH"
fi
conda run --no-capture-output -n eor-models python "$PKG/batch130.py" run \
 --batch "$BATCH" --execute 2>&1 | tee -a "$BATCH/generation.log"
conda run --no-capture-output -n eor-models python "$PKG/batch130.py" review \
 --batch "$BATCH" 2>&1 | tee -a "$BATCH/review.log"
printf '\nSend this file:\n%s/seven130_input_review.zip\n' "$BATCH"
