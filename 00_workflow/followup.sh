#!/usr/bin/env bash
set -euo pipefail
BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${EOR_DATA_ROOT:?Set EOR_DATA_ROOT to your writable simulation/output directory}"
RUNS="$DATA_ROOT/fixed_heating_130_followup_runs"
TAG="${ANALYSIS_TAG:-paperpack_v1}"
CODE="$BUNDLE/04_heating_followup/code"
stage="${1:-help}"
case "$stage" in
 check) python -c "import py21cmfast as p; assert p.__version__ == '3.4.0', p.__version__; print(p.__version__)"; python "$CODE/self_test.py" ;;
 status) for seed in 23456 34567; do for model in F_control Mplus_LX40.00 Mplus_LX40.25; do f="$RUNS/seed$seed/production/130MHz/$model/real/COMPLETE.json"; if [[ -f "$f" ]]; then echo "PRESENT $f"; else echo "MISSING $f"; fi; done; done ;;
 plan|run) seed="${2:?Supply 23456 or 34567}"; [[ "$seed" == 23456 || "$seed" == 34567 ]] || exit 2; args=(); [[ "$stage" == run ]] && args+=(--execute); mkdir -p "$DATA_ROOT/logs"; python "$CODE/run_scan.py" --out "$RUNS" --seed "$seed" --frequencies 130 --threads 4 "${args[@]}" 2>&1 | tee -a "$DATA_ROOT/logs/fixed130_${seed}_${stage}.log" ;;
 analyze) for seed in 23456 34567; do python "$CODE/analyze_scan.py" --runs "$RUNS" --seed "$seed" --frequencies 130 --out "$DATA_ROOT/fixed_heating_130_s${seed}_analysis_$TAG"; done ;;
 pack) for seed in 23456 34567; do d="$DATA_ROOT/fixed_heating_130_s${seed}_analysis_$TAG"; python -c 'import json,sys; assert json.load(open(sys.argv[1]))["status"]=="PASS_FRESH_CHECKS"' "$d/status.json"; mkdir -p "$d/logs"; for log in "$DATA_ROOT"/fixed_heating_130_s${seed}*.log "$DATA_ROOT"/logs/fixed130_${seed}_*.log "$HOME"/eor_model_test/fixed_heating_130_s${seed}*.log; do [[ ! -f "$log" ]] || cp --backup=numbered "$log" "$d/logs/"; done; done; python "$CODE/package_results.py" --runs "$RUNS" --results "$DATA_ROOT/fixed_heating_130_s23456_analysis_$TAG" "$DATA_ROOT/fixed_heating_130_s34567_analysis_$TAG" --out "$DATA_ROOT/fixed_heating_130_review_$TAG.zip" ;;
 *) echo 'Usage: bash 00_workflow/followup.sh check|status|plan SEED|run SEED|analyze|pack';;
esac
