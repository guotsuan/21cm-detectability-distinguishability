#!/usr/bin/env bash
set -euo pipefail
ROOT="${EOR_PROJECT_ROOT:-$HOME/eor_model_test}"
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIF="${OSKAR_SIF:-$ROOT/oskar_tools/OSKAR-2.13.0-Python3.sif}"
OUT="$ROOT/matched_response130_140_v1"
MANIFEST="$ROOT/matched_response_inputs.json"
run_python() {
 apptainer exec --cleanenv --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 \
  --bind "$ROOT:$ROOT" --bind "$PKG:$PKG" "$SIF" python3 "$@"
}
# Changes are explicit stages. Calling without an argument prints the sequence.
STAGE="${1:-help}"
if [[ "$STAGE" == help ]]; then
 printf 'Run sequentially: bash eor_matched_response_v1/run.sh check\n  then discover, prepare, beams, analyze, package\n'
 exit 0
fi
[[ -f "$SIF" ]] || { printf 'Missing container: %s\n' "$SIF"; exit 1; }
if [[ "$STAGE" == check ]]; then
 run_python -c 'import oskar,numpy,scipy,h5py,astropy,shutil; print("Python dependencies OK"); assert shutil.which("oskar_sim_beam_pattern"), "Missing oskar_sim_beam_pattern"'
 run_python "$PKG/self_test.py"
 exit 0
fi
if [[ "$STAGE" == discover && -e "$MANIFEST" ]]; then
 printf 'Manifest already exists: %s\nInspect it; remove or rename it explicitly before rediscovery.\n' "$MANIFEST"; exit 1
fi
EXTRA=()
if [[ "$STAGE" == beams ]]; then
 BP="$ROOT/physical140_prepared_20261018T083333_3986957"
 BEAMS="$ROOT/physical140_beams_20261018T083333_3986957"
 if [[ ! -d "$BP" && -f "${BP}_review.zip" ]]; then BP="${BP}_review.zip"; fi
 if [[ ! -d "$BEAMS" && -f "${BEAMS}_review.zip" ]]; then BEAMS="${BEAMS}_review.zip"; fi
 if [[ -e "$BP" && -e "$BEAMS" ]]; then
  EXTRA+=(--reuse140-prepared "$BP" --reuse140-beams "$BEAMS")
  printf '140 MHz beam reuse requested; exact geometry and file hashes will be checked.\n'
 else
  printf 'Original 140 MHz beam files not found at recorded paths; recomputing both frequencies.\n'
 fi
fi
run_python "$PKG/pipeline.py" "$STAGE" --root "$ROOT" --manifest "$MANIFEST" --out "$OUT" "${EXTRA[@]}" \
  2>&1 | tee "$ROOT/matched_response_${STAGE}_$(date -u +%Y%m%dT%H%M%S).log"
