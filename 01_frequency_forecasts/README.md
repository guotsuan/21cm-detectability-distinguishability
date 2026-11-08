# Frequency-dependent EoR forecasts — staged execution

This package continues the seven-model 96^3/128^3, DIM=384 analysis. It reads existing coeval cubes; it does not run 21cmFAST or change simulation outputs. All commands below run in Bash on SKA-MILAN from `~/eor_model_test`.

## Scope and limits

- Extract all 434 spectra: 31 frequencies × seven models × two resolutions. Compare complete physical bins within 0.1–0.3 Mpc^-1; retain surrounding bins for interpolation.
- Rebuild SKA weights for 29 analysis bands centred 55–195 MHz, width 8 MHz. This keeps every full band within 50–200 MHz. Spectra at 50 and 200 MHz remain useful evolution diagnostics.
- Reproduce both the original native 130 MHz instrument weights and the final seven-model 96/128 forecasts before extending across frequency.
- Compute detection and 21 pair distances per band, three masks, three kmax cuts. Use the same saved 21cmSense 3.0.0 public covariance convention: total standard error = thermal standard error + alpha × signal. No extra factor of two, polarization factor, days factor or simulation mode count is applied.
- Zero signal gives zero SNR when modes exist. No retained modes gives an empty metric. No positive artificial power floor is introduced. Fully zero spectra stay zero; positive spectra use the previous log-log interpolation; mixed zero/positive spectra use linear interpolation across all populated centres and are labelled. Ratios with zero denominators are blank.
- The map uses 2|P128−P96|/(P128+P96), not the high/low ratio. The CSV also supplies that ratio. Both-zero bins get a zero symmetric difference with a `both_zero` label; this is not evidence of convergence. No bins are silently excluded based on a chosen tolerance.
- Forecasts retain the real-space, coeval, binary-mask assumptions. Resolution across frequency is checked in power; density-phase evolution is not independently revalidated by this package.
- Instrument calibration here means code-level regression and internal algebra checks. It is NOT an independent physical sensitivity benchmark, measured instrumental response, or validated cross-frequency cosmological covariance.
- Each band uses the original observing-time allocation (normally 1,000 hours). This is a conditional per-band assumption, not a claim that the telescope observes the full 50–200 MHz simultaneously for 1,000 hours.
- No combined-band significance is produced. The response/overlap experiment is illustrative and must not be inserted into the forecast covariance.

## Existing inputs

Expected relative paths (overridable via CLI):

```
matched_130MHz_resolution/res96d384_seed12345_real/<model>/
matched_130MHz_resolution/res128d384_seed12345_real/<model>/
    requested_config.json, COMPLETE.json, history.csv, coeval_z*.npz
forecast_130MHz_full/
    forecast_metadata.json, forecast_native.npz, detection_snr.csv
grid_comparison_96_128_all7/
    detection_comparison.csv, pairwise_comparison.csv
<original SKA reference directory>/
    reference_observation.json, station_positions_native_xyz_m.npy
```

Models: F, Hminus, Hplus, Mminus, Mplus, Iminus, Iplus. Cubes use the verified existing keys `brightness_temp`, `xH_box`, `redshift`, `box_len_cMpc`, `random_seed`, and `coordinate_space`. The scripts check configuration, exact stored redshift, geometry, source checksums, history/cube consistency and Parseval normalization. Requested astrophysical settings must match the original forecast metadata. Both grids must retain DIM=384.

The SKA script automatically locates one matching reference directory beneath the working directory. If no unique match exists, add `--reference /exact/path/to/reference_directory` to BOTH the preflight and calculation commands. Do not recreate the station file or substitute a different array layout.

## 1. Extract and test the package

Assuming the downloaded ZIP is in `~/Downloads`:

```bash
conda activate eor-models
cd ~/eor_model_test
unzip ~/Downloads/eor_frequency_pipeline.zip -d .
python -I eor_frequency_pipeline/self_test.py
```

NumPy and Matplotlib are required. The self-test uses temporary synthetic cubes and removes them afterward; it needs neither 21cmFAST nor 21cmSense. It includes Fourier normalization, zero-power handling, extrapolation rejection, exact spectral-response checks, analytic covariance metrics, command-line extraction/forecast runs and resume rejection after changed inputs.

## 2. Check existing simulations without calculating FFTs

```bash
python -I eor_frequency_pipeline/extract_spectra.py --preflight
```

Expected: all 434 cube paths found, matching configurations, exact frequency schedules. Fix missing or mismatched inputs before continuing. This is a path/configuration check; actual array checks follow during extraction.

## 3. Extract spectra and produce the resolution map

```bash
mkdir -p frequency_logs
/usr/bin/time -v -o frequency_logs/extraction_resources.txt \
    python -I -u eor_frequency_pipeline/extract_spectra.py \
    > frequency_logs/extraction.log 2>&1
status=$?
echo "Exit status: $status"
tail -n 20 frequency_logs/extraction.log
```

Continue only with exit status 0. In another terminal, `tail -f frequency_logs/extraction.log` monitors progress. Ctrl+C in the tail terminal stops only monitoring.

If interrupted, reuse the exact command with `--resume` and a different log name if you want to keep the old log. The scripts never skip blindly: source cube and output checksums must agree. Changes in inputs, code or settings require a new output directory. There is no overwrite/force option.

Review:

```
frequency_analysis/spectra/resolution_map.png
frequency_analysis/spectra/resolution_map.csv
frequency_analysis/spectra/snapshot_summary.csv
frequency_analysis/spectra/spectra.csv
```

An all-zero spectrum alone does not prove full ionization; inspect `xHI_volume`. Use the map to identify sensitive epochs, especially heating transitions. A small grid difference is not a numerical error bound. Current default forecasts still include the stated k range for diagnostics; do not claim uniform precision across it.

## 4. Validate the original SKA environment and reference

```bash
conda activate ska-forecast
cd ~/eor_model_test
export PYTHONNOUSERSITE=1
python -I eor_frequency_pipeline/build_noise.py --preflight
```

Use your existing pinned environment. The code requires 21cmSense 3.0.0 and compares installed station coordinates against the saved file. It prints installed versions. Do not upgrade packages to bypass a failed check. If a dependency is missing, restore the previously working environment first.

## 5. Rebuild ONLY the 130 MHz instrument weights

```bash
python -I -u eor_frequency_pipeline/build_noise.py --frequencies 130
```

Expected `frequency_analysis/noise/130MHz/COMPLETE.json` contains `regression_130.status: PASS`. The test compares physical native coordinates, thermal errors, sample coefficients and all masks against your saved original forecast. Duplicate native coordinate entries are preserved. Constant and scale-dependent test spectra additionally verify sample-error scaling and the saved public combined-error convention at each frequency.

A mismatch stops the run. Send the full traceback; do not relax the tolerance or continue to other bands. Matching the old implementation is a reproducibility test, not an independent proof of its physical normalization.

## 6. Reproduce the final seven-model 130 MHz results

```bash
python -I -u eor_frequency_pipeline/forecast_bands.py --frequencies 130
```

This compares both resolutions' detection metrics and all exported pairwise covariance conventions against `grid_comparison_96_128_all7`. Expected `frequency_analysis/forecasts/130MHz/COMPLETE.json` contains `regression_130.status: PASS`.

At kmax=0.30 Mpc^-1 and horizon+buffer, the previously reviewed 128^3 F SNR was approximately 12.737, F–Mplus distance 1.938 and Mminus–Iminus distance 1.743. The CSV comparison, not these rounded numbers, is the actual pass criterion.

Optional checkpoint for review:

```bash
python -I eor_frequency_pipeline/package_results.py \
    --stage pilot --out frequency_130MHz_regression_review.zip
```

## 7. Recalculate the 29 observing bands

Both 130 MHz checks should pass before this step. Execute sequentially:

```bash
BANDS=( {55..195..5} )
/usr/bin/time -v -o frequency_logs/noise_resources.txt \
    python -I -u eor_frequency_pipeline/build_noise.py \
    --frequencies "${BANDS[@]}" --resume \
    > frequency_logs/noise_all_bands.log 2>&1
status=$?
echo "Exit status: $status"
tail -n 20 frequency_logs/noise_all_bands.log
```

The existing 130 MHz calculation is checked and reused. Each other band receives a fresh observation, frequency-specific UV weights, temperature, cosmological conversion, beam boundary and horizon buffer. Calculation is serial to limit peak memory; resources depend on frequency. There is no runtime promise.

If interrupted, repeat with `--resume`. Completed bands are checksum-checked and skipped; incomplete bands are recalculated. Do not run two writers against the same output directory.

## 8. Compute all per-band forecasts

```bash
BANDS=( {55..195..5} )
python -I -u eor_frequency_pipeline/forecast_bands.py \
    --frequencies "${BANDS[@]}" --resume
```

Main outputs:

```
frequency_analysis/forecasts/detection_all_bands.csv
frequency_analysis/forecasts/pairs_all_bands.csv
frequency_analysis/forecasts/template_shifts_all_bands.csv
frequency_analysis/forecasts/scale_information_all_bands.csv
frequency_analysis/forecasts/detection_vs_frequency.png
```

Expect 1,827 detection rows and 5,481 pair rows (29 × 63 and 29 × 189). Pairwise columns label F, model-a, model-b, and common-low-resolution-F covariance. SNR is conditional template-amplitude precision; pair distances are not calibrated sigma exclusions or Bayes factors. The numerical-template shift remains a diagnostic, not a fitted error covariance.

`scale_information_all_bands.csv` records squared-SNR contributions in common k intervals. These contributions add within a given band/case/model under the adopted diagonal convention. They do not justify adding overlapping bands.

## 9. Run the spectral-window and band-overlap diagnostic

```bash
python -I eor_frequency_pipeline/response_diagnostics.py
```

Compares rectangular, Hann, and hypothetically flagged Hann windows. It exports the complex Fourier transfer matrix, power mixing, white-noise covariance, Monte Carlo comparison and leakage of a unit-normalized smooth spectral function. Window RMS normalization is recorded; it is not a signal-loss correction. This test does not model full interferometric chromaticity or a calibrated foreground sky.

`band_overlap.csv` shows shared bandwidth and correlations for a toy white complex-Gaussian amplitude/periodogram. Adjacent 8 MHz bands spaced by 5 MHz share 3 MHz, but their true cosmological bandpower covariance is NOT determined by that fraction. Do not apply this toy correlation to the forecast uncertainties. The per-band forecasts continue to use the original instrument convention and are not automatically modified by this window experiment.

What remains for a valid joint forecast: define one actual estimator/window, propagate thermal noise through shared channels, simulate signal evolution and common sky/volume correlations, measure/filter response, and form consistent off-diagonal covariance. Validate with enough independent realizations and held-out injections. Reusing the same seed at several frequencies does not supply an ensemble covariance.

An optional `validate_joint_covariance.py` is supplied for a future externally justified covariance. It requires an NPZ with `delta`, `covariance`, and unique Unicode `labels`, plus a provenance JSON documenting units, bands, window and covariance model. It checks symmetry, positive definiteness and conditioning before evaluating a fixed-covariance distance. It does not infer a missing covariance or calibrate significance. See its `--help`; it is NOT a required step now.

## 10. Send results back

```bash
python -I eor_frequency_pipeline/package_results.py \
    --stage full --out frequency_review.zip
```

Upload `frequency_review.zip` and, if useful, the two resource reports in `frequency_logs`. The archive includes derived tables, maps, diagnostics and the small 130 MHz native products, not the simulation cubes. We can then choose trustworthy frequency/k ranges and prioritize physics/response validation before revising scientific claims.

## Validation status and provenance

The distributed scripts passed synthetic tests for extraction, numerical normalization, covariance metrics, 130 MHz regression logic, zero signals, multi-band CLI output, exact response operations, and resumability. The real SKA calculation was not executed in the authoring environment, which lacks 21cmSense and Astropy. Steps 4–6 therefore provide required real-environment validation before production. No scientific results for the new frequency bands are bundled.

The numerical core and original covariance checker are vendored from your previously supplied `resolution_followup_130MHz.zip` / `eor_next_diagnostics.zip`. The observation API and masks follow your `forecast_130MHz.py` and `ska_low_mode_loss.py`. This preserves the prior covariance convention; a per-cell joint-weighting alternative remains a separate methodological test.

Official API/background sources checked during preparation:
- https://21cmsense.readthedocs.io/en/v3.0.0/tutorials/understanding_21cmsense.html
- https://21cmsense.readthedocs.io/en/v3.0.0/reference/_autosummary/observatory/py21cmsense.observatory.Observatory.html
- https://21cmsense.readthedocs.io/en/latest/_modules/py21cmsense/sensitivity.html

Documentation pages can show development revisions. Runtime checks intentionally pin the same 3.0.0 implementation as your existing forecast rather than silently switching conventions.
