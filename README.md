# Detectability versus Distinguishability in 21-cm Model Comparisons

Code and selected derived data accompanying **Detectability versus Distinguishability: How Foreground Mode Loss and Heating Degeneracies Shape 21-cm Model Comparisons**, by Shunwu Xiao, Quan Guo and Qian Zheng (manuscript in preparation for submission).

This initial release is a curated research-code snapshot, not a complete one-command reproduction of every paper figure or of the full physical simulations. `FIGURE_COVERAGE.csv` gives the exact plotting coverage. The reference manuscript source at curation is Overleaf Git commit `0cbab1f` (the subsequent data-availability edit does not change calculations).

## Contents

| Path | Purpose |
| --- | --- |
| `00_workflow/` | Explicit fixed-heating follow-up stages; set `EOR_DATA_ROOT` before use |
| `01_frequency_forecasts/` | Spectrum extraction, 29-band instrument weights and template-distance forecasts |
| `02_seven_model_simulations/` | Matched seven-model 130 MHz production runner and reference configurations |
| `03_matched_response/` | Two-frequency beam/noise/visibility processing; reference inputs and synthetic tests |
| `04_heating_followup/code/` | Fixed heating-candidate simulation and analysis code |
| `06_visibility_distance/` | Full conditional covariance F–Mplus output-distance analysis and reference tables |
| `paper_snapshot/` | Selected final-figure scripts with their required compact derived inputs |
| `derived_tables/` | Additional spectra, forecast and interpretation tables |
| `SOURCE_MANIFEST.json` | Original source locations and SHA-256 before portability edits |
| `CHECKSUMS.json` | Checksums of the curated release files |

`matched_response_review.zip` is distributed separately as a **GitHub Release asset**. It contains six frequency/seed result sets, processed visibilities, powers, conditional covariance products, noise diagnostics and frozen estimator code. It is required to recompute `06_visibility_distance`.

## Lightweight checks and selected-figure replay

Use a dedicated Python environment. The light analysis needs NumPy and Matplotlib; `requirements-analysis-tested.txt` records the versions tested during curation, not the production simulation environment.

```bash
python -m pip install -r requirements-analysis-tested.txt
python 01_frequency_forecasts/self_test.py
python 02_seven_model_simulations/self_test.py
python 03_matched_response/self_test.py
python 04_heating_followup/code/self_test.py
python reproduce_selected_figures.py --out /tmp/21cm-paper-figures
```

The replay command requires a new output directory, copies the compact input snapshot there, and writes five figures. It does not run simulations or edit your manuscript. Historical script names contain old figure numbers; use `FIGURE_COVERAGE.csv` and the generated publication filenames instead.

To recompute the conditional output distances, download `matched_response_review.zip` from the release page, then run:

```bash
python 06_visibility_distance/analyze_visibility.py \
  --archive /path/to/matched_response_review.zip \
  --out /tmp/21cm-output-distances
```

## Physical production environments and external inputs

The manuscript uses **21cmFAST 3.4.0**, **21cmSense 3.0.0**, **OSKAR 2.13.0**, and **ska-ost-array-config 5.0.0** in separate production stages. Do not substitute the analysis-only environment for those environments. The code checks versions and reference inputs. Compiled libraries and the OSKAR container are not included.

A full rerun additionally needs the original coeval/lightcone cubes, beam products, saved forecast/reference-array inputs and the station sensitivity HDF5 table. These are not all included. In particular, the upstream release of the historical station-sensitivity table was not recorded, so downloading an arbitrary current table is not a verified replacement. See `REPRODUCIBILITY.md`.

Module README files retain historical execution examples and may mention earlier project folders or unfinished work. Their historical status statements do **not** describe the final manuscript: the fixed 130 MHz heating results are now included in `paper_snapshot`, and the later output-distance analysis is in `06_visibility_distance`. Use the present README and `VALIDATION.json` for this release's coverage. Original wrapper defaults may require explicit local paths; inspect them before production execution.

## Interpretation

The forecast and visibility routes have different geometry, exposure and estimators. Distances are conditional comparisons, not calibrated exclusion probabilities. Three seeds describe repeatability, not an estimated cosmological covariance. The release does not add a multi-band likelihood, RSD convergence, or an independently calibrated absolute SKA sensitivity.

## Citation and rights

See `CITATION.cff`. No journal DOI has been assigned. Author-contributed software is released under MIT. Third-party software and externally sourced instrument products retain their own rights and must be obtained under their upstream terms; see `THIRD_PARTY_NOTICES.md`.
