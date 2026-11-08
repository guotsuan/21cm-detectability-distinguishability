# 130 MHz extension: seven models, three matched seeds

This package performs stage A of the matched physical-frequency extension: generate and audit the 130 MHz cosmological inputs and export the signal geometry for the next instrumental stage. It does not yet compute130-MHz visibilities or a130-versus140 discrimination ranking.

## Run

Save the ZIP to ~/eor_model_test and execute:

```bash
cd "$HOME/eor_model_test"
unzip eor_seven130_v1.zip
bash eor_seven130_v1/01_run.sh
```

Uses `conda run -n eor-models` and the previously used py21cmfast3.4.0. No installation or environment upgrade is required. No OSKAR command runs at this stage.

The script runs configuration self-tests and a lightcone geometry preflight, creates42 jobs, executes them sequentially, audits all outputs and exports per-seed signal ZIPs. The output root is `~/eor_model_test/seven_models130_three_seeds_v1`. If interrupted, rerun the same command: completed files are checked against their configuration and checksum before skipping, and partial simulations reuse the worker's isolated caches. Do not run concurrent copies. Existing140-MHz results are untouched. This package deliberately generates a fresh controlled130-MHz batch and does not search for or automatically reuse older130-MHz runs whose configurations may differ.

## Optional preflight only

This creates no physical cubes:

```bash
cd "$HOME/eor_model_test"
conda run --no-capture-output -n eor-models python eor_seven130_v1/batch130.py preflight \
 --out seven_models130_three_seeds_v1/geometry_preflight.json
```

The preflight uses the installed3.4.0 RectilinearLightconer to inspect actual physical slice distances and coeval slice selection. All128 target channel centres must be covered (126.03125–133.96875MHz), with increasing comoving distances and no repeated coeval slices. If it fails, preserve/send the complete terminal error; do not enlarge the box, allow extrapolation, or bypass the check. API execution cannot be checked here because the3.4.0 environment is on your machine.

## Fixed configuration

- Models: F,Hminus,Hplus,Mminus,Mplus,Iminus,Iplus.
- Seeds:12345,23456,34567, each shared among all seven models.
- Both coeval and real-space lightcone outputs:7x3x2=42 jobs.
- Same150cMpc box,HII_DIM128,DIM384 and four threads as the accepted140-MHz batch.
- Same cosmology, astrophysical templates, real-space flags, evolution anchor and lightcone construction procedure.
- Changes from the archived140-MHz requested configurations are only centre frequency, centre redshift and the declared seed. The130-MHz redshift is computed as1420.40575177/130-1 (not silently matched to a rounded redshift from unrelated older runs).
-128 channels of62.5kHz over126–134MHz.

Coeval density fingerprints must match among all seven models WITHIN each seed and must differ between the three seeds. Density arrays at130 and140MHz are not required to be byte-identical: they correspond to different redshifts. The same seed and initial-grid configuration preserve the intended phase matching.

## Outputs to return

Send only:

```
~/eor_model_test/seven_models130_three_seeds_v1/seven130_input_review.zip
```

It contains input audit reports, requested/resolved configurations, environment metadata, provenance, the geometry preflight and `signal_geometry130.npz`. The latter contains actual signal directions, solid angles,128 frequencies, source/target comoving distances and source-slice indices.

Keep the large cubes and three `seed*/seven_models_signal_inputs.zip` archives locally. Those signal archives are needed for the subsequent forward calculation and are not included in the small review.

## What follows this stage

After reviewing the geometry and all42 outputs:

1. Recompute the Jones/station response at130-MHz channel frequencies and these signal directions, together with the same foreground catalogue and observing times.
2. Obtain/validate130-MHz SEFD/noise inputs for the same observation, keeping polarization and complex-variance conventions explicit. Reusing140-MHz variances would be a fixed-noise diagnostic, not the requested physical frequency comparison.
3. Preserve the mean/Hann/FFT, fitting prescriptions, baseline selection, bin definitions and reporting masks. Recompute the wavelength-dependent fringe and averaging response.
4. Project all seven models, compare all21 pairs across matched seeds, and compare the two frequencies using an explicitly stated angular/comoving support convention.

A fixed150cMpc transverse box subtends different angles at130 and140MHz. The new review reports the angular span; do not call the two angular windows identical. A fair frequency interpretation must either show common angular support or label the fixed-comoving-window comparison and its limitations. Matching delay bins also does not automatically match cosmological k_parallel across frequencies. No SNR values from overlapping bands should be combined as independent observations.

Do not send these130-MHz signals into the old140-MHz prepare/analyze tools by merely changing filenames. This package intentionally does not attach stale140-MHz beam/noise arrays. The final physical-frequency comparison remains incomplete until the response and noise stages are finished.

## Validation before delivery

The numerical worker is the same one used successfully for the140-MHz seven-model batch. The generic verifier was changed only to compute coverage endpoints from each job's frequency rather than hard-code140MHz. Self-tests verify that only the declared frequency/redshift/seed changes occur, all42 planned configurations are distinct and correct, dry-run scheduling works, incomplete-band and duplicate-slice cases are rejected, and the vendor verifier accepts a synthetic130-MHz band. Python compilation and shell syntax passed. No actual130-MHz simulation or installed3.4.0 geometry preflight was executed here. The shell command runs that preflight before any expensive job.

No paper source is changed by this package. Results and frequency-ranking conclusions should be added only after the full130-MHz forward comparison is reviewed.
