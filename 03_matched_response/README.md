# Matched 130/140 MHz beam, noise and seven-model response

## What was reviewed

The uploaded seven130_input_review.zip reports 42 completed jobs: seven models, three seeds, coeval and real-space lightcone. Independent checks of its metadata and geometry passed: all requested configurations match the declared reference; within-seed coeval density fingerprints agree; the three seed fingerprints differ; 128 channel centres (126.03125–133.96875 MHz) are covered by 127 physical radial slices without coeval wrapping. The geometry archive hash and target-distance interpolation were checked. The original large signal arrays were not uploaded; preparation verifies their archived array hashes locally.

This package performs the next instrumental stage. It does NOT rerun 21cmFAST, modify the paper, or modify older calculations.

## Commands: run sequentially

Put this ZIP in ~/eor_model_test, then:

```bash
cd "$HOME/eor_model_test"
unzip eor_matched_response_v1.zip
bash eor_matched_response_v1/run.sh check
bash eor_matched_response_v1/run.sh discover
```

Discovery reads your existing seven_models_signal_inputs*.zip exports and verifies all 14 requested configurations in each. It requires all six frequency/seed combinations. Identical duplicate exports are accepted; conflicting versions cause an error. The existing ska_station_sensitivity_AAVS2.h5 must be somewhere under the project directory. No ripgrep command is needed. Inspect the printed paths and matched_response_inputs.json.

Next calculate the noise curves and prepare all six sets of seven models:

```bash
bash eor_matched_response_v1/run.sh prepare
```

This uses the large exports already on your disk. In particular, retain:

- seven_models130_three_seeds_v1/seed12345/seven_models_signal_inputs.zip
- seven_models130_three_seeds_v1/seed23456/seven_models_signal_inputs.zip
- seven_models130_three_seeds_v1/seed34567/seven_models_signal_inputs.zip
- the existing 140 MHz seed12345 seven-model export
- seven_seeds140_v1/seed23456/seven_models_signal_inputs.zip
- seven_seeds140_v1/seed34567/seven_models_signal_inputs.zip

The small review ZIP alone does not contain the signal cubes.

Then evaluate the beam responses:

```bash
bash eor_matched_response_v1/run.sh beams
```

At 130 MHz this runs OSKAR on the exact new 16,384 signal directions plus the original 24 foreground directions and pointing centre, at three original time midpoints and all 128 channel frequencies. It retains raw Jones outputs locally as well as normalized unpolarized power. It computes 768 dipole/isotropic jobs and nine Gaussian regression jobs per newly computed frequency. Gaussian power is evaluated analytically at the remaining channels. This stage may take substantial time; a one-channel OSKAR call can be silent until it finishes.

The wrapper requests reuse of the original 140 MHz beams if both recorded beam/prepared paths exist. Reuse requires exact equality of directions, solid angles, frequency and radial geometry, matching provenance, all 777 records and verified hashes. If those files have moved, edit BP and BEAMS in run.sh. If unavailable, the wrapper computes 140 MHz beams again. A validation failure stops the run; it does not silently accept a near match.

Project all seven models at both frequencies, apply the same estimator, and calculate covariance/bias diagnostics:

```bash
bash eor_matched_response_v1/run.sh analyze
bash eor_matched_response_v1/run.sh package
```

Send:

```
~/eor_model_test/matched_response130_140_v1/matched_response_review.zip
```

Keep physical.npz, beam_power.npz, individual Jones files and original signal exports locally. The review archive omits these bulky inputs but includes visibilities, powers, all21-pair contrasts, conditional covariances, noise realizations, compact geometry, configurations and hashes.

## Safe resumption

Run one copy at a time. Completed preparation, beams and result stages have checks before reuse. Completed OSKAR jobs are reused individually only if saved settings and response hashes match. A partial prepare directory is not silently overwritten: preserve the failed directory and choose a new OUT path in run.sh, or diagnose the failure before explicitly removing only that partial prepared directory. An existing discovery manifest is retained; do not rerun discover unless input locations need changing. Capture any traceback and the relevant stage log if a check fails.

If two different versions of a signal export are found, do not choose the newest by timestamp. Use the version from the reviewed calculation. You can make a directory of symlinks to only the six intended archives and the intended SEFD file, and call pipeline.py discover --root that_directory using the same container. The script validates configurations again when preparing.

## Definitions preserved between frequencies

- Seven models: F, Hminus, Hplus, Mminus, Mplus, Iminus, Iplus.
- Seeds: 12345, 23456, 34567; 14 coeval/lightcone fields per seed.
- The same 16 selected station pairs, four groups, three MJD60000 snapshots, fixed celestial catalogue and phase centre.
- 128 channels, each 62.5 kHz; 10-second integration per independent noise split.
- A fixed physical reference-station layout, fixed dipole length c/(2×140 MHz), and Gaussian reference width 4 degrees at 140 MHz. Changing observing frequency does not resize the physical dipole.
- Linear radial interpolation; transverse mean removed per slice; full 128²-pixel, 150 cMpc square; no tiling or extrapolation.
- Fixed catalogue, weighted degree1 and degree2 scaled-log-frequency fits; then unweighted channel mean removal, RMS-normalized Hann and unitary FFT.
- Absolute-delay edges [0, .0625, .25, .5, 1, 2, 4, 8.01] microseconds. The final bin remains available but is flagged outside the previous six-bin support.
- Equal final averaging of records and modes; inverse noise variance is used for the polynomial fits. Adding frequency-dependent inverse-variance averaging to the final power bins would define a different estimator and is not done here.

## Beam and noise conventions

This continues the existing shared S8-1 reference-station calculation. All selected stations use the same beam shape; dipole raw Jones responses are converted to trace(JJ†)/2, normalized to the pointing centre. It is NOT a station-specific SKALA response or a calibrated full SKA observation. Frequency-dependent response is evaluated at the actual signal directions and catalogue positions. Beam evolution within each integration/channel remains neglected; the existing fringe sinc factors are retained.

The local HDF5 lookup follows the supplied SKAO v11.1.1 prescription, selecting one nearest LST/Az/ZA entry at each original midpoint and interpolating X/Y SEFD across frequency. It retains separate time curves:

    SEFD_I = 0.5 sqrt(SEFD_X^2 + SEFD_Y^2)
    Var_complex(I) = SEFD_I^2 / (delta_nu * t_split)
    Var_real = Var_imag = Var_complex / 2
    fitting_weight = 1 / Var_complex

There is no additional polarization factor. These are single-station-equivalent values for identical stations, not array SEFDs; there is no division by station count. The table's AAVS2 model and the normalized reference beam are a conditional hybrid sensitivity model, not an independent absolute beam/sensitivity calibration.

Source: https://developer.skao.int/projects/ska-ost-senscalc/en/11.1.1/_modules/ska_ost_senscalc/low/sefd_lookup.html

## Outputs and interpretation

- noise/noise.npz: per-time/per-channel X/Y and Stokes-I SEFD, complex variance, inverse-variance fit weights.
- noise/configuration.json: table hash, selected table bins, actual pointing/LST, units and normalization.
- beams/: locally retained raw Jones files, normalized response, settings and hashes.
- seed*/results/results.npz: visibilities; power[beam,field,pipeline,28]; all21 pair power differences; foreground residual bias including the signal/residual cross term; thermal plus fixed-sky mixed covariance.
- noise_validation/: 2,048 independent noise-pair realizations per frequency, analytic non-diagonal covariance, empirical validation diagnostics and operators. Realizations are independent within each split pair; the same draws across pipeline variants are paired comparisons.
- all_pair_contrasts.csv: absolute and normalized contrasts for all frequencies, seeds, beams, estimator variants and bins. The 1e-3 power-sum reporting floor is a descriptive support convention, not a detection threshold.

The Gaussian-model residual uses dipole truth minus Gaussian-model foregrounds, preserving the previous mismatch diagnostic. It is a known finite synthetic catalogue without diffuse emission or gain solving. Analytic covariance includes thermal noise and fixed-sky mixed terms; it is not cosmological sample variance or a foreground ensemble covariance. Independent noise between baselines/times/channels is assumed before applying the estimator; filtering/windowing create within-group delay-bin correlations.

The two bands share a fixed comoving square, not exactly the same angular window. Equal delay bins also correspond to different cosmological k_parallel. Report this as a fixed-comoving-window comparison; do not claim a common-angular-support frequency optimum without a further matched-support check. Do not add these bands' SNR² as independent observations. No distance, covariance inversion, likelihood or 1000-hour scaling is performed.

## Checks performed before delivery

- All42 uploaded requested configurations and the supplied130 geometry were independently checked.
- Local SEFD computation ran on the supplied HDF5 table for both frequencies. The140 complex variances reproduce the archived values with max relative error 4.44e-16.
- Station-equivalent Stokes-I SEFD ranges:130 MHz 1957.37–3669.87 Jy;140 MHz 1864.95–3547.90 Jy across time/channels.
- 2,048 noise-pair realizations were run per frequency. Empirical/analytic diagonal covariance ratios across the three estimator variants span approximately0.921–1.082 at130 and0.950–1.056 at140. These finite-Monte-Carlo diagnostics require review; they do not validate the full sky covariance.
- Self-tests cover operator identity, exact averaging conventions, noise and mixed covariance scaling, PSD and frequency-dependent response.
- A synthetic full-size seven-model130 preparation reproduced the uploaded geometry exactly and removed slice means as intended. Synthetic fixtures are not included as physical results.
- A full128-channel synthetic seven-model projection reproduced identical results across beam branches supplied with identical test beams.
- Actual OSKAR runs and physical seven-model projection are to be executed locally; they were not performed here. The140 physical results must pass a reproduction gate for every seed using the archived powers and foregrounds.

local_validation/ contains the actual table-based noise calculations and diagnostic summaries from these tests. They are reference checks, not a replacement for the new physical beam/model runs.
