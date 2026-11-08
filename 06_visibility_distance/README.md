# Same-operator conditional visibility distances

This supplementary analysis uses the returned `matched_response_review.zip` arrays, without modifying the manuscript. It compares the original F and Mplus **dipole-beam lightcones**, separately at 130 and 140 MHz and for seeds 12345, 23456 and 34567.

## Definition and scope

For each frequency, seed and processing treatment, let `delta = P_F - P_Mplus` be the difference between the two separately processed signal powers. This is not the power of the difference visibility. The reference covariance is always the archive's **F_lightcone** covariance, propagated separately through the same fixed-catalogue, linear-fit or quadratic-fit operator as its corresponding signal. The covariance is held common between the two model means within each comparison; it is not the sum of the two model covariances.

The primary result uses `matched_catalogue` (residual index 0): the dipole signal with zero foreground residual. Define the conditional distance `D = sqrt(delta.T @ inverse(C_F) @ delta)`. This is a distance using a fixed covariance, rather than the full likelihood ratio between model-dependent covariances.

There are four baseline groups and seven delay bins. The last bin, [4, 8.01) microseconds, is excluded in every group. The retained covariance is the 24-by-24 principal submatrix of the archived 28-by-28 covariance **before** inversion: excluded outputs are marginalized, not conditioned upon. All retained off-diagonal terms are kept.

Each noise split has 10 seconds per sampled geometry. Covariance contains thermal and fixed-sky sky–noise mixed terms only. There is no cosmological sample variance, gain/calibration uncertainty, foreground-error ensemble covariance, or 1000-hour extrapolation. Seeds are separate fixed-sky cases, not an estimated covariance ensemble. No rejection probability or general information-retention fraction follows. The two frequencies have different angular and cosmological-mode windows, so comparisons below concern treatments **within** a frequency and seed; they do not establish an observing-frequency ranking.

## Primary results

Ranges below are the minimum and maximum over the three matched seeds, not confidence intervals. Ratios are computed seed by seed before taking their range.

| MHz | treatment | conditional D | D / fixed D | D squared / fixed D squared |
| --- | --- | --- | --- | --- |
| 130 | fixed catalogue | 1.1062e-6–1.7837e-6 | 1 | 1 |
| 130 | linear | 0.7053e-6–1.6963e-6 | 0.6376–0.9510 | 0.4066–0.9044 |
| 130 | quadratic | 0.6020e-6–0.6924e-6 | 0.3390–0.6259 | 0.1149–0.3917 |
| 140 | fixed catalogue | 1.1769e-6–1.9726e-6 | 1 | 1 |
| 140 | linear | 0.9394e-6–1.7684e-6 | 0.7982–0.8965 | 0.6371–0.8036 |
| 140 | quadratic | 0.5480e-6–0.8925e-6 | 0.3918–0.4657 | 0.1535–0.2168 |

Under this conditional metric, both fits reduce the F–Mplus distance in every tested case; quadratic fitting reduces it further than linear fitting. These are not statements about the best treatment when an uncertain foreground is present. Changes in output amplitude and changes in propagated noise partly compensate, so a raw power-difference ratio cannot replace the conditional distance ratio.

For example, at 130 MHz, the median over four groups of the lowest-delay absolute difference spans 1.6921e-7–7.1413e-7 Jy² for fixed-catalogue processing, versus 7.2437e-9–1.1549e-8 Jy² after quadratic fitting. Yet the corresponding median absolute **normalized** contrast, `abs(P_F-P_Mplus)/(P_F+P_Mplus)`, spans 0.2948–0.4918 and 0.2586–0.4822, respectively. These are distinct diagnostics: the normalization itself changes under fitting. Detailed powers, sums, differences, normalized contrasts and diagonal error scales are in `per_output.csv`. No reporting floor is needed for the distance calculation; its output differences enter continuously.

## Covariance and inversion checks

All 36 frequency/seed/treatment/residual cases are positive definite and have rank 24 at relative eigenvalue cutoff 1e-12. Approximate primary covariance condition numbers (the seed dependence is negligible at the quoted precision) are:

| MHz | fixed catalogue | linear | quadratic |
| --- | --- | --- | --- |
| 130 | 8.38 | 9.54 | 666.13 |
| 140 | 8.41 | 9.45 | 708.56 |

The eigensolve and independent Cholesky solve agree in D to maximum relative difference 4.44e-16. Relative pseudoinverse cutoffs from 1e-14 through 1e-3 retain all 24 modes and leave D unchanged. The deliberately aggressive 1e-2 cutoff drops four modes in quadratic cases (rank 20), so those values are a truncation stress test, not the adopted result. `pseudoinverse_stability.csv` retains every cutoff and distance, rather than silently choosing a pseudoinverse threshold.

For every input archive the stored results hash is checked. Powers are reconstructed from F/Mplus visibilities using the returned Cleaner and estimator code; F covariances are independently reconstructed from the local noise variance and transformed reference sky, including the residual where applicable. Maximum relative covariance reconstruction error is 4.27e-16. This validates array use and algebra under the stipulated operator; it does not validate missing physical covariance.

## Secondary residual case

The CSVs also retain `dipole_truth_gaussian_catalogue` (residual index 1). Here `C_F` includes the same deterministic dipole-minus-Gaussian foreground residual used in the archive. `distance` continues to weight the signal-power difference alone. `observed_mean_distance` instead uses `(P_F + b_F) - (P_Mplus + b_Mplus)`, where the archived bias includes foreground residual power and signal–residual cross terms. A common residual does not imply identical model bias because the signal–residual cross term depends on the model. Neither distance marginalizes over uncertain foregrounds or estimates a foreground-error ensemble. Thus this case must not be advertised as robustness to foreground uncertainty.

## Reproduce and files

From the repository root:

```sh
python output/mnras_revision_20261104/analysis_visibility/analyze_visibility.py
```

Only NumPy and the Python standard library are required. The script reads directly from the zip and temporarily loads only its returned `estimator.py` and `cleaner.py` for reconstruction. All destinations can be changed with `--archive` and `--out`.

- `distances.csv`: 36 rows of distances, absolute/normalized contrasts, reference covariance diagnostics and treatment ratios.
- `per_output.csv`: 864 rows with all 24 retained output powers and contrasts per case.
- `pseudoinverse_stability.csv`: 288 rows covering eight eigenvalue cutoffs per case.
- `validation.csv`: archive hash, power/covariance reconstruction and positivity checks.
- `provenance.json`: archive SHA-256, software version and precise interpretation limits.

## Suggested manuscript text

As an additional conditional diagnostic, we evaluated F–Mplus lightcone differences through the same grouped-delay estimator and its propagated F-reference covariance, retaining the first six delay bins in each of four baseline groups. For dipole beams with matched-catalogue subtraction, quadratic fitting leaves 0.339–0.626 of the fixed-catalogue distance at 130 MHz and 0.392–0.466 at 140 MHz across three seeds; linear fitting leaves 0.638–0.951 and 0.798–0.896, respectively. The retained 24-dimensional covariances are full rank, and the results are stable to relative eigenvalue cutoffs through 1e-3. These very small distances use the stipulated short-split thermal and fixed-sky mixed covariance; they are not 1000-hour forecasts or model-exclusion significances. Their ratios differ from raw output-contrast ratios because fitting changes both the model difference and covariance. The frequency windows and omitted sky/foreground covariance prevent an observational frequency ranking.
