#!/usr/bin/env python3
"""Same-operator conditional F--Mplus diagnostic from the returned archive.

Run from the repository root. Requires Python 3 and NumPy only.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import zipfile

import numpy as np


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('matched_response_review/matched_response_review.zip'))
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary, outputs, stability, checks = [], [], [], []
    keep = np.arange(28).reshape(4, 7)[:, :6].ravel()
    cuts = [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2]
    with zipfile.ZipFile(args.archive) as archive, tempfile.TemporaryDirectory() as tmp:
        # Reconstruct powers/covariance from archived visibilities, separately
        # from the new eigensolve; use exactly the returned estimator code.
        for name in ('estimator', 'cleaner'):
            Path(tmp, name + '.py').write_bytes(archive.read('code/' + name + '.py'))
        sys.path.insert(0, tmp)
        import estimator as est
        from cleaner import Cleaner
        for nu in (130, 140):
            noise = np.load(io.BytesIO(archive.read(f'{nu}MHz/noise/noise.npz')))
            var = noise['complex_variance_I_Jy2']
            ops = [Cleaner(var, noise['frequencies_hz'], degree) for degree in (-1, 1, 2)]
            for seed in (12345, 23456, 34567):
                prefix = f'{nu}MHz/seed{seed}/results/'
                raw = archive.read(prefix + 'results.npz')
                complete = json.loads(archive.read(prefix + 'complete.json'))
                assert hashlib.sha256(raw).hexdigest() == complete['results_sha256']
                data = np.load(io.BytesIO(raw))
                config = json.loads(archive.read(prefix + 'configuration.json'))
                assert config['beam_order'][0] == 'dipole'
                fi = config['field_order'].index('F_lightcone')
                mi = config['field_order'].index('Mplus_lightcone')
                for j, pipeline in enumerate(config['pipeline_order']):
                    op = ops[j]
                    pf, pm = data['signal_power_Jy2'][0, [fi, mi], j]
                    for field, power in ((fi, pf), (mi, pm)):
                        rec = op.power(data['visibilities_Jy'][0, field])
                        np.testing.assert_allclose(rec, power, rtol=1e-11, atol=1e-24)
                    sf = op.apply(data['visibilities_Jy'][0, fi])
                    residual = op.apply(data['foreground_Jy'][0] - data['foreground_Jy'][2])
                    for case, label in enumerate(('matched_catalogue', 'dipole_truth_gaussian_catalogue')):
                        r = np.zeros_like(residual) if case == 0 else residual
                        thermal, mixed = est.theory(var, op.L, sf + r)
                        archived = data['conditional_covariance_Jy4'][fi, j, case]
                        cov_error = np.linalg.norm(thermal + mixed - archived) / np.linalg.norm(archived)
                        np.testing.assert_allclose(thermal + mixed, archived, rtol=1e-10, atol=1e-20)
                        c = archived[np.ix_(keep, keep)]
                        assert np.all(np.isfinite(c))
                        np.testing.assert_allclose(c, c.T, rtol=1e-12, atol=1e-20)
                        c = (c + c.T) / 2
                        eigen, vectors = np.linalg.eigh(c)
                        assert eigen[0] > 0
                        delta = (pf - pm)[keep]
                        total = (pf + pm)[keep]
                        bias_difference = (data['bias_Jy2'][fi, j, case] - data['bias_Jy2'][mi, j, case])[keep]
                        observed_delta = delta + bias_difference
                        coordinates = vectors.T @ delta
                        d2 = np.sum(coordinates ** 2 / eigen)
                        chol = np.linalg.cholesky(c)
                        independent_d2 = np.sum(np.linalg.solve(chol, delta) ** 2)
                        np.testing.assert_allclose(independent_d2, d2, rtol=1e-11, atol=0)
                        context = dict(frequency_MHz=nu, seed=seed, beam='dipole', sky='lightcone', pipeline=pipeline, covariance_case=label)
                        row = dict(**context, n_outputs=len(keep), reference_field='F_lightcone', distance=np.sqrt(d2), distance_squared=d2,
                                   observed_mean_distance=np.sqrt(np.sum((vectors.T @ observed_delta)**2/eigen)),
                                   rank_rtol_1e_minus12=int(np.sum(eigen > eigen[-1]*1e-12)), eigen_min_Jy4=eigen[0], eigen_max_Jy4=eigen[-1], condition_2=eigen[-1]/eigen[0],
                                   l2_absolute_difference_Jy2=np.linalg.norm(delta), median_absolute_difference_Jy2=np.median(abs(delta)),
                                   median_absolute_normalized_contrast=np.median(abs(delta/total)),
                                   lowest_delay_median_absolute_difference_Jy2=np.median(abs(delta.reshape(4,6)[:,0])),
                                   lowest_delay_median_absolute_normalized_contrast=np.median(abs((delta/total).reshape(4,6)[:,0])),
                                   covariance_reconstruction_relative_error=cov_error, cholesky_distance_relative_error=abs(np.sqrt(independent_d2/d2)-1))
                        summary.append(row)
                        for cutoff in cuts:
                            mask = eigen > eigen[-1] * cutoff
                            value = np.sqrt(np.sum(coordinates[mask]**2/eigen[mask]))
                            stability.append(dict(**context, relative_eigen_cutoff=cutoff, effective_rank=int(mask.sum()), distance=value, relative_change_from_full_rank=value/np.sqrt(d2)-1))
                        for k, index in enumerate(keep):
                            outputs.append(dict(**context, group=index//7, delay_bin=index%7, delay_lo_us=config['delay_edges_us'][index%7], delay_hi_us=config['delay_edges_us'][index%7+1],
                                                F_power_Jy2=pf[index], Mplus_power_Jy2=pm[index], signed_difference_Jy2=delta[k], absolute_difference_Jy2=abs(delta[k]),
                                                power_sum_Jy2=total[k], signed_normalized_contrast=delta[k]/total[k], conditional_standard_deviation_Jy2=np.sqrt(c[k,k]),
                                                difference_over_diagonal_std=delta[k]/np.sqrt(c[k,k]), foreground_bias_difference_Jy2=bias_difference[k]))
                        checks.append(dict(**context, hash_verified=True, powers_reconstructed=True, covariance_reconstructed=True, positive_definite=True, covariance_relative_error=cov_error))
    # Comparisons concern treatments within each frequency and seed only.
    for row in summary:
        fixed = next(x for x in summary if x['frequency_MHz']==row['frequency_MHz'] and x['seed']==row['seed'] and x['covariance_case']==row['covariance_case'] and x['pipeline']=='fixed_catalogue')
        row['distance_ratio_to_fixed'] = row['distance']/fixed['distance']
        row['squared_distance_ratio_to_fixed'] = row['distance_squared']/fixed['distance_squared']
    for name, rows in [('distances', summary), ('per_output', outputs), ('pseudoinverse_stability', stability), ('validation', checks)]:
        write_csv(args.out / (name + '.csv'), rows)
    provenance = dict(archive=str(args.archive), archive_sha256=hashlib.sha256(args.archive.read_bytes()).hexdigest(), numpy_version=np.__version__,
                      primary='dipole lightcones; matched_catalogue covariance (no foreground residual); F reference sky, separately propagated per pipeline',
                      secondary='dipole truth minus Gaussian catalogue residual; same F reference sky; observed_mean_distance also includes differential signal-residual bias',
                      metric='sqrt((P_F-P_Mplus)^T C_F^{-1} (P_F-P_Mplus)); original 28-dimensional covariance marginalized to retained 24 outputs before inversion',
                      support='first six delay bins in each of four baseline groups; excludes 4--8.01 us final bin',
                      scope='Fixed-sky, 10-second-per-split thermal plus sky-noise mixed covariance. No cosmological sample variance, calibration uncertainty, residual ensemble or 1000-hour rescaling. Not a rejection significance, information fraction, or cross-frequency observing ranking.',
                      validation='archive hashes; visibility-power and conditional-covariance reconstruction; eigensolve versus Cholesky; cutoff sweep')
    (args.out / 'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    for row in summary:
        if row['covariance_case']=='matched_catalogue':
            print(row['frequency_MHz'], row['seed'], row['pipeline'], f"D={row['distance']:.9g}", f"D/fixed={row['distance_ratio_to_fixed']:.6f}", f"condition={row['condition_2']:.3f}")


if __name__ == '__main__':
    main()
