#!/usr/bin/env python3
"""Reforecast matched grids using saved, audited public-convention instrument weights.

For fixed geometry, Gaussian sample error = alpha(kperp,kpar)*Delta2_model.
Validate alpha across ALL saved models and total=thermal+sample before reuse.
This preserves 21cmSense 3.0.0 public error ordering, not per-cell joint weighting.
Repeated radial entries are preserved; their UV provenance must already be audited.
"""
import argparse,csv,hashlib,json,runpy
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
C=runpy.run_path(str(HERE/'signal_core.py'))
H=runpy.run_path(str(HERE/'compare_histories.py'))
CASES=['no_wedge','beam_wedge','horizon_buffer']
EDGES=np.array([.05,.075,.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35])

def check(ok,message):
    if not ok:raise ValueError(message)

def close(x,y,label,rtol=1e-6,atol=1e-10):
    check(np.shape(x)==np.shape(y) and np.allclose(x,y,rtol=rtol,atol=atol),label)

def interp(k,x,y):
    check(np.isfinite(x).all() and np.isfinite(y).all() and np.all(y>0),'Need positive finite spectra; no artificial floor')
    check(np.all(np.diff(x)>0),'Unordered spectrum')
    check(np.all((k>=x[0])&(k<=x[-1])),'Extrapolation forbidden')
    return np.exp(np.interp(np.log(k),np.log(x),np.log(y)))

def snr(s,e,mask):return float(np.linalg.norm(s[mask]/e[mask]))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def saved_weights(folder):
    meta=json.loads((folder/'forecast_metadata.json').read_text())
    check(meta.get('thermal_only') is False,'Need full thermal + sample forecast')
    sm=meta['signal_metadata']
    check(sm['units']==dict(k='comoving Mpc^-1',P='mK^2 Mpc^3',Delta2='mK^2'),'Signal units mismatch')
    check(sm['coordinate_space']=='real' and sm['mean_removed'],'Need real-space, mean-subtracted spectra')
    check(abs(sm['frequency_MHz']-130)<1e-4,'Not a 130 MHz forecast')
    with np.load(folder/'forecast_native.npz',allow_pickle=False) as z:d={k:z[k] for k in z.files}
    k=d['k_Mpc_inv'];t=d['thermal_sigma_mK2']
    check(k.ndim==1 and len(k)>0,'Invalid native samples')
    close(k,np.hypot(d['kperp_Mpc_inv'],d['kpar_Mpc_inv']),'Native k geometry mismatch',rtol=1e-10)
    check(t.shape==k.shape and np.all(np.isfinite(t)&(t>0)),'Invalid thermal errors')
    for case in CASES:
        check(d[case].dtype==bool and d[case].shape==k.shape,'Expected Boolean native masks')
    check(np.all(~d['horizon_buffer']|d['beam_wedge']) and np.all(~d['beam_wedge']|d['no_wedge']),'Masks not nested')
    models=sorted(v[:-len('_signal_mK2')] for v in d if v.endswith('_signal_mK2'))
    check('F' in models,'Saved F signal needed')
    coeff=[]
    for m in models:
        s=d[m+'_signal_mK2'];ss=d[m+'_sample_sigma_mK2'];total=d[m+'_total_sigma_mK2']
        check(s.shape==ss.shape==total.shape==k.shape,'Native array shape mismatch')
        check(np.all(np.isfinite(s)&(s>0)) and np.all(np.isfinite(ss)&(ss>0)),'Non-positive signal/sample error')
        close(total,t+ss,'Saved total is not thermal+sample public convention')
        coeff.append(ss/s)
    alpha=coeff[models.index('F')]
    for m,a in zip(models,coeff):close(a,alpha,'Sample error / signal is not geometry-only: '+m)
    # Reproduce the published full-range SNR, catching stale or mixed inputs.
    with (folder/'detection_snr.csv').open() as f:rows=list(csv.DictReader(f))
    for m in models:
        for case in CASES:
            r=[r for r in rows if r['model']==m and r['case']==case]
            check(len(r)==1,'Missing or repeated saved detection row')
            close(snr(d[m+'_signal_mK2'],d[m+'_total_sigma_mK2'],d[case]),float(r[0]['forecast_snr']),
                  'Cannot reproduce saved SNR '+m+'/'+case)
    print('PASS: saved SNRs, public covariance identity, and model-independent sample weights',flush=True)
    return meta,d,alpha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--low-root',type=Path,required=True);p.add_argument('--high-root',type=Path,required=True)
    p.add_argument('--forecast',type=Path,default=Path('forecast_130MHz_full'))
    p.add_argument('--models',nargs='+',choices=C['MODELS'],default=['F','Mminus','Iminus','Mplus'])
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();models=a.models
    check(len(set(models))==len(models) and len(models)>=2 and 'F' in models,'Need unique models including F and at least one alternative')
    if a.out.exists() and any(a.out.iterdir()):raise FileExistsError('Choose new output directory')
    meta,d,alpha=saved_weights(a.forecast);k=d['k_Mpc_inv'];t=d['thermal_sigma_mK2']
    krange=meta['k_range_Mpc_inv'];check(krange[0]<=.1 and krange[1]>=.3,'Saved forecast must cover 0.1--0.3')
    use=(k>=.1)&(k<=.3);check(use.any(),'No usable modes')
    # Identical native rows/masks, with duplicate-coordinate groups kept intact.
    k=k[use];t=t[use];alpha=alpha[use];masks={c:d[c][use] for c in CASES}
    spec={64:{},128:{}};provenance={};rows=[];configs={}
    ref_common=meta['signal_metadata']['common_configuration']
    requested_astro=meta['signal_metadata']['requested_astro_inputs']
    common=None
    for m in models:
        info_pair={}
        for n,root in [(64,a.low_root),(128,a.high_root)]:
            check((root/m/'COMPLETE.json').exists(),'Incomplete simulation '+str(root/m))
            cube,L,cfg,info=C['load_model'](root,m,130.)
            configs[n,m]=cfg;info_pair[str(n)]=info
            check(cfg['user_params']['HII_DIM']==n and cfg['user_params']['DIM']==384,'Grid mismatch')
            for key in ['cosmo_params','flag_options','random_seed']:
                check(cfg[key]==ref_common[key],'Original forecast configuration mismatch: '+key)
            check(cfg['astro_params_log10_inputs']==requested_astro[m],'Original model physics mismatch: '+m)
            for key,value in ref_common['user_params'].items():
                if key not in ('HII_DIM','DIM','N_THREADS'):
                    check(cfg['user_params'].get(key)==value,'Original user-parameter mismatch '+key)
            base={key:cfg[key] for key in ['cosmo_params','flag_options','random_seed','redshifts']}
            check(common is None or base==common,'Models have inconsistent common histories/configuration')
            common=base
            check(EDGES[-1]<np.pi*n/L,'Bins reach Nyquist')
            kk,pp,mean,vr,vf=C['fft_power'](cube,L)
            km,pm,de,count,_=C['spherical'](kk,pp,EDGES);pop=count>0
            spec[n][m]=interp(k,km[pop],de[pop])
            info.update(mean_Tb_mK=mean,variance_mK2=vr,parseval_fft_variance_mK2=vf)
            for i in range(len(count)):
                rows.append(dict(model=m,N=n,k_low_Mpc_inv=EDGES[i],k_high_Mpc_inv=EDGES[i+1],
                    k_mean_Mpc_inv=km[i],fft_cells=int(count[i]),P_mK2_Mpc3=pm[i],Delta2_mK2=de[i]))
        H['compare_configs'](configs[64,m],configs[128,m],m)
        provenance[m]=info_pair
        print('PASS: matched spectra and model configuration',m,flush=True)
    a.out.mkdir(parents=True,exist_ok=True)
    H['write_csv'](a.out/'matched_signal_spectra.csv',rows)
    np.savetxt(a.out/'k_bin_edges_Mpc_inv.txt',EDGES,fmt='%.17g')
    errors={n:{m:t+alpha*spec[n][m] for m in models} for n in spec}
    det=[];pairs=[];shifts=[]
    for cut in [.2,.25,.3]:
        for case in CASES:
            mask=masks[case]&(k<=cut)
            for m in models:
                lo=snr(spec[64][m],errors[64][m],mask);hi=snr(spec[128][m],errors[128][m],mask)
                det.append(dict(case=case,kmax_Mpc_inv=cut,model=m,snr_64=lo,snr_128=hi,
                    ratio_128_to_64=hi/lo if lo>0 else '',native_entries=int(mask.sum())))
                shifts.append(dict(case=case,kmax_Mpc_inv=cut,model=m,
                    resolution_shift_distance_using_64_covariance=snr(spec[128][m]-spec[64][m],errors[64][m],mask)))
            for i,m in enumerate(models):
                for other in models[i+1:]:
                    row=dict(case=case,kmax_Mpc_inv=cut,model_a=m,model_b=other)
                    for n in [64,128]:
                        diff=spec[n][m]-spec[n][other]
                        for name,ref in [('fixed_F','F'),('if_a_true',m),('if_b_true',other)]:
                            row[f'distance_{n}_{name}']=snr(diff,errors[n][ref],mask)
                        # Additional diagnostic separates template shift from covariance change.
                        row[f'distance_{n}_common_64_F_covariance']=snr(diff,errors[64]['F'],mask)
                    lo=row['distance_64_fixed_F'];hi=row['distance_128_fixed_F']
                    row['ratio_fixed_F_128_to_64']=hi/lo if lo>0 else ''
                    pairs.append(row)
    H['write_csv'](a.out/'resolution_detection_comparison.csv',det)
    H['write_csv'](a.out/'resolution_pairwise_comparison.csv',pairs)
    H['write_csv'](a.out/'resolution_template_shifts.csv',shifts)
    for n in [64,128]:
        np.savez_compressed(a.out/f'forecast_native_N{n}.npz',k_Mpc_inv=k,
            kperp_Mpc_inv=d['kperp_Mpc_inv'][use],kpar_Mpc_inv=d['kpar_Mpc_inv'][use],
            thermal_sigma_mK2=t,sample_coefficient=alpha,**masks,
            **{m+'_signal_mK2':s for m,s in spec[n].items()},
            **{m+'_sample_sigma_mK2':alpha*s for m,s in spec[n].items()},
            **{m+'_total_sigma_mK2':s for m,s in errors[n].items()})
    for table,name,key,label in [(det,'resolution_detection','snr','Template amplitude SNR'),
                                 (pairs,'resolution_pairs','distance','Fixed-F template distance')]:
        fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        for ax,case in zip(axes,CASES):
            rr=[r for r in table if r['case']==case and r['kmax_Mpc_inv']==.3]
            xx=np.arange(len(rr))
            for off,n in [(-.18,64),(.18,128)]:
                col=f'snr_{n}' if key=='snr' else f'distance_{n}_fixed_F'
                ax.bar(xx+off,[r[col] for r in rr],width=.36,label=f'{n}³')
            labels=[r['model'] if key=='snr' else r['model_a']+' / '+r['model_b'] for r in rr]
            ax.set_xticks(xx,labels,rotation=65,ha='right');ax.set_title(case);ax.set_ylabel(label);ax.grid(axis='y',alpha=.2);ax.legend()
        fig.suptitle('130 MHz | matched DIM=384 | 0.1–0.3 Mpc⁻¹ | public covariance convention')
        fig.savefig(a.out/f'{name}.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    for m in models:
        rr=[r for r in rows if r['N']==64 and r['model']==m and r['fft_cells']>0]
        rh=[r for r in rows if r['N']==128 and r['model']==m and r['fft_cells']>0]
        ax.plot([r['k_mean_Mpc_inv'] for r in rr],[h['Delta2_mK2']/l['Delta2_mK2'] for l,h in zip(rr,rh)],'o-',label=m)
    ax.axhline(1,color='black');ax.axhspan(.9,1.1,alpha=.15,color='gray');ax.set(xlim=(.1,.3),xlabel='k [Mpc⁻¹]',ylabel='128³ / 64³ power');ax.legend();ax.grid(alpha=.2)
    fig.savefig(a.out/'resolution_power_ratios.png',dpi=180);plt.close(fig)
    report=dict(status='PASS algebraic weight reuse checks',models=models,frequency_MHz=130,
        covariance='21cmSense 3.0.0 public: total_sigma = thermal_sigma + alpha * model_Delta2',
        reference_forecast=str(a.forecast.resolve()),
        reference_sha256={f:digest(a.forecast/f) for f in ['forecast_native.npz','forecast_metadata.json','detection_snr.csv']},
        source_cubes=provenance,k_range_Mpc_inv=[.1,.3],edges_Mpc_inv=EDGES.tolist(),
        notes=['No SKA geometry or thermal errors are recalculated or rescaled.',
               'Sample coefficient checked against every positive saved model, not inferred from Fourier counts.',
               'Preserves saved duplicate radial entries and Gaussian diagonal covariance assumptions; prior UV provenance audit remains required.',
               'Not the per-cell joint thermal+sample covariance alternative.',
               'Fixed-F covariance changes with resolution; common-64-F diagnostic also exported.',
               'Distances are template separations, not Bayes factors or calibrated rejection significances.',
               'Both resolutions use DIM=384; not a direct correction of original DIM=192 results.',
               'Identical bin edges and interpolation for both grids. No signal extrapolation.',
               'Two resolutions and one seed do not establish convergence or realization variance.'])
    (a.out/'reforecast_checks.json').write_text(json.dumps(report,indent=2))
    for r in det:
        if r['case']=='horizon_buffer' and r['kmax_Mpc_inv']==.3:
            print(f"{r['model']:7s} horizon SNR: {r['snr_64']:.5g} -> {r['snr_128']:.5g}")
    print('Finished:',a.out.resolve())

if __name__=='__main__':main()
