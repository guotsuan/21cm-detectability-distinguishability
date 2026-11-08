#!/usr/bin/env python3
"""Per-band forecasts at both resolutions; no summation across bands."""
import argparse,csv,json,runpy
from pathlib import Path
import numpy as np
G=runpy.run_path(str(Path(__file__).resolve().parent/'common.py'));globals().update({k:v for k,v in G.items() if not k.startswith('__')})
def evaluate(nu,native,signals):
    k=native['k_Mpc_inv'];t=native['thermal_sigma_mK2'];alpha=native['sample_coefficient'];err={n:{m:t+alpha*s for m,s in ss.items()} for n,ss in signals.items()}
    det=[];pairs=[];shifts=[];scale=[]
    for cut in [.2,.25,.3]:
        for case in CASES:
            mask=native[case]&(k>=.1)&(k<=cut)
            for m in MODELS:
                lo=metric(signals[96][m],err[96][m],mask);hi=metric(signals[128][m],err[128][m],mask)
                det.append(dict(frequency_MHz=nu,case=case,kmax_Mpc_inv=cut,model=m,N_low=96,N_high=128,snr_low=lo,snr_high=hi,ratio=hi/lo if lo is not None and lo>0 else None,native_entries=int(mask.sum()),thermal_only_snr_high=metric(signals[128][m],t,mask),status='no_retained_modes' if not mask.any() else ('zero_signal_both' if lo==hi==0 else 'ok')))
                shifts.append(dict(frequency_MHz=nu,case=case,kmax_Mpc_inv=cut,model=m,distance_using_low_covariance=metric(signals[128][m]-signals[96][m],err[96][m],mask)))
            for i,a in enumerate(MODELS):
                for b in MODELS[i+1:]:
                    r=dict(frequency_MHz=nu,case=case,kmax_Mpc_inv=cut,model_a=a,model_b=b,N_low=96,N_high=128)
                    for n,label in [(96,'low'),(128,'high')]:
                        delta=signals[n][a]-signals[n][b]
                        for lab,ref in [('fixed_F','F'),('if_a_true',a),('if_b_true',b)]:r[f'distance_{label}_{lab}']=metric(delta,err[n][ref],mask)
                        r[f'distance_{label}_common_low_F_covariance']=metric(delta,err[96]['F'],mask)
                    r['status']='no_retained_modes' if not mask.any() else 'ok';pairs.append(r)
    for case in CASES:
        for left,right in zip([.1,.125,.15,.175,.2,.225,.25,.275],[.125,.15,.175,.2,.225,.25,.275,.3]):
            mask=native[case]&(k>=left)&((k<=right) if right==.3 else (k<right))
            for m in MODELS:
                v=metric(signals[128][m],err[128][m],mask)
                scale.append(dict(frequency_MHz=nu,case=case,model=m,k_low_Mpc_inv=left,k_high_Mpc_inv=right,detection_snr2_high=0. if v is None else v*v,native_entries=int(mask.sum())))
    return det,pairs,shifts,scale,err
def regression(rows,path,key_fields,metrics):
    with Path(path).open() as f:old=list(csv.DictReader(f))
    check(len(rows)==len(old),'130 MHz comparison row count differs: '+str(path))
    def key(r):return tuple(str(float(r[k])) if k=='kmax_Mpc_inv' else str(r[k]) for k in key_fields)
    mapped={key(r):r for r in old};check(len(mapped)==len(old),'Duplicate comparison keys')
    largest=0.
    for r in rows:
        o=mapped[key(r)]
        for field in metrics:
            x=float(o[field]);y=r[field];check(y is not None and np.isclose(y,x,rtol=1e-6,atol=1e-9),'130 MHz result changed: '+str(key(r))+'/'+field)
            largest=max(largest,abs(y-x)/max(abs(x),1e-12))
    return largest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectra',type=Path,default=Path('frequency_analysis/spectra'));p.add_argument('--noise',type=Path,default=Path('frequency_analysis/noise'))
    p.add_argument('--comparison',type=Path,default=Path('grid_comparison_96_128_all7'))
    p.add_argument('--frequencies',type=float,nargs='+',default=[130]);p.add_argument('--out',type=Path,default=Path('frequency_analysis/forecasts'));p.add_argument('--resume',action='store_true')
    a=p.parse_args();check(len(set(a.frequencies))==len(a.frequencies),'Duplicate centres');check((a.spectra/'COMPLETE.json').exists(),'Run extraction to completion first')
    job=dict(spectra_job=sha(a.spectra/'job.json'),noise_job=sha(a.noise/'job.json'),code=code_hash(),comparison_hashes={f:sha(a.comparison/f) for f in ['detection_comparison.csv','pairwise_comparison.csv']})
    a.out.mkdir(parents=True,exist_ok=True);jp=a.out/'job.json'
    if jp.exists():check(a.resume and json.loads(jp.read_text())==job,'Changed inputs or --resume missing; choose new output')
    else:check(not any(a.out.iterdir()),'Nonempty output');json_write(jp,job)
    for nu in sorted(a.frequencies,key=lambda v:(v!=130,v)):
        folder=a.out/tag(nu)
        if nu!=130:check((a.out/'130MHz/COMPLETE.json').exists(),'First reproduce 130 MHz results')
        nd=a.noise/tag(nu);ni=json.loads((nd/'COMPLETE.json').read_text());check(sha(nd/'native.npz')==ni['native_sha256'],'Noise checksum mismatch')
        check(abs(ni['frequency_MHz']-nu)<1e-7,'Noise frequency mismatch')
        inputs={'noise':sha(nd/'native.npz')}
        for n in [96,128]:
            for m in MODELS:
                sp=signal_file(a.spectra,nu,m,n);inputs[f'{m}_{n}']=sha(sp)
                si=json.loads(sp.with_suffix('.json').read_text());check(inputs[f'{m}_{n}']==si['spectrum_sha256'],'Spectrum corrupt')
        if (folder/'COMPLETE.json').exists():
            check(a.resume,'Use --resume');done=json.loads((folder/'COMPLETE.json').read_text());check(done['inputs']==inputs,'Inputs changed')
            for f,v in done['outputs'].items():check(sha(folder/f)==v,'Forecast checkpoint corrupt')
            print('Verified resume:',nu,flush=True);continue
        native=columns(nd/'native.npz');use=(native['k_Mpc_inv']>=.1)&(native['k_Mpc_inv']<=.3);native={key:v[use] for key,v in native.items()}
        signal={n:{m:signal_at(a.spectra,nu,m,n,native['k_Mpc_inv']) for m in MODELS} for n in [96,128]}
        det,pairs,shifts,scale,err=evaluate(nu,native,signal)
        reg=None
        if nu==130:
            reg=dict(status='PASS',detection_max_relative_difference=regression(det,a.comparison/'detection_comparison.csv',['case','kmax_Mpc_inv','model'],['snr_low','snr_high']),pair_max_relative_difference=regression(pairs,a.comparison/'pairwise_comparison.csv',['case','kmax_Mpc_inv','model_a','model_b'],[f'distance_{side}_{c}' for side in ['low','high'] for c in ['fixed_F','if_a_true','if_b_true','common_low_F_covariance']]))
        folder.mkdir(parents=True,exist_ok=True)
        outputs={}
        for name,rows in [('detection.csv',det),('pairs.csv',pairs),('template_shifts.csv',shifts),('scale_information.csv',scale)]:csv_write(folder/name,rows);outputs[name]=sha(folder/name)
        for n in [96,128]:
            file=folder/f'native_N{n}.npz';np.savez_compressed(file,**native,**{m+'_signal_mK2':v for m,v in signal[n].items()},**{m+'_total_sigma_mK2':v for m,v in err[n].items()});outputs[file.name]=sha(file)
        json_write(folder/'COMPLETE.json',dict(frequency_MHz=nu,inputs=inputs,outputs=outputs,regression_130=reg,covariance='public 3.0.0 thermal_sigma + alpha * signal; diagonal',cross_frequency_combination=False));print('Forecasted',nu,'MHz',flush=True)
    # Rebuild summaries from only requested bands, checking each completed checkpoint.
    for name in ['detection','pairs','template_shifts','scale_information']:
        rr=[]
        for nu in sorted(a.frequencies):
            with (a.out/tag(nu)/(name+'.csv')).open() as f:rr.extend(list(csv.DictReader(f)))
        csv_write(a.out/(name+'_all_bands.csv'),rr)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with (a.out/'detection_all_bands.csv').open() as f:rr=list(csv.DictReader(f))
    fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for ax,c in zip(axs,CASES):
        for m in MODELS:
            r=[r for r in rr if r['case']==c and r['model']==m and float(r['kmax_Mpc_inv'])==.3]
            ax.plot([float(x['frequency_MHz']) for x in r],[float(x['snr_high']) if x['snr_high'] else np.nan for x in r],'.-',label=m)
        ax.set(title=c,xlabel='Frequency [MHz]',ylabel='Conditional SNR');ax.grid(alpha=.2)
    axs[0].legend(fontsize=8);fig.savefig(a.out/'detection_vs_frequency.png',dpi=160);plt.close(fig)
    json_write(a.out/'summary.json',dict(frequencies=sorted(a.frequencies),models=MODELS,notes=['No joint-band SNR or significance is calculated.','Zero signal yields zero SNR if modes exist; no retained modes yields blank, not zero.','Full k range retained for diagnostics; consult resolution_map before scientific claims.','Binary masks and coeval real-space signal; no measured response or foreground leakage.']))
if __name__=='__main__':main()
