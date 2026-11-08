#!/usr/bin/env python3
"""Extract all saved coeval spectra without rerunning 21cmFAST. Serial, checkpointed."""
import argparse,csv,json,runpy
from pathlib import Path
import numpy as np
G=runpy.run_path(str(Path(__file__).resolve().parent/'common.py'));globals().update({k:v for k,v in G.items() if not k.startswith('__')})
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--low-root',type=Path,default=Path('matched_130MHz_resolution/res96d384_seed12345_real'))
    p.add_argument('--high-root',type=Path,default=Path('matched_130MHz_resolution/res128d384_seed12345_real'))
    p.add_argument('--forecast',type=Path,default=Path('forecast_130MHz_full'))
    p.add_argument('--frequencies',type=float,nargs='+',default=list(range(50,201,5)))
    p.add_argument('--out',type=Path,default=Path('frequency_analysis/spectra'))
    p.add_argument('--resume',action='store_true');p.add_argument('--preflight',action='store_true')
    a=p.parse_args();ref=ref_meta(a.forecast);roots={96:a.low_root,128:a.high_root};configs={};paths={};hist={}
    check(len(set(a.frequencies))==len(a.frequencies),'Repeated frequencies')
    common=None
    for m in MODELS:
        for n,root in roots.items():
            folder=root/m;check((folder/'COMPLETE.json').exists(),'Missing COMPLETE: '+str(folder))
            cfg=json.loads((folder/'requested_config.json').read_text());config_checks(cfg,m,n,ref);configs[n,m]=cfg
            cc={k:cfg[k] for k in ['cosmo_params','flag_options','random_seed','redshifts']}
            check(common is None or cc==common,'Common schedules/physics differ');common=cc
            with (folder/'history.csv').open() as f:hh=list(csv.DictReader(f))
            check(len(hh)==len(cfg['redshifts']),'Incomplete history')
            for nu in a.frequencies:
                z=NU21/nu-1;check(min(abs(z-float(v)) for v in cfg['redshifts'])<1e-6,'Frequency absent in config')
                found=[f for f in folder.glob('coeval_z*.npz') if abs(float(f.stem.split('_z')[1])-z)<.00051]
                check(len(found)==1,f'Need exactly one cube: {m} N{n} {nu} MHz');paths[n,m,nu]=found[0]
                rr=[r for r in hh if abs(float(r['frequency_MHz'])-nu)<1e-4];check(len(rr)==1,'Missing/duplicate history row');hist[n,m,nu]=rr[0]
        matched(configs[96,m],configs[128,m])
    print(f'PASS preflight: {len(paths)} cube paths; all seven models; fixed DIM=384.',flush=True)
    if a.preflight:return
    job=dict(roots={str(n):str(r.resolve()) for n,r in roots.items()},frequencies=sorted(a.frequencies),configs={f'{n}_{m}':v for (n,m),v in configs.items()},edges=EDGES.tolist(),reference_sha256=sha(a.forecast/'forecast_metadata.json'),code=code_hash())
    a.out.mkdir(parents=True,exist_ok=True);manifest=a.out/'job.json'
    if manifest.exists():check(a.resume and json.loads(manifest.read_text())==job,'Output config changed or --resume missing; choose new --out')
    else:check(not any(a.out.iterdir()),'Nonempty output without manifest');json_write(manifest,job)
    allrows=[];stats=[];compare=[]
    for nu in sorted(a.frequencies):
        for m in MODELS:
            spectra={}
            for n,root in roots.items():
                target=signal_file(a.out,nu,m,n);jp=target.with_suffix('.json');path=paths[n,m,nu];digest=sha(path)
                if a.resume and target.exists() and jp.exists():
                    info=json.loads(jp.read_text());check(info['cube_sha256']==digest,'Source cube changed; choose new output')
                    check(info['spectrum_sha256']==sha(target),'Checkpoint corrupt');z=columns(target)
                else:
                    cube,L,cfg,info=CORE['load_model'](root,m,nu)
                    check(EDGES[-1]<np.pi*n/L,'Nyquist too low')
                    kk,pp,mean,vr,vf=CORE['fft_power'](cube,L)
                    km,pm,de,count,sums=CORE['spherical'](kk,pp,EDGES);check(np.all(count>0),'Empty bins')
                    with np.load(path,allow_pickle=False) as dd:
                        x=np.asarray(dd['xH_box'],float);check(x.shape==cube.shape and np.isfinite(x).all(),'Bad neutral fraction cube');xm=float(x.mean())
                    rr=hist[n,m,nu]
                    for value,key in [(mean,'dTb_mean_mK'),(np.sqrt(vr),'dTb_std_mK'),(xm,'xHI_volume')]:
                        check(np.isclose(value,float(rr[key]),rtol=2e-5,atol=1e-6),'History/cube mismatch: '+key)
                    z=dict(k=km,power=pm,delta2=de,count=count)
                    target.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(target,**z)
                    info.update(mean_Tb_mK=mean,xHI_volume=xm,variance_mK2=vr,variance_fft_mK2=vf,interpolation=status(de),spectrum_sha256=sha(target))
                    json_write(jp,info)
                    del cube,kk,pp,x
                spectra[n]=z;stats.append(dict(frequency_MHz=nu,model=m,N=n,mean_Tb_mK=info['mean_Tb_mK'],xHI_volume=info['xHI_volume'],variance_mK2=info['variance_mK2'],interpolation=info['interpolation']))
                for i,k in enumerate(z['k']):allrows.append(dict(frequency_MHz=nu,model=m,N=n,k_low_Mpc_inv=EDGES[i],k_high_Mpc_inv=EDGES[i+1],k_mean_Mpc_inv=k,Delta2_mK2=z['delta2'][i],P_mK2_Mpc3=z['power'][i],fft_cells=int(z['count'][i])))
            lo,hi=spectra[96],spectra[128];check(np.allclose(lo['k'],hi['k'],rtol=1e-10) and np.array_equal(lo['count'],hi['count']),'Shell sampling differs')
            for i,k in enumerate(lo['k']):
                x,y=lo['delta2'][i],hi['delta2'][i];ratio=y/x if x>0 else None
                symmetric=2*abs(y-x)/(x+y) if x+y>0 else 0.
                compare.append(dict(frequency_MHz=nu,model=m,k_low_Mpc_inv=EDGES[i],k_high_Mpc_inv=EDGES[i+1],k_mean_Mpc_inv=k,Delta2_96_mK2=x,Delta2_128_mK2=y,power_ratio=ratio,symmetric_fractional_difference=symmetric,status='both_zero' if x==y==0 else ('one_zero' if min(x,y)==0 else 'positive'),complete_science_bin=bool(EDGES[i]>=.1-1e-12 and EDGES[i+1]<=.3+1e-12)))
            print(f'Extracted/resumed {nu:g} MHz {m}',flush=True)
    csv_write(a.out/'spectra.csv',allrows);csv_write(a.out/'snapshot_summary.csv',stats);csv_write(a.out/'resolution_map.csv',compare)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(4,2,figsize=(11,12),layout='constrained')
    for ax,m in zip(axs.flat,MODELS):
        data=np.array([[r['symmetric_fractional_difference'] for r in compare if r['model']==m and r['frequency_MHz']==nu and r['complete_science_bin']] for nu in sorted(a.frequencies)])
        im=ax.imshow(100*data,origin='lower',aspect='auto',extent=[.1,.3,0,len(a.frequencies)],vmin=0,vmax=30,cmap='magma')
        tick=np.arange(0,len(a.frequencies),max(1,len(a.frequencies)//6));ax.set_yticks(tick+.5,np.array(sorted(a.frequencies))[tick]);ax.set(title=m,xlabel='k [Mpc^-1]',ylabel='Frequency [MHz]');fig.colorbar(im,ax=ax,label='Symmetric difference [%]',extend='max')
    axs.flat[-1].axis('off');fig.savefig(a.out/'resolution_map.png',dpi=160);plt.close(fig)
    json_write(a.out/'COMPLETE.json',dict(snapshots=len(stats),frequencies=sorted(a.frequencies),rows=len(allrows),note='Output-grid sensitivity only; no convergence or phase-matching claim across frequency. Both-zero bins have difference zero, not evidence of accuracy.'))
    print('Finished:',a.out)
if __name__=='__main__':main()
