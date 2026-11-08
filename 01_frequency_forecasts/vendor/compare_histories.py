#!/usr/bin/env python3
"""Compare matched histories; cross-check exact 130 MHz rows against saved cubes."""
import argparse, csv, hashlib, json, runpy
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CORE = runpy.run_path(str(Path(__file__).with_name('signal_core.py')))
FIELDS = {'xHI_volume':'xH_box', 'Tk_mean_K':'Tk_box',
          'Ts_mean_K':'Ts_box', 'dTb_mean_mK':'brightness_temp',
          'dTb_std_mK':'brightness_temp'}

def write_csv(path, rows):
    with Path(path).open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def compare_configs(low, high, model):
    for key in ('model','required_version','random_seed','redshifts','coordinate_space',
                'cosmo_params','flag_options','astro_params_log10_inputs'):
        if key not in low or key not in high or low[key]!=high[key]:
            raise ValueError(f'{model}: mismatch or missing config {key}')
    if low['model']!=model or low['coordinate_space']!='real':
        raise ValueError('Model/space mismatch')
    ul,uh=low['user_params'],high['user_params']
    if (ul['HII_DIM'],uh['HII_DIM'],ul['DIM'],uh['DIM'])!=(64,128,384,384):
        raise ValueError('Expected 64/128 grids at fixed DIM=384')
    for k in set(ul)|set(uh):
        if k not in ('HII_DIM','N_THREADS') and ul.get(k)!=uh.get(k):
            raise ValueError(f'{model}: mismatched user parameter {k}')

def history(folder):
    if not (folder/'COMPLETE.json').exists():
        raise ValueError(f'Missing COMPLETE.json: {folder}')
    cfg=json.loads((folder/'requested_config.json').read_text())
    with (folder/'history.csv').open() as f: rows=list(csv.DictReader(f))
    if not rows: raise ValueError('Empty history')
    required=['redshift','frequency_MHz',*FIELDS]
    arr={k:np.array([float(r[k]) for r in rows]) for k in required}
    if any(not np.isfinite(v).all() for v in arr.values()):raise ValueError('Non-finite history')
    order=np.argsort(arr['frequency_MHz']);arr={k:v[order] for k,v in arr.items()}
    nu=arr['frequency_MHz']
    if np.any(np.diff(nu)<=0):raise ValueError('Duplicate/unordered frequency')
    if not np.allclose(nu,CORE['NU21']/(1+arr['redshift']),rtol=0,atol=1e-4):
        raise ValueError('Frequency/redshift mismatch')
    if len(rows)!=len(cfg['redshifts']) or not np.allclose(sorted(arr['redshift']),sorted(cfg['redshifts']),rtol=0,atol=1e-7):
        raise ValueError('History does not match requested redshift schedule')
    cube,L,c,info=CORE['load_model'](folder.parent,folder.name,130.)
    ix=np.flatnonzero(np.abs(nu-130)<1e-4)
    if len(ix)!=1:raise ValueError('Need exactly one 130 MHz history row (no interpolation)')
    ix=int(ix[0])
    with np.load(info['cube_path'],allow_pickle=False) as d:
        for key,field in FIELDS.items():
            x=np.asarray(d[field],dtype=float)
            if x.shape!=cube.shape or not np.isfinite(x).all():raise ValueError('Invalid '+field)
            val=x.std() if key=='dTb_std_mK' else x.mean()
            if not np.isclose(val,arr[key][ix],rtol=2e-5,atol=1e-6):
                raise ValueError(f'{folder}: cube/history mismatch for {key}')
            if field=='xH_box' and (x.min() < -1e-5 or x.max()>1.00001):raise ValueError('xHI bounds')
            if field in ('Tk_box','Ts_box') and x.min()<=0:raise ValueError('Temperature bounds')
    return cfg,arr,info

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--low-root',type=Path,required=True);p.add_argument('--high-root',type=Path,required=True)
    p.add_argument('--models',nargs='+',default=['F']);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if len(set(a.models))!=len(a.models):p.error('Duplicate models')
    if a.out.exists() and any(a.out.iterdir()):raise FileExistsError('Choose a new output directory')
    a.out.mkdir(parents=True,exist_ok=True)
    target=[];allrows=[];sources={}
    for model in a.models:
        cl,lo,il=history(a.low_root/model);ch,hi,ih=history(a.high_root/model)
        compare_configs(cl,ch,model)
        if lo['frequency_MHz'].shape!=hi['frequency_MHz'].shape or not np.allclose(lo['frequency_MHz'],hi['frequency_MHz'],rtol=0,atol=1e-4):
            raise ValueError('Histories use different frequency schedules')
        for i,nu in enumerate(lo['frequency_MHz']):
            for field in FIELDS:
                l=float(lo[field][i]);h=float(hi[field][i])
                row=dict(model=model,frequency_MHz=float(nu),quantity=field,value_64=l,value_128=h,
                         difference_128_minus_64=h-l,
                         fractional_difference=(h-l)/abs(l) if abs(l)>1e-12 else '')
                allrows.append(row)
                if abs(nu-130)<1e-4:target.append(row)
        fig,axs=plt.subplots(2,3,figsize=(13,7),layout='constrained')
        for ax,field in zip(axs.ravel(),FIELDS):
            ax.plot(lo['frequency_MHz'],lo[field],'o--',ms=3,label='64³, DIM=384')
            ax.plot(hi['frequency_MHz'],hi[field],'o-',ms=3,label='128³, DIM=384')
            ax.axvline(130,color='gray',lw=1);ax.set(xlabel='Frequency [MHz]',ylabel=field);ax.grid(alpha=.2)
            if field in ('Tk_mean_K','Ts_mean_K'):ax.set_yscale('log')
        axs[0,0].legend(fontsize=8);axs.ravel()[-1].axis('off')
        fig.suptitle(f'{model}: matched histories (no smoothing/interpolation)')
        fig.savefig(a.out/f'history_comparison_{model}.png',dpi=180);plt.close(fig)
        sources[model]={'low_cube':il,'high_cube':ih}
        print(f'PASS {model}: configuration and 130 MHz cube/history consistency')
    write_csv(a.out/'history_comparison_all_frequencies.csv',allrows)
    write_csv(a.out/'history_comparison_130MHz.csv',target)
    for r in target:
        print(f"{r['model']:7s} {r['quantity']:16s} 64={r['value_64']:.7g} 128={r['value_128']:.7g} delta={r['difference_128_minus_64']:+.6g}")
    (a.out/'history_checks.json').write_text(json.dumps(dict(status='PASS',sources=sources,
        note='Differences are diagnostics, not an automatic convergence verdict. Fractional differences near zero brightness are unstable.'),indent=2))

if __name__=='__main__':main()
