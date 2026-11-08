"""Shared strict numerical helpers; compatible with python -I via runpy."""
import csv, hashlib, json, runpy
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
CORE=runpy.run_path(str(HERE/'vendor/signal_core.py'))
OLD=runpy.run_path(str(HERE/'vendor/reforecast_resolution.py'))
MODELS=CORE['MODELS']; CASES=OLD['CASES']; EDGES=OLD['EDGES']; NU21=CORE['NU21']
def check(ok,msg):
    if not ok: raise ValueError(msg)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def json_write(p,x):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def csv_write(p,rows):
    check(bool(rows),'No rows: '+str(p))
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def tag(nu):return f'{nu:g}MHz'
def code_hash():return {str(p.relative_to(HERE)):sha(p) for p in sorted(HERE.rglob('*.py'))}
def interp(k,x,y):
    k,x,y=map(np.asarray,(k,x,y))
    check(np.all(np.isfinite(x)) and np.all(np.diff(x)>0),'Invalid k support')
    check(np.all(np.isfinite(y)&(y>=0)),'Negative/nonfinite signal')
    check(np.all((k>=x[0])&(k<=x[-1])),'Signal extrapolation forbidden')
    if np.all(y==0):return np.zeros_like(k,dtype=float)
    if np.all(y>0):return np.exp(np.interp(np.log(k),np.log(x),np.log(y)))
    # No arbitrary floor and no bridging across zero bins with a log fit.
    return np.interp(k,x,y)
def status(y):return 'zero_power' if np.all(y==0) else ('positive_loglog' if np.all(y>0) else 'mixed_zero_linear')
def metric(s,e,mask):
    check(np.all(np.isfinite(e)&(e>0)),'Invalid error')
    return float(np.linalg.norm(s[mask]/e[mask])) if np.any(mask) else None
def columns(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def matched(a,b):
    for k in ['model','required_version','random_seed','redshifts','coordinate_space','cosmo_params','flag_options','astro_params_log10_inputs']:
        check(a[k]==b[k],'Grid config mismatch: '+k)
    for k in set(a['user_params'])|set(b['user_params']):
        if k not in ['HII_DIM','N_THREADS']:check(a['user_params'].get(k)==b['user_params'].get(k),'Grid parameter mismatch: '+k)
def signal_file(root,nu,m,n):return Path(root)/'snapshots'/tag(nu)/f'{m}_N{n}.npz'
def signal_at(root,nu,m,n,k):
    p=signal_file(root,nu,m,n);j=json.loads(p.with_suffix('.json').read_text())
    check(j['spectrum_sha256']==sha(p),'Spectrum checksum mismatch: '+str(p))
    z=columns(p);return interp(k,z['k'],z['delta2'])
def config_checks(cfg,m,n,ref):
    check(cfg['model']==m and cfg['coordinate_space']=='real' and cfg['required_version']=='3.4.0','Wrong simulation type')
    check(cfg['user_params']['HII_DIM']==n and cfg['user_params']['DIM']==384,'Expected N=96/128 and DIM=384')
    old=ref['signal_metadata'];common=old['common_configuration']
    for k in ['cosmo_params','flag_options','random_seed']:check(cfg[k]==common[k],'Reference config mismatch: '+k)
    for k,v in common['user_params'].items():
        if k not in ['HII_DIM','DIM','N_THREADS']:check(cfg['user_params'].get(k)==v,'Reference numerical mismatch: '+k)
    check(cfg['astro_params_log10_inputs']==old['requested_astro_inputs'][m],'Model astrophysics mismatch: '+m)
def ref_meta(folder):
    p=Path(folder)/'forecast_metadata.json';r=json.loads(p.read_text())
    check(r.get('thermal_only') is False,'Need original full forecast')
    check(abs(r['instrument']['reference']['frequency_MHz']-130)<1e-5,'Need 130 MHz reference')
    return r
