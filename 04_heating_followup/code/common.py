from pathlib import Path
import json,hashlib,sys,platform
import numpy as np
import pandas as pd
CASES=('no_wedge','beam_wedge','horizon_buffer')

def require(ok,msg):
    if not ok: raise ValueError(msg)
def output(path):
    path=Path(path);require(not path.exists() or not any(path.iterdir()),f'Output is not empty: {path}. Use a new RUN_TAG.')
    path.mkdir(parents=True,exist_ok=True);return path
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def native(root,nu):
    p=Path(root)/'forecasts'/f'{nu:g}MHz'/'native_N128.npz'
    require(p.is_file(),f'Missing native forecast: {p}')
    with np.load(p,allow_pickle=False) as f:d={k:f[k] for k in f.files}
    return d,p
def mask(d,case):return d[case].astype(bool)&(d['k_Mpc_inv']>=.1)&(d['k_Mpc_inv']<=.3)
def norm(x,s):return float(np.linalg.norm(x/s))
def fit(f,m,s):
    require(np.isfinite([f,m,s]).all() and np.all(s>0),'Invalid signals or sigma')
    den=np.dot(m/s,m/s);require(den>0,'Zero template norm: amplitude unidentifiable')
    a=max(0.,float(np.dot(m/s,f/s)/den))
    return a,norm(f-m,s),norm(f-a*m,s)
def finish(out,inputs,notes):
    data=dict(command=sys.argv,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,inputs={str(p):sha(p) for p in dict.fromkeys(inputs)},notes=notes)
    (out/'provenance.json').write_text(json.dumps(data,indent=2)+'\n')
