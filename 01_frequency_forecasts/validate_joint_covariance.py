#!/usr/bin/env python3
"""Optional: assess a USER-SUPPLIED joint bandpower covariance and template difference.
Input NPZ: delta (n,), covariance (n,n), labels (n,) Unicode; no pickle.
A JSON provenance file must specify origin, units, bands_MHz, window, covariance_model.
This checks algebra and conditioning, not physical validity. No covariance is invented.
"""
import argparse,json,runpy
from pathlib import Path
import numpy as np
G=runpy.run_path(str(Path(__file__).resolve().parent/'common.py'));globals().update({k:v for k,v in G.items() if not k.startswith('__')})
def joint(delta,C):
    check(C.shape==(len(delta),len(delta)) and np.isfinite(C).all() and np.isfinite(delta).all(),'Invalid joint arrays')
    check(np.allclose(C,C.T,rtol=1e-10,atol=1e-12*max(abs(C).max(),1e-300)),'Covariance not symmetric')
    check(np.all(np.diag(C)>0),'Nonpositive variance');scale=np.sqrt(np.diag(C));R=C/scale[:,None]/scale[None,:];cond=float(np.linalg.cond(R))
    check(cond<1e12,'Ill-conditioned correlation; revise estimator, do not silently pseudoinvert')
    L=np.linalg.cholesky(R);whitened=np.linalg.solve(L,delta/scale)
    return float(np.linalg.norm(whitened)),cond

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--provenance',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    provenance=json.loads(a.provenance.read_text())
    for k in ['origin','units','bands_MHz','window','covariance_model']:check(k in provenance and bool(provenance[k]),'Missing provenance '+k)
    z=columns(a.input);check(len(set(z['labels'].tolist()))==len(z['delta']),'Labels missing/duplicated');D,cond=joint(z['delta'],z['covariance'])
    check(not a.out.exists(),'Choose a new output file');json_write(a.out,dict(template_distance=D,correlation_condition_number=cond,provenance=provenance,input_sha256=sha(a.input),note='Algebraic covariance check only. Not a calibrated significance or a physical validation of supplied covariance.'))
if __name__=='__main__':main()
