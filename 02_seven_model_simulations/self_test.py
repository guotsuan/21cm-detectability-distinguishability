#!/usr/bin/env python3
"""Tests exact configuration changes, plan and band-coverage rejection; no simulation."""
import tempfile,json,subprocess,sys
from pathlib import Path
import numpy as np
from batch130 import config,plan,checked_plans,geometry_check,ROOT,MODELS,SEEDS,NU21
old=json.loads((ROOT/'reference_140_configs.json').read_text())
for seed in SEEDS:
 for m in MODELS:
  for sp in ['real','lightcone_real']:
   c=config(m,sp,seed);assert c['frequency_MHz']==130 and c['random_seed']==seed
   assert c['redshift']==NU21/130-1
   c['frequency_MHz']=old[m][sp]['frequency_MHz'];c['redshift']=old[m][sp]['redshift'];c['random_seed']=old[m][sp]['random_seed'];assert c==old[m][sp]
with tempfile.TemporaryDirectory() as t:
 p=Path(t);plan(p);assert sum(len(j) for _,_,j in checked_plans(p))==42
 subprocess.run([sys.executable,str(ROOT/'batch130.py'),'run','--batch',str(p)],check=True,stdout=subprocess.DEVNULL)
 d=9000+np.arange(128)*150/128;r=NU21/np.linspace(134,126,128)-1;idx=np.arange(128)
 assert geometry_check(d,r,idx)['unique_coeval_slices']
 for rr,ii in [(NU21/np.linspace(141,139,128)-1,idx),(r,np.zeros(128,int))]:
  try:geometry_check(d,rr,ii)
  except ValueError:pass
  else:raise AssertionError('Bad coverage or duplicate indices accepted')
 # The vendor coverage verifier must use130, not hard-coded140 endpoints.
 from seven_models import verify,write,sha
 c=config('F','lightcone_real',12345);out=p/'fixture';out.mkdir()
 np.savez_compressed(out/'cube.npz',brightness_temp=np.zeros((128,128,128),np.float32),box_len_cMpc=150.,cell_size_cMpc=150/128,random_seed=12345,los_axis=2,coordinate_space='lightcone_real',slice_distance_Mpc=d,slice_redshifts=r,coeval_slice_indices=idx)
 write(out/'requested.json',c);write(out/'COMPLETE.json',{'cube_sha256':sha(out/'cube.npz')});verify(out/'cube.npz',c)
print('PASS: only frequency/redshift/seed changed,42 jobs,dry run,130-MHz band validation,no duplicate slices. No physical simulation executed.')
