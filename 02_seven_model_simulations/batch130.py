#!/usr/bin/env python3
"""130-MHz cosmological input generation and audit, before new instrument response."""
import argparse,copy,hashlib,io,json,os,subprocess,sys,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
from seven_models import verify,sha,write
NU21=1420.40575177
MODELS=['F','Hminus','Hplus','Mminus','Mplus','Iminus','Iplus']
SEEDS=[12345,23456,34567]

def config(m,space,seed):
 c=copy.deepcopy(json.loads((ROOT/'reference_140_configs.json').read_text())[m][space])
 c.update(frequency_MHz=130,redshift=NU21/130-1,random_seed=seed)
 return c

def geometry_check(dist,red,idx):
 dist=np.asarray(dist);red=np.asarray(red);idx=np.asarray(idx)
 if not np.isfinite(dist).all() or not np.isfinite(red).all():raise ValueError('Nonfinite coordinates')
 if len(dist)!=len(red) or len(dist)!=len(idx):raise ValueError('Coordinate lengths differ')
 np.testing.assert_allclose(np.diff(dist),150/128,rtol=0,atol=1e-7)
 freq=NU21/(1+red)
 if not np.all(np.diff(freq)<0):raise ValueError('Unexpected frequency ordering')
 if freq.min()>126.03125 or freq.max()<133.96875:raise ValueError('Lightcone does not cover all 128 channel centres; do not extrapolate')
 if not np.all(idx==idx.astype(int)) or len(set(idx.tolist()))!=len(idx) or idx.min()<0 or idx.max()>=128:raise ValueError('Repeated/wrapped coeval slices: do not change box or wrap silently')
 return {'source_slices':len(dist),'source_frequency_MHz':[float(freq.min()),float(freq.max())],'channel_centres_MHz':[126.03125,133.96875],'unique_coeval_slices':True}

def preflight(path):
 import py21cmfast as p21
 import astropy.units as u
 from astropy.cosmology import z_at_value
 if p21.__version__!='3.4.0':raise ValueError('Use existing eor-models environment with py21cmfast exactly3.4.0')
 c=config('F','lightcone_real',12345);cell=150/128*u.Mpc;cosmo=p21.CosmoParams(c['cosmo_params']).cosmo
 lc=p21.RectilinearLightconer.with_equal_cdist_slices(min_redshift=NU21/134-1,max_redshift=NU21/126-1,resolution=cell,cosmo=cosmo,quantities=('brightness_temp',),index_offset=0)
 pixels=list(lc.get_lc_distances_in_pixels(cell))
 distances=np.asarray(lc.lc_distances.to_value(u.Mpc),dtype=float)
 if len(distances)!=len(pixels):raise ValueError('Lightconer physical/pixel coordinate lengths differ')
 marker=np.arange(128)[None,None,:]
 idx=np.array([int(lc.coeval_subselect(v,marker,cell).item()) for v in pixels])
 red=np.array([float(z_at_value(cosmo.comoving_distance,x*u.Mpc,zmin=0,zmax=50)) for x in distances])
 report=geometry_check(distances,red,idx)
 report.update(status='PASS_GEOMETRY_PREFLIGHT',py21cmfast=p21.__version__,source_distance_Mpc=distances.tolist(),source_redshifts=red.tolist(),coeval_slice_indices=idx.tolist(),scope='No simulation run; exact installed lightconer slice mapping inspected')
 write(path,report);print(json.dumps({k:v for k,v in report.items() if not isinstance(v,list) or len(v)<4},indent=2))

def plan(out):
 if (out/'batch.json').exists():raise ValueError('Batch already planned; use run/review')
 jobs={}
 for seed in SEEDS:
  entries=[];folder=out.resolve()/f'seed{seed}'
  for m in MODELS:
   for sp in ['real','lightcone_real']:
    entries.append(dict(model=m,space=sp,config=config(m,sp,seed),action='generate',cube=str(folder/'production/130MHz'/m/sp/'cube.npz'),verified=None))
  file=folder/'plan.json';write(file,dict(schema=1,scope='130 MHz, 150 cMpc,128/384, matched real-space templates',jobs=entries,candidate_errors=[]));jobs[str(seed)]=str(file)
 write(out/'batch.json',{'seeds':SEEDS,'models':MODELS,'frequency_MHz':130,'plans':jobs,'jobs':42,'reference_140_configs_sha256':sha(ROOT/'reference_140_configs.json')})
 print('Planned42 jobs; no simulation executed:',out)

def checked_plans(out):
 b=json.loads((out/'batch.json').read_text())
 if b['seeds']!=SEEDS or b['models']!=MODELS or b['frequency_MHz']!=130 or b['reference_140_configs_sha256']!=sha(ROOT/'reference_140_configs.json'):raise ValueError('Batch changed')
 for seed in SEEDS:
  p=Path(b['plans'][str(seed)]);jobs=json.loads(p.read_text())['jobs']
  if len(jobs)!=14 or {(j['model'],j['space']) for j in jobs}!={(m,s) for m in MODELS for s in ['real','lightcone_real']}:raise ValueError('Unexpected job set')
  for j in jobs:
   if j['config']!=config(j['model'],j['space'],seed):raise ValueError('Plan altered; fixed reference required')
  yield seed,p,jobs

def run(out,execute):
 if execute:
  preflight(out/'geometry_preflight.json')
 for seed,p,jobs in checked_plans(out):
  cmd=[sys.executable,str(ROOT/'vendor/seven_models.py'),'run','--plan',str(p)]
  if execute:cmd+=['--execute']
  subprocess.run(cmd,check=True)

def review(out):
 reference=None;rows=[];density={}
 for seed,p,jobs in checked_plans(out):
  cmd=[sys.executable,str(ROOT/'vendor/seven_models.py'),'review','--plan',str(p)]
  signals=p.parent/'seven_models_signal_inputs.zip'
  if not signals.exists():cmd+=['--export-signals']
  subprocess.run(cmd,check=True)
  with zipfile.ZipFile(signals) as z:prov={(r['model'],r['kind']):r for r in json.loads(z.read('provenance.json'))}
  density[seed]=set()
  for j in jobs:
   v=verify(j['cube'],j['config'])
   if v['sha256']!=prov[j['model'],j['space']]['cube_sha256']:raise ValueError('Stale signal archive; preserve old export elsewhere and rerun review')
   if j['space']=='real':density[seed].add(v['density_npy_sha256'])
   else:
    with np.load(j['cube'],allow_pickle=False) as d:g=[d[k].copy() for k in ['slice_distance_Mpc','slice_redshifts','coeval_slice_indices']]
    row=geometry_check(*g);rows.append(dict(seed=seed,model=j['model'],**row))
    if reference is None:reference=g
    else:
     for x,y in zip(g,reference):np.testing.assert_allclose(x,y,rtol=0,atol=1e-9)
  if len(density[seed])!=1:raise ValueError('Within-seed density differs')
 if len({next(iter(d)) for d in density.values()})!=3:raise ValueError('Different seeds unexpectedly have identical density fingerprints')
 nu=126e6+(np.arange(128)+.5)*62500;dist,red,idx=reference;freq=NU21*1e6/(1+red);target=np.interp(nu,freq[::-1],dist[::-1]);X=float(np.interp(130e6,freq[::-1],dist[::-1]));dl=(150/128)/X
 axis=(np.arange(128)+.5-64)*dl;l,m=np.meshgrid(axis,axis,indexing='ij');n=np.sqrt(1-l*l-m*m);directions=np.c_[l.ravel(),m.ravel(),n.ravel()]
 np.savez_compressed(out/'signal_geometry130.npz',signal_directions_lmn=directions,signal_solid_angle_sr=dl*dl/directions[:,2],frequencies_hz=nu,target_distance_Mpc=target,source_distance_Mpc=dist,source_redshifts=red,source_slice_indices=idx)
 report={'status':'PASS_SEVEN_MODEL_THREE_SEED_130_INPUT_AUDIT','frequency_MHz':130,'seeds':SEEDS,'models':MODELS,'jobs':42,'density_sha256':{str(s):next(iter(h)) for s,h in density.items()},'geometry_rows':rows,'central_distance_Mpc':X,'signal_direction_cosine_span_deg':float(np.rad2deg(150/X)),'geometry_sha256':sha(out/'signal_geometry130.npz'),'window_note':'Same comoving box does not give identical angular support at130 and140MHz; record the difference and choose common angular support explicitly before any frequency ranking','next':'Evaluate 130 MHz Jones response on these directions plus the unchanged catalogue at the original times, and obtain 130 MHz noise weights before forward comparison. Do not reuse 140 MHz beam/noise arrays.','limits':'Input generation only; finite real-space150cMpc patch. No noise, visibilities, foreground fit, detection/distance, likelihood or cross-frequency ranking result.'}
 write(out/'review130.json',report)
 target=out/'seven130_input_review.zip'
 with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
  for name in ['batch.json','geometry_preflight.json','review130.json','signal_geometry130.npz']:z.write(out/name,name)
  for seed in SEEDS:z.write(out/f'seed{seed}/seven_models_review.zip',f'seed{seed}_review.zip')
  for f in ROOT.rglob('*'):
   if f.is_file() and '__pycache__' not in f.parts:z.write(f,'code/'+str(f.relative_to(ROOT)))
 print('Send this small archive:',target)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='cmd',required=True)
 q=sub.add_parser('preflight');q.add_argument('--out',type=Path,required=True)
 for command in ['plan','run','review']:
  q=sub.add_parser(command);q.add_argument('--batch',type=Path,required=True)
  if command=='run':q.add_argument('--execute',action='store_true')
 a=p.parse_args()
 if a.cmd=='preflight':preflight(a.out)
 elif a.cmd=='plan':plan(a.batch)
 elif a.cmd=='run':run(a.batch,a.execute)
 else:review(a.batch)
