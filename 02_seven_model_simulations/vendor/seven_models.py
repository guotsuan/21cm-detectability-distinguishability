#!/usr/bin/env python3
"""Seven-model physical-input audit, resumable generation and review export.
Stage 1 only: no OSKAR response or model-significance calculation.
"""
import argparse, copy, hashlib, io, json, os, subprocess, sys, zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
VARIANTS={'F':{},'Hminus':{'L_X':39.0},'Hplus':{'L_X':41.0},
          'Mminus':{'M_TURN':8.0},'Mplus':{'M_TURN':9.5},
          'Iminus':{'F_ESC10':-1.5},'Iplus':{'F_ESC10':-.5}}

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(p)
def config(model,space,threads=4):
 c=json.loads((ROOT/'reference_F140.json').read_text())
 c['model']=model;c['space']=space;c['user_params']['N_THREADS']=threads
 c['astro_params'].update(VARIANTS[model]);return c
def signature(c):
 c=copy.deepcopy(c);c.get('user_params',{}).pop('N_THREADS',None)
 return c

def verify(path,c):
 import numpy as np
 path=Path(path)
 req=json.loads((path.parent/'requested.json').read_text())
 if signature(req)!=signature(c):raise ValueError('configuration mismatch')
 done=json.loads((path.parent/'COMPLETE.json').read_text())
 digest=sha(path)
 if digest!=done.get('cube_sha256'):raise ValueError('COMPLETE checksum mismatch')
 with np.load(path,allow_pickle=False) as d:
  b=d['brightness_temp'];n=c['user_params']['HII_DIM'];L=c['user_params']['BOX_LEN']
  if b.ndim!=3 or b.shape[:2]!=(n,n) or not np.isfinite(b).all():raise ValueError('invalid temperature grid')
  if float(d['box_len_cMpc'])!=L or float(d['cell_size_cMpc'])!=L/n:raise ValueError('box/cell mismatch')
  if int(d['random_seed'])!=c['random_seed'] or int(d['los_axis'])!=2:raise ValueError('seed/LOS mismatch')
  if str(d['coordinate_space'])!=c['space']:raise ValueError('space mismatch')
  info={'shape':list(b.shape),'temperature_mean_mK':float(b.mean()),'temperature_std_mK':float(b.std())}
  if c['space']=='real':
   if b.shape!=(n,n,n) or d['density'].shape!=b.shape or not np.isfinite(d['density']).all():raise ValueError('invalid coeval/density')
   with zipfile.ZipFile(path) as z:info['density_npy_sha256']=hashlib.sha256(z.read('density.npy')).hexdigest()
  else:
   dist=d['slice_distance_Mpc'];red=d['slice_redshifts'];idx=d['coeval_slice_indices']
   freq=1420.40575177/(1+red)
   if len(dist)!=b.shape[2] or not np.isfinite(dist).all() or not np.isfinite(red).all():raise ValueError('invalid coordinates')
   if not np.allclose(np.diff(dist),L/n,rtol=0,atol=1e-7) or not np.all(np.diff(freq)<0):raise ValueError('radial ordering/spacing')
   if len(idx)!=len(dist) or not np.all(idx==idx.astype(int)) or len(set(idx.tolist()))!=len(idx) or idx.min()<0 or idx.max()>=n:raise ValueError('invalid/repeated coeval slice indices')
   if freq.min()>c['frequency_MHz']-4+0.03125 or freq.max()<c['frequency_MHz']+4-0.03125:raise ValueError('lightcone does not cover target channels')
   info.update(frequency_min_MHz=float(freq.min()),frequency_max_MHz=float(freq.max()))
 return dict(path=str(path.resolve()),sha256=digest,**info)

def plan(a):
 out=a.out.resolve()
 if (out/'plan.json').exists():raise ValueError('Plan already exists. Reuse it with run/review, or choose a new output.')
 candidates=[];errors=[]
 for root in a.search:
  if not root.is_dir():raise ValueError('Missing search directory: '+str(root))
  for base,dirs,files in os.walk(root):
   dirs[:]=[d for d in dirs if not d.startswith('.') and d not in {'cache','__pycache__','node_modules'}]
   if 'requested.json' not in files or 'cube.npz' not in files or 'COMPLETE.json' not in files:continue
   p=Path(base)/'cube.npz'
   try:c=json.loads((p.parent/'requested.json').read_text());candidates.append((p,c))
   except Exception as e:errors.append({'path':str(p),'reason':str(e)})
 jobs=[]
 for m in VARIANTS:
  for sp in ['real','lightcone_real']:
   c=config(m,sp,a.threads);match=None
   for p,found in candidates:
    if signature(found)!=signature(c):continue
    try:match=verify(p,c);break
    except Exception as e:errors.append({'path':str(p),'reason':str(e)})
   dest=out/'production'/'140MHz'/m/sp/'cube.npz'
   jobs.append({'model':m,'space':sp,'config':c,'action':'reuse' if match else 'generate','cube':match['path'] if match else str(dest),'verified':match})
   print(jobs[-1]['action'].upper(),m,sp,flush=True)
 write(out/'plan.json',{'schema':1,'scope':'140 MHz, seven original templates, real space, seed12345','jobs':jobs,'candidate_errors':errors})
 print('Plan:',out/'plan.json')
 print('Generations required:',sum(j['action']=='generate' for j in jobs),'of 14')

def run(a):
 p=json.loads(a.plan.read_text())
 for j in p['jobs']:
  c=j['config'];cube=Path(j['cube'])
  if j['action']=='reuse' or (cube.parent/'COMPLETE.json').exists():
   v=verify(cube,c)
   if j['action']=='reuse' and v['sha256']!=j['verified']['sha256']:raise ValueError('Reused file changed since planning')
   print('VERIFIED',j['model'],j['space'],flush=True);continue
  if not a.execute:print('WOULD GENERATE',j['model'],j['space']);continue
  env=os.environ.copy();env.update(OMP_NUM_THREADS=str(c['user_params']['N_THREADS']),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  subprocess.run([sys.executable,str(ROOT/'physics_worker.py'),'--worker',json.dumps(c),'--out',str(cube.parent)],env=env,check=True)
  verify(cube,c)
 print('Run finished.' if a.execute else 'Dry run only. Add --execute for simulation.')

def review(a):
 import numpy as np
 p=json.loads(a.plan.read_text());rows=[];ref=None;densities=[];errors=[]
 for j in p['jobs']:
  try:
   v=verify(j['cube'],j['config'])
   if j['action']=='reuse' and v['sha256']!=j['verified']['sha256']:raise ValueError('Reused file changed since planning')
   if j['space']=='real':densities.append(v['density_npy_sha256'])
   else:
    with np.load(j['cube'],allow_pickle=False) as d:
     geom=[d[k].copy() for k in ['slice_distance_Mpc','slice_redshifts','coeval_slice_indices']]
    if ref is None:ref=geom
    else:
     for x,y in zip(ref,geom):np.testing.assert_allclose(x,y,rtol=0,atol=1e-9)
   rows.append({'model':j['model'],'space':j['space'],**v})
  except Exception as e:errors.append({'model':j['model'],'space':j['space'],'reason':str(e)})
 if len(densities)==7 and len(set(densities))!=1:errors.append({'reason':'Density fingerprints differ across seven models. Investigate; do not disable this check.'})
 report={'status':'PASS_INPUT_AUDIT' if not errors and len(rows)==14 else 'INCOMPLETE_OR_FAILED','rows':rows,'errors':errors,'note':'Input consistency only; not numerical convergence, instrumental response or scientific significance.'}
 target=a.plan.parent/'seven_models_review.zip'
 with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
  z.writestr('review.json',json.dumps(report,indent=2));z.write(a.plan,'plan.json')
  for f in ROOT.glob('*.py'):z.write(f,'code/'+f.name)
  z.write(ROOT/'reference_F140.json','code/reference_F140.json')
  for j in p['jobs']:
   for name in ['requested.json','resolved.json','environment.json','COMPLETE.json']:
    f=Path(j['cube']).parent/name
    if f.exists():z.write(f,j['model']+'/'+j['space']+'/'+name)
 print(json.dumps(report,indent=2));print('Send this small archive:',target)
 if errors or len(rows)!=14:raise SystemExit(2)
 if a.export_signals:
  target=a.plan.parent/'seven_models_signal_inputs.zip'
  if target.exists():raise ValueError('Signal export exists; preserve or move it before exporting again')
  provenance=[]
  keys=['brightness_temp','box_len_cMpc','cell_size_cMpc','redshift','random_seed','coordinate_space','los_axis']
  with zipfile.ZipFile(target,'w',zipfile.ZIP_STORED) as z:
   for j,row in zip(p['jobs'],rows):
    names=keys+(['slice_distance_Mpc','slice_redshifts','coeval_slice_indices'] if j['space']=='lightcone_real' else [])
    with np.load(j['cube'],allow_pickle=False) as d:
     buf=io.BytesIO();np.savez_compressed(buf,**{k:d[k] for k in names})
    raw=buf.getvalue();base=j['model']+'/'+j['space']
    z.writestr(base+'/signal.npz',raw)
    z.write(Path(j['cube']).parent/'requested.json',base+'/requested.json')
    with zipfile.ZipFile(io.BytesIO(raw)) as q:hashes={name:hashlib.sha256(q.read(name)).hexdigest() for name in q.namelist()}
    provenance.append({'model':j['model'],'kind':j['space'],'source':j['cube'],'cube_sha256':row['sha256'],'array_sha256':hashes,'density_npy_sha256':row.get('density_npy_sha256')})
   z.writestr('provenance.json',json.dumps(provenance,indent=2))
  print('Large local signal archive (keep locally):',target)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='cmd',required=True)
 q=sub.add_parser('plan');q.add_argument('--search',type=Path,nargs='+',required=True);q.add_argument('--out',type=Path,required=True);q.add_argument('--threads',type=int,default=4)
 q=sub.add_parser('run');q.add_argument('--plan',type=Path,required=True);q.add_argument('--execute',action='store_true')
 q=sub.add_parser('review');q.add_argument('--plan',type=Path,required=True);q.add_argument('--export-signals',action='store_true')
 a=p.parse_args()
 if getattr(a,'threads',4)<1:p.error('threads must be positive')
 {'plan':plan,'run':run,'review':review}[a.cmd](a)
