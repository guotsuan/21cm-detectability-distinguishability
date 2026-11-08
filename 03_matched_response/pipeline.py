#!/usr/bin/env python3
"""Matched 130/140 MHz response pipeline. No cosmological simulations are launched."""
import argparse,copy,hashlib,io,json,os,shutil,sys,zipfile,csv
from pathlib import Path
import numpy as np
from common import ROOT,sha,dump,load,load_reference,fringe,gaussian,C
SEEDS=[12345,23456,34567];MODELS=['F','Hminus','Hplus','Mminus','Mplus','Iminus','Iplus']
class Source:
 def __init__(self,p):
  self.p=Path(p);self.z=zipfile.ZipFile(p) if self.p.is_file() else None
 def read(self,n):return self.z.read(n) if self.z else (self.p/n).read_bytes()
 def obj(self,n):return json.loads(self.read(n))
def identify(p):
 ref=json.loads((ROOT/'inputs/reference_140_configs.json').read_text())
 with zipfile.ZipFile(p) as z:
  c=json.loads(z.read('F/real/requested.json'));nu=int(c['frequency_MHz']);seed=int(c['random_seed'])
  if nu not in [130,140] or seed not in SEEDS:raise ValueError('Not a requested frequency/seed')
  for m in MODELS:
   for k in ['real','lightcone_real']:
    expected=copy.deepcopy(ref[m][k]);expected['random_seed']=seed
    if nu==130:expected.update(frequency_MHz=130,redshift=1420.40575177/130-1)
    if json.loads(z.read(m+'/'+k+'/requested.json'))!=expected:raise ValueError('Unmatched config '+m+'/'+k)
  prov=json.loads(z.read('provenance.json'))
  if len(prov)!=14 or {(r['model'],r['kind']) for r in prov}!={(m,k) for m in MODELS for k in ['real','lightcone_real']}:raise ValueError('Incomplete archive')
  signature=hashlib.sha256(json.dumps(sorted([(r['model'],r['kind'],r['array_sha256']) for r in prov]),sort_keys=True).encode()).hexdigest()
 return nu,seed,signature

def discover(a):
 found={};errors=[]
 for p in sorted(a.root.rglob('seven_models_signal_inputs*.zip')):
  try:nu,seed,sig=identify(p)
  except Exception as e:errors.append(dict(path=str(p),error=str(e)));continue
  found.setdefault((nu,seed),[]).append((p,sig))
 entries=[]
 for nu in [130,140]:
  for seed in SEEDS:
   candidates=found.get((nu,seed),[])
   if not candidates:raise ValueError(f'Missing matched {nu} MHz seed {seed} signal archive. Keep seven_models_signal_inputs.zip exports on disk.')
   if len({v[1] for v in candidates})!=1:raise ValueError(f'Different signal versions for {nu}/{seed}: '+str([str(v[0]) for v in candidates])+'; specify one by removing unneeded exports from the search root or use a directory containing symlinks to the six intended exports.')
   p,sig=min(candidates,key=lambda v:(len(str(v[0])),str(v[0])))
   entries.append(dict(frequency_MHz=nu,seed=seed,signals=str(p.resolve()),array_manifest_sha256=sig))
 tables=list(a.root.rglob('ska_station_sensitivity_AAVS2.h5'))
 if not tables:raise ValueError('Missing ska_station_sensitivity_AAVS2.h5 under project root')
 if len({sha(p) for p in tables})!=1:raise ValueError('Multiple different SEFD tables; choose one in manifest explicitly')
 db=min(tables,key=lambda p:len(str(p)))
 result=dict(schema=1,signals=entries,database=str(db.resolve()),database_sha256=sha(db),rejected_signal_candidates=errors,settings=dict(interpolation='linear',angular_support='fixed 150 cMpc square; different angular span at each frequency',mc=2048))
 dump(a.manifest,result);print('Created',a.manifest)
 for r in entries:print(r['frequency_MHz'],r['seed'],r['signals'])

def setup_context(out,nu):
 folder=out/f'{nu}MHz';os.environ['EOR_MATCHED_CONTEXT']=str((folder/'noise').resolve());return folder

def prepare_stage(a,cfg):
 from noise import make
 import prepare
 if sha(Path(cfg['database']))!=cfg['database_sha256']:raise ValueError('SEFD table changed after discovery')
 for nu in a.frequencies:
  folder=setup_context(a.out,nu);folder.mkdir(parents=True,exist_ok=True)
  noise=folder/'noise'
  if (noise/'summary.json').exists():
   if json.loads((noise/'summary.json').read_text())['status']!='COMPLETE_LOCAL_TABLE_NOISE' or json.loads((noise/'configuration.json').read_text())['database_sha256']!=cfg['database_sha256']:raise ValueError('Incomplete/stale noise output')
  elif noise.exists():raise ValueError('Incomplete noise directory; preserve it and choose fresh --out')
  else:make(Path(cfg['database']),nu,noise)
  reference=None
  for entry in [r for r in cfg['signals'] if r['frequency_MHz']==nu]:
   signals=Path(entry['signals']);actual=identify(signals)
   if actual!=(nu,entry['seed'],entry['array_manifest_sha256']):raise ValueError('Signal archive/config changed')
   dest=folder/f"seed{entry['seed']}"/'prepared'
   if (dest/'metadata.json').exists():
    meta=json.loads((dest/'metadata.json').read_text())
    if meta['signal_archive_sha256']!=sha(signals) or meta['physical_sha256']!=sha(dest/'physical.npz'):raise ValueError('Prepared resume hash mismatch')
   elif dest.exists():raise ValueError('Partial prepared directory; preserve it and choose fresh --out')
   else:prepare.main(argparse.Namespace(case='baseline128',signals=signals,out=dest))
   with np.load(dest/'physical.npz') as z:
    keys=['directions','solid_angle_sr','frequencies_hz','target_distance_Mpc','source_slice_indices','source_distance_Mpc','source_redshifts']
    geom={k:z[k].copy() for k in keys}
   if reference is None:reference=geom
   else:
    for k in keys:np.testing.assert_array_equal(geom[k],reference[k],err_msg='Across-seed geometry mismatch '+k)
  np.savez_compressed(folder/'geometry.npz',**reference)
  dump(folder/'prepared_complete.json',dict(status='COMPLETE',frequency_MHz=nu,seeds=SEEDS,geometry_sha256=sha(folder/'geometry.npz')))

def collect_beams(source,dirs,out):
 lookup={};records=source.obj('manifest.json')
 for rec in records:
  raw=source.read(rec['file'])
  if hashlib.sha256(raw).hexdigest()!=rec['sha256']:raise ValueError('Beam hash mismatch')
  key=rec['time'],rec['channel'],rec['model']
  if key in lookup:raise ValueError('Duplicate beam key')
  with np.load(io.BytesIO(raw)) as z:v=z['normalised_I_power'].copy()
  if v.shape!=(len(dirs),) or not np.isfinite(v).all() or np.any(v<0):raise ValueError('Bad beam response')
  lookup[key]=v
 expected={(t,ch,m) for t in range(3) for ch in range(128) for m in ['dipole','isotropic']}|{(t,ch,'gaussian') for t in range(3) for ch in [0,63,127]}
 if set(lookup)!=expected:raise ValueError('Expected all 777 beam entries')
 p=np.array([[[lookup[t,ch,m] for ch in range(128)] for t in range(3)] for m in ['dipole','isotropic']])
 np.savez_compressed(out/'beam_power.npz',power=p)

def beam_stage(a):
 import run_beams
 for nu in a.frequencies:
  folder=setup_context(a.out,nu);out=folder/'beams';out.mkdir(exist_ok=True)
  geom=np.load(folder/'geometry.npz');dirs=geom['directions'];prep=folder/'seed12345/prepared';meta=json.loads((prep/'metadata.json').read_text())
  if (out/'complete.json').exists():
   done=json.loads((out/'complete.json').read_text())
   if done['geometry_sha256']!=sha(folder/'geometry.npz') or done['power_sha256']!=sha(out/'beam_power.npz'):raise ValueError('Beam resume hash mismatch')
   print('Reuse checked complete beams:',nu);continue
  if nu==140 and a.reuse140_beams:
   if not a.reuse140_prepared:raise ValueError('--reuse140-prepared also required')
   old=Source(a.reuse140_prepared);om=old.obj('metadata.json');raw=old.read('physical.npz')
   if hashlib.sha256(raw).hexdigest()!=om['physical_sha256']:raise ValueError('Old preparation hash mismatch')
   with np.load(io.BytesIO(raw)) as z:
    for key in geom.files:np.testing.assert_array_equal(geom[key],z[key],err_msg='140 beam geometry reuse mismatch '+key)
   source=Source(a.reuse140_beams);bc=source.obj('configuration.json');pr,_,_,_=load_reference()
   if om['original_input_provenance']!=pr:raise ValueError('Old prepared observing provenance differs')
   if bc['site']!=pr['site'] or bc['dipole_length_m']!=C/(2*140e6):raise ValueError('Old beam station settings differ')
   if bc['input_sha256']!={**pr['input_sha256'],'physical.npz':om['physical_sha256']} or not bc['full']:raise ValueError('Old beam provenance mismatch')
   if source.obj('summary.json')['status']!='BEAMS_COMPLETE_REVIEW_REQUIRED':raise ValueError('Old beam run incomplete')
   collect_beams(source,dirs,out);dump(out/'configuration.json',bc)
   origin=dict(type='exact-geometry authenticated 140 MHz reuse',source=str(a.reuse140_beams))
  else:
   pr,d,cat,_=load();pr['input_sha256']={**pr['input_sha256'],'physical.npz':meta['physical_sha256']}
   run_beams.load=lambda:(pr,d,cat,dirs)
   run_beams.main(argparse.Namespace(full=True),out)
   collect_beams(Source(out),dirs,out);origin=dict(type='new OSKAR evaluation',frequency_MHz=nu)
  dump(out/'complete.json',dict(status='COMPLETE_REFERENCE_STATION_BEAMS',geometry_sha256=sha(folder/'geometry.npz'),power_sha256=sha(out/'beam_power.npz'),origin=origin,scope='Shared reference S8-1 array; fixed physical dipole length c/(2*140MHz); on-axis-normalised unpolarized power shape; not actual per-station SKALA beams'))
  print('Beam response complete',nu,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['discover','prepare','beams','analyze','package']);p.add_argument('--root',type=Path,default=Path.home()/'eor_model_test');p.add_argument('--manifest',type=Path,default=Path('matched_response_inputs.json'));p.add_argument('--out',type=Path,default=Path('matched_response130_140_v1'));p.add_argument('--frequencies',type=int,nargs='+',choices=[130,140],default=[130,140]);p.add_argument('--reuse140-beams',type=Path);p.add_argument('--reuse140-prepared',type=Path);p.add_argument('--mc',type=int,default=2048);a=p.parse_args();a.out=a.out.resolve()
 if a.stage=='discover':discover(a)
 elif a.stage=='prepare':prepare_stage(a,json.loads(a.manifest.read_text()))
 elif a.stage=='beams':beam_stage(a)
 elif a.stage=='analyze':
  from analyze import run
  if a.mc<128:p.error('--mc must be at least 128; use 2048 for review')
  for nu in a.frequencies:run(setup_context(a.out,nu),a.mc)
 else:
  from report import package
  package(a.out)
