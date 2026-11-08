#!/usr/bin/env python3
"""Dry-run first. Sequential matched real/RSD coevals or selected-band real lightcones.
Pinned 21cmFAST 3.4.0 only. Original files are never modified.
"""
import argparse,hashlib,json,os,sys,time,importlib.metadata
from pathlib import Path
NU21=1420.40575177

def jobs(stage):
 if stage=='evolution120':return [(m,nu,'real') for m in ['F','Mplus'] for nu in [116,118,120,122,124]]
 if stage=='lightcone120':return [(m,120,sp) for m in ['F','Mplus'] for sp in ['real','lightcone_real']]
 if stage=='rsd':return [(m,nu,space) for m in ['F','Mplus'] for nu in [130,140] for space in ['real','redshift']]
 if stage=='evolution':return [(m,nu,'real') for m in ['F','Mplus'] for nu in [126,128,130,132,134,136,138,140,142,144]]
 if stage=='lightcone':return [(m,nu,'lightcone_real') for m in ['F','Mplus'] for nu in [130,140]]

def config(model,nu,space,threads,smoke=False):
 ref=json.loads((Path(__file__).parent/'reference_configs.json').read_text())[model]
 ref=json.loads(json.dumps(ref));ref['user_params']['N_THREADS']=threads
 if smoke:ref['user_params'].update(HII_DIM=32,DIM=96)
 ref['flag_options'].update(APPLY_RSDS=space=='redshift',SUBCELL_RSD=space=='redshift')
 # Use original stored centre redshifts to avoid introducing a redshift regression offset.
 candidates=[z for z in ref['redshifts'] if abs(NU21/(1+z)-nu)<1e-4]
 redshift=candidates[0] if candidates else NU21/nu-1
 return dict(model=model,frequency_MHz=nu,redshift=redshift,space=space,required_version='3.4.0',random_seed=12345,user_params=ref['user_params'],cosmo_params=ref['cosmo_params'],flag_options=ref['flag_options'],astro_params=ref['astro_params_log10_inputs'],los_axis=2,smoke=smoke,evolution_anchor_redshift=min(ref['redshifts']))

def run_one(c,dest):
 import numpy as np
 import py21cmfast as p21
 if p21.__version__!='3.4.0':raise RuntimeError('Requires exactly py21cmfast 3.4.0; no v4 substitution')
 flags=p21.FlagOptions(c['flag_options'])
 for k in ['APPLY_RSDS','SUBCELL_RSD']:
  if bool(getattr(flags,k))!=c['flag_options'][k]:raise RuntimeError('RSD flag was not accepted: '+k)
 dest.mkdir(parents=True,exist_ok=True);manifest=dest/'requested.json'
 if manifest.exists() and json.loads(manifest.read_text())!=c:raise RuntimeError('Existing run configuration differs: '+str(dest))
 result=dest/'cube.npz'
 if (dest/'COMPLETE.json').exists():
  done=json.loads((dest/'COMPLETE.json').read_text())
  if not result.exists() or hashlib.sha256(result.read_bytes()).hexdigest()!=done['cube_sha256']:raise RuntimeError('Completed cube missing or changed')
  print('Already complete',dest,flush=True);return
 manifest.write_text(json.dumps(c,indent=2));start=time.monotonic()
 (dest.parents[2]/'cache'/c['model']/c['space']).mkdir(parents=True,exist_ok=True)
 kwargs=dict(user_params=c['user_params'],cosmo_params=c['cosmo_params'],astro_params=c['astro_params'],flag_options=c['flag_options'],random_seed=c['random_seed'],direc=str((dest.parents[2]/'cache'/c['model']/c['space']).resolve()),write=True,regenerate=False)
 if c['space']=='lightcone_real':
  import astropy.units as u
  lcdef=p21.RectilinearLightconer.with_equal_cdist_slices(min_redshift=NU21/(c['frequency_MHz']+4)-1,max_redshift=NU21/(c['frequency_MHz']-4)-1,resolution=c['user_params']['BOX_LEN']/c['user_params']['HII_DIM']*u.Mpc,cosmo=p21.CosmoParams(c['cosmo_params']).cosmo,quantities=('brightness_temp',),index_offset=0)
  result_obj=p21.run_lightcone(redshift=c['evolution_anchor_redshift'],lightconer=lcdef,**kwargs)
  data={'brightness_temp':np.asarray(result_obj.brightness_temp,dtype=np.float32),'slice_redshifts':np.asarray(result_obj.lightcone_redshifts),'slice_distance_Mpc':np.asarray(result_obj.lightcone_distances.to_value(u.Mpc))}
  cell=c['user_params']['BOX_LEN']/c['user_params']['HII_DIM']*u.Mpc
  marker=np.arange(c['user_params']['HII_DIM'])[None,None,:]
  data['coeval_slice_indices']=np.array([int(lcdef.coeval_subselect(d,marker,cell).item()) for d in lcdef.get_lc_distances_in_pixels(cell)])
 else:
  result_obj=p21.run_coeval(redshift=[c['redshift'],c['evolution_anchor_redshift']],**kwargs)[0]
  data={k:np.asarray(getattr(result_obj,k),dtype=np.float32) for k in ['brightness_temp','density','xH_box','Tk_box','Ts_box','velocity']}
 for k,v in data.items():
  if not np.isfinite(v).all():raise RuntimeError('Nonfinite output '+k)
 n=c['user_params']['HII_DIM'];L=c['user_params']['BOX_LEN']
 if data['brightness_temp'].shape[:2]!=(n,n):raise RuntimeError('Unexpected grid')
 native={k:dict(getattr(result_obj,k).defining_dict) for k in ['user_params','cosmo_params','astro_params','flag_options']};native['global_params']=result_obj.global_params
 (dest/'resolved.json').write_text(json.dumps(native,indent=2,default=str))
 with (dest/'cube.tmp').open('wb') as stream:np.savez_compressed(stream,**data,box_len_cMpc=L,cell_size_cMpc=L/n,redshift=c['redshift'],random_seed=c['random_seed'],coordinate_space=c['space'],los_axis=2)
 (dest/'cube.tmp').replace(result)
 env={d.metadata['Name']:d.version for d in importlib.metadata.distributions() if d.metadata.get('Name')}
 (dest/'environment.json').write_text(json.dumps(dict(python=sys.version,packages=env),indent=2))
 (dest/'COMPLETE.json').write_text(json.dumps(dict(elapsed_seconds=time.monotonic()-start,cube_sha256=hashlib.sha256(result.read_bytes()).hexdigest(),status='COMPLETE_SIMULATION',physics_validated=False),indent=2))
 print('Completed',dest,flush=True)

def main(a):
 # One subprocess per snapshot bounds Python-held cube memory and permits resumable execution.
 import subprocess
 selected=jobs(a.stage)
 if a.smoke:selected=[('F',130,'real'),('F',130,'redshift')]
 for m,nu,space in selected:
  c=config(m,nu,space,a.threads,a.smoke);dest=a.out/('smoke' if a.smoke else 'production')/f'{nu}MHz'/m/space
  if a.stage=='lightcone120':c['user_params'].update(BOX_LEN=168.75,HII_DIM=144,DIM=432)
  print(m,nu,space,'->',dest,flush=True)
  if a.execute:
   env=os.environ.copy();env.update(OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
   subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',json.dumps(c),'--out',str(dest)],check=True,env=env)
 print('Done.' if a.execute else 'Plan only. Add --execute to run; --smoke is installation testing only.')
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['rsd','evolution','lightcone','evolution120','lightcone120'],default='rsd');p.add_argument('--out',type=Path,required=True);p.add_argument('--threads',type=int,default=4);p.add_argument('--execute',action='store_true');p.add_argument('--smoke',action='store_true');p.add_argument('--worker',help=argparse.SUPPRESS)
 a=p.parse_args()
 if a.threads<1:p.error('threads must be positive')
 if a.worker:run_one(json.loads(a.worker),a.out)
 else:main(a)
