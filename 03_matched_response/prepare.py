#!/usr/bin/env python3
"""Prepare matched 130/140-MHz small-patch fluctuations without box tiling."""
import argparse,io,json,zipfile,hashlib,shutil
from pathlib import Path
import numpy as np
from sensitivity import CASES,cubic4,window_field
from common import load,dump,sha,ROOT,archive
MODELS=['F','Hminus','Hplus','Mminus','Mplus','Iminus','Iplus']

def interpolate(cube,dist,target):
 j=np.searchsorted(dist,target,side='right')-1;j=np.clip(j,0,len(dist)-2)
 if target.min()<dist[0] or target.max()>dist[-1]:raise ValueError('Target outside supplied radial support')
 w=(target-dist[j])/(dist[j+1]-dist[j])
 return np.moveaxis(cube[...,j]*(1-w)+cube[...,j+1]*w,-1,0)

def main(a):
 method,size=CASES[a.case]
 out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);pr,d,cat,olddirs=load();fields=[];means=[];reference=None
 with zipfile.ZipFile(a.signals) as z:
  prov=json.loads(z.read('provenance.json'))
  density=[r['density_npy_sha256'] for r in prov if r['kind']=='real']
  if len(density)!=7 or len(set(density))!=1:raise ValueError('Coeval density fingerprints differ')
  configs={};records={}
  for record in prov:
   key=(record['model'],record['kind'])
   if key in records:raise ValueError('Duplicate model/kind')
   raw=z.read('/'.join(key)+'/signal.npz')
   with zipfile.ZipFile(io.BytesIO(raw)) as inner:
    for name,h in record['array_sha256'].items():
     if hashlib.sha256(inner.read(name)).hexdigest()!=h:raise ValueError('Array hash mismatch')
   with np.load(io.BytesIO(raw),allow_pickle=False) as n:records[key]={k:n[k] for k in n.files}
   configs[key]=json.loads(z.read('/'.join(key)+'/requested.json'))
  if set(records)!={(m,k) for m in MODELS for k in ['real','lightcone_real']}:raise ValueError('Expected exactly seven matched model pairs')
  reference_configs=json.loads((ROOT/'inputs/reference_140_configs.json').read_text())
  seed=int(configs['F','real']['random_seed']);centre=int(round(np.mean(d['frequencies_hz'])/1e6))
  for model in MODELS:
   for kind in ['real','lightcone_real']:
    expected=dict(reference_configs[model][kind]);expected['random_seed']=seed
    if centre==130:expected.update(frequency_MHz=130,redshift=1420.40575177/130-1)
    if configs[model,kind]!=expected:raise ValueError('Unmatched reference configuration '+model+' '+kind)
  for m in MODELS:
   co=records[m,'real'];lc=records[m,'lightcone_real'];cfg=configs[m,'real']
   for k in ['user_params','cosmo_params','astro_params','flag_options','random_seed']:
    if cfg[k]!=configs[m,'lightcone_real'][k]:raise ValueError('Coeval/lightcone configuration mismatch '+k)
   if cfg['flag_options']['APPLY_RSDS'] or int(lc['los_axis'])!=2:raise ValueError('Expected real-space LOS axis 2')
   L=float(co['box_len_cMpc']);cell=float(co['cell_size_cMpc']);dist=lc['slice_distance_Mpc'];red=lc['slice_redshifts'];idx=lc['coeval_slice_indices'].astype(int);freq=1420.40575177e6/(1+red)
   if co['brightness_temp'].shape!=(128,128,128) or lc['brightness_temp'].shape!=(128,128,len(dist)):raise ValueError('Unexpected shape')
   if not np.allclose(np.diff(dist),cell) or not np.all(np.diff(freq)<0) or len(set(idx))!=len(idx) or idx.min()<0 or idx.max()>=128:raise ValueError('Invalid slice geometry')
   if reference is None:reference=(dist.copy(),red.copy(),idx.copy(),L,cell)
   else:
    for x,y in zip(reference,(dist,red,idx,L,cell)):np.testing.assert_allclose(x,y,rtol=0,atol=1e-10)
   target_nu=d['frequencies_hz']
   if target_nu.min()<freq.min() or target_nu.max()>freq.max():raise ValueError('Frequency centres outside LC coverage; no extrapolation allowed')
   target_dist=np.interp(target_nu,freq[::-1],dist[::-1])
   X=float(np.interp(np.mean(target_nu),freq[::-1],dist[::-1]))
   for kind,cube in [('coeval',co['brightness_temp'][:,:,idx]),('lightcone',lc['brightness_temp'])]:
    cube=np.asarray(cube,float)
    if not np.isfinite(cube).all():raise ValueError('Nonfinite temperatures')
    f=interpolate(cube,dist,target_dist) if method=='linear' else cubic4(cube,dist,target_dist)
    f,mean=window_field(f,size);means.append(mean);fields.append(f.reshape(128,-1).astype('float32'))
  # Both models use exactly the same geometry and cosmology.
  for k in ['user_params','cosmo_params','flag_options','random_seed']:
   for model in MODELS:
    if configs['F','real'][k]!=configs[model,'real'][k]:raise ValueError('Model configuration mismatch '+model+' '+k)
 dl=cell/X;axis=(np.arange(128)+.5-64)*dl;l,m=np.meshgrid(axis,axis,indexing='ij');n=np.sqrt(1-l*l-m*m);pix=np.c_[l.ravel(),m.ravel(),n.ravel()]
 dirs=np.vstack([olddirs[:25],pix]);solid=dl*dl/pix[:,2]
 np.savez_compressed(out/'physical.npz',temperatures_mK=np.array(fields),slice_means_mK=np.array(means),directions=dirs,solid_angle_sr=solid,frequencies_hz=d['frequencies_hz'],target_distance_Mpc=target_dist,source_slice_indices=idx,source_distance_Mpc=dist,source_redshifts=red)
 meta={'status':'PREPARED_PHYSICAL_SMALL_PATCH','signal_archive_sha256':sha(a.signals),'field_order':[m+'_'+kind for m in MODELS for kind in ['coeval','lightcone']],'model_order':MODELS,'physical_sha256':sha(out/'physical.npz'),'original_input_provenance':pr,'box_Mpc':L,'pixel_Mpc':cell,'central_distance_Mpc':X,'direction_cosine_span_deg':float(np.rad2deg(L/X)),'input_frequency_range_MHz':[float(freq.min()/1e6),float(freq.max()/1e6)],'target_frequency_range_MHz':[float(target_nu.min()/1e6),float(target_nu.max()/1e6)],'source_configs':{m:configs[m,'real'] for m in MODELS},'conventions':['brightness_temp interpreted as mK from existing 21cmFAST workflow','slice transverse means removed after identical radial interpolation for both controls; no global-signal claim','distance interpolated from supplied slice redshift/distance coordinates; no frequency reversal shortcut','linear interpolation in radial distance, identical for coeval/LC, no extrapolation','central-distance fixed angular mapping and pixel centres; not a varying-distance cone','zero outside 150-Mpc square, no tiling, new angular response required','one shared seed, real space, no sample-variance significance or full-FOV forecast']}
 meta['sensitivity']={'case':a.case,'radial_interpolation':method,'active_pixels_per_axis':size,'active_span_cMpc':L*size/128,'mean_definition':'arithmetic transverse mean over active pixels, after radial interpolation','angular_window':'central square crop; zeros outside; pixel coordinates and quadrature unchanged','normalization':'raw Jy2; no area or volume renormalization','scope':'window/interpolation sensitivity, not full-sky or numerical convergence'}
 meta['conventions']=[v for v in meta['conventions'] if not v.startswith('linear interpolation') and not v.startswith('zero outside')]
 dump(out/'metadata.json',meta);shutil.copy2(__file__,out/'prepare.py');print('Prepared:',out,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',choices=CASES,required=True);p.add_argument('--signals',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
