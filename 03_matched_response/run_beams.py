#!/usr/bin/env python3
"""Evaluate common reference station beams on the original FIXED celestial sky."""
import argparse,csv,json,shutil,subprocess,traceback,sys
from pathlib import Path
import numpy as np
from common import *

def main(a,out):
 import oskar
 exe=shutil.which('oskar_sim_beam_pattern')
 if not exe:raise RuntimeError('oskar_sim_beam_pattern not in PATH')
 pr,d,cat,dirs=load();rd=radec(dirs,pr['phase_centre_ra_deg'],pr['phase_centre_dec_deg'])
 sky=out/'directions_ra_dec.txt';np.savetxt(sky,rd,fmt='%.17g')
 tm=out/'reference_station.tm';(tm/'station000').mkdir(parents=True,exist_ok=True)
 site=pr['site'];(tm/'position.txt').write_text(' '.join(str(site[k]) for k in ['longitude_deg','latitude_deg','height_m'])+'\n');(tm/'layout.txt').write_text('0 0 0\n')
 rows=list(csv.DictReader((ROOT/'inputs/station_antennas.csv').open()));xyz=np.array([[float(r[k]) for k in ['east','north','up']] for r in rows]);assert xyz.shape==(256,3)
 if not all(float(r['amp'])==1. and float(r['phase'])==0. for r in rows):raise ValueError('Expected uniform unit element amplitudes and zero phase offsets')
 np.savetxt(tm/'station000/layout.txt',xyz,fmt='%.17g')
 for name in ['feed_angle_x.txt','feed_angle_y.txt']:np.savetxt(tm/'station000'/name,np.zeros(256))
 channels=list(range(128)) if a.full else [0,63,127];manifest=[];gauss_errors=[]
 dump(out/'configuration.json',{'full':a.full,'channels':channels,'source_count':len(dirs),'source_order':'centre,24 fixed foreground sources,prepared physical signal pixels','input_sha256':pr['input_sha256'],'oskar_version':str(oskar.oskar_version_string()),'dipole_length_m':C/(2*140e6),'common_reference_station':pr['station_mapping'],'beam_normalisation':'raw Jones retained; scalar power shape divided by on-axis trace(J J^H)/2. No Jones inverse or polarized calibration. Noise kept in the original conditional Stokes-I units.','time_handling':'zero-duration beam evaluation at original sample midpoints; visibility fringes use archived UVW and original sinc averaging','site':site})
 for t,mjd in enumerate(pr['mjd_times']):
  for ch in channels:
   nu=float(d['frequencies_hz'][ch]);models=['isotropic','dipole']+(['gaussian'] if ch in [0,63,127] else [])
   for model in models:
    job=out/f't{t}_c{ch:03d}_{model}';job.mkdir(exist_ok=True);ini=job/'settings.ini';ini.touch()
    p={'simulator':{'double_precision':True,'use_gpus':False,'num_devices':1,'keep_log_file':False,'max_sources_per_chunk':32768},
       'observation':{'mode':'Tracking','start_time_utc':mjd,'length':0.,'num_time_steps':1,'phase_centre_ra_deg':pr['phase_centre_ra_deg'],'phase_centre_dec_deg':pr['phase_centre_dec_deg'],'start_frequency_hz':nu,'num_channels':1},
       'telescope':{'input_directory':str(tm),'station_type':'Gaussian beam' if model=='gaussian' else 'Aperture array','pol_mode':'Full' if model=='dipole' else 'Scalar','normalise_beams_at_phase_centre':False,'aperture_array/array_pattern/enable':True,'aperture_array/array_pattern/normalise':True,'aperture_array/element_pattern/enable_numerical':False,'aperture_array/element_pattern/normalise':False,'aperture_array/element_pattern/functional_type':'Dipole' if model=='dipole' else 'Isotropic (unpolarised)','aperture_array/element_pattern/dipole_length':C/(2*140e6),'aperture_array/element_pattern/dipole_length_units':'Metres','gaussian_beam/fwhm_deg':4.,'gaussian_beam/ref_freq_hz':140e6},
       'beam_pattern':{'station_ids':'0','coordinate_frame':'Equatorial','coordinate_type':'Sky model','sky_model/file':str(sky),'root_path':str(job/'beam'),'output/separate_time_and_channel':True,'output/average_single_axis':'None','station_outputs/text_file/raw_complex':True}}
    cache=job/'record.json'
    if cache.exists():
     rec=json.loads(cache.read_text())
     if json.loads((job/'requested_settings.json').read_text())!=p or sha(job/'beam.npz')!=rec['sha256']:raise ValueError('Resume cache changed; choose a new run directory')
     with np.load(job/'beam.npz') as h:
      if h['normalised_I_power'].shape!=(len(dirs),):raise ValueError('Resume source count changed')
     manifest.append(rec);dump(out/'manifest.json',manifest);continue
    dump(job/'requested_settings.json',p);s=oskar.SettingsTree('oskar_sim_beam_pattern',str(ini));s.from_dict(p)
    with (job/'oskar.log').open('w') as log:subprocess.run([exe,str(ini)],cwd=job,stdout=log,stderr=subprocess.STDOUT,check=True)
    files=list(job.glob('*RAW_COMPLEX.txt'))
    if len(files)!=1:raise ValueError('Unexpected raw output files')
    x=np.loadtxt(files[0],ndmin=2);ncol=8 if model=='dipole' else 2
    if x.shape!=(len(dirs),ncol) or not np.isfinite(x).all():raise ValueError(f'Unexpected beam output {x.shape}, expected {(len(dirs),ncol)}')
    z=x[:,0::2]+1j*x[:,1::2]
    power=np.sum(abs(z)**2,axis=1)/(2 if model=='dipole' else 1)
    if power[0]<=1e-12:raise ValueError('Nonpositive/vanishing on-axis response')
    raw0=float(power[0]);power=power/raw0
    if model=='gaussian':
     err=float(max(abs(power-gaussian(nu,dirs))));gauss_errors.append(err)
     if err>2e-6:raise ValueError(f'Equatorial Gaussian regression failed: {err}')
    np.savez_compressed(job/'beam.npz',raw_complex=z,normalised_I_power=power,raw_on_axis_power=raw0)
    # Preserve exact text-file hash; losslessly retain parsed numbers in compressed output.
    rawhash=sha(files[0]);files[0].unlink()
    manifest.append({'time':t,'channel':ch,'model':model,'file':str((job/'beam.npz').relative_to(out)),'sha256':sha(job/'beam.npz'),'raw_text_sha256':rawhash})
    dump(job/'record.json',manifest[-1]);dump(out/'manifest.json',manifest)
   print(f'Completed pointing {t+1}/3 channel {ch+1}/128',flush=True)
 dump(out/'summary.json',{'status':'BEAMS_COMPLETE_REVIEW_REQUIRED','full':a.full,'entries':len(manifest),'max_gaussian_power_regression_error':max(gauss_errors,default=0.),'not_validated':['station-specific layouts','physical SKALA response','intra-integration beam evolution','sky covariance']})
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--full',action='store_true');a=p.parse_args();out=a.out.expanduser().resolve()
 if out.exists():p.error('Choose a new --out directory')
 out.mkdir(parents=True);shutil.copy2(__file__,out/'run_beams.py');shutil.copy2(ROOT/'common.py',out/'common.py');shutil.copy2(ROOT/'inputs/provenance.json',out/'input_provenance.json')
 code=0
 try:main(a,out)
 except Exception:
  code=1;err=traceback.format_exc();(out/'error.txt').write_text(err);dump(out/'summary.json',{'status':'ERROR','error':err});print(err,file=sys.stderr)
 finally:archive(out)
 sys.exit(code)
