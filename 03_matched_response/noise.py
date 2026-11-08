"""Local documented SKAO station SEFD lookup; retain three time-dependent curves."""
from pathlib import Path
import numpy as np
from common import load_reference,sha,dump
SOURCE='https://developer.skao.int/projects/ska-ost-senscalc/en/11.1.1/_modules/ska_ost_senscalc/low/sefd_lookup.html'
def make(database,centre,out):
 import h5py
 from scipy.interpolate import UnivariateSpline
 import astropy.units as u
 from astropy.coordinates import EarthLocation,SkyCoord,AltAz
 from astropy.time import Time
 from astropy.utils import iers
 iers.conf.auto_download=False
 pr,ref,cat,dirs=load_reference();site=pr['site']
 location=EarthLocation.from_geodetic(site['longitude_deg']*u.deg,site['latitude_deg']*u.deg,site['height_m']*u.m)
 times=Time(pr['mjd_times'],format='mjd',scale='utc',location=location)
 sky=SkyCoord(pr['phase_centre_ra_deg']*u.deg,pr['phase_centre_dec_deg']*u.deg,frame='icrs')
 altaz=sky.transform_to(AltAz(obstime=times,location=location,pressure=0*u.hPa));lst=times.sidereal_time('apparent').hour
 if np.any(altaz.alt.deg<=0):raise ValueError('Pointing below horizon')
 nu=(centre-4)*1e6+(np.arange(128)+.5)*62500;curves=[];selections=[]
 with h5py.File(database,'r') as h:
  dims=h['dimensions'];f=np.asarray(dims['frequency']);az=np.asarray(dims['azimuth']);za=np.asarray(dims['zenith_angle']);lt=np.asarray(dims['lst']);pol=np.asarray(dims['polarisation'])
  if h['sefd'].shape!=(len(lt),len(az),len(za),len(f),2):raise ValueError('Unexpected SEFD table axes')
  if not np.all(np.diff(f)>0) or f[0]>nu.min()/1e6 or f[-1]<nu.max()/1e6:raise ValueError('Frequency coverage/unit mismatch')
  for t in range(3):
   zen=90-float(altaz.alt.deg[t]);zi=0 if zen<=2.5 else int(np.argmin(abs(za-zen)));ai=0 if zen<=2.5 else int(np.argmin(abs(az-float(altaz.az.deg[t]))));ti=int(np.argmin(abs(lt-lst[t])))
   xy=np.asarray(h['sefd'][ti,ai,zi,:,:]);fine=np.array([UnivariateSpline(f,xy[:,p],s=0)(nu/1e6) for p in range(2)])
   if not np.isfinite(fine).all() or np.any(fine<=0) or not np.isfinite(xy).all() or np.any(xy<=0):raise ValueError('Invalid SEFD or spline overshoot')
   curves.append(fine);selections.append(dict(time_index=t,mjd=float(times.mjd[t]),az_deg=float(altaz.az.deg[t]),el_deg=float(altaz.alt.deg[t]),lst_hour=float(lst[t]),table_indices=[ti,ai,zi],selected_lst=float(lt[ti]),selected_az=float(az[ai]),selected_za=float(za[zi])))
 xy=np.array(curves);sefd=.5*np.sqrt(xy[:,0]**2+xy[:,1]**2)
 var=sefd**2/(62500*10.) # E|n_I|^2, per 10-second split; no extra polarization factor.
 out.mkdir(parents=True,exist_ok=False)
 np.savez_compressed(out/'noise.npz',frequencies_hz=nu,sefd_XY_Jy=xy,station_stokes_I_sefd_Jy=sefd,complex_variance_I_Jy2=var,fit_inverse_variance_Jy_minus2=1/var)
 regression=None
 if centre==140:
  regression=float(np.max(abs(var/ref['complex_variance_I_Jy2']-1)))
  if regression>1e-7:raise ValueError('140 MHz table/geometry noise regression differs: '+str(regression)+'; preserve results, do not silently replace weights')
 dump(out/'configuration.json',dict(frequency_MHz=centre,database_path=str(database.resolve()),database_sha256=sha(database),source=SOURCE,method='Single nearest LST/Az/ZA entry at each exact original midpoint; cubic frequency interpolation per X/Y, then Stokes I combination',selection=selections,polarisation_labels=[v.decode() if isinstance(v,bytes) else str(v) for v in pol],units='SEFD Jy; frequency MHz in table; complex variance Jy^2',variance='E|n_I|^2=SEFD_I^2/(62500 Hz * 10 s); Re/Im variance each half; no additional npol factor',station_assumption='Identical sensitivity at every selected station; table default AAVS2 model is not measured station-specific AA4 sensitivity',split_seconds=10,original_140_noise_max_relative_error=regression,fit_weights='inverse complex variance; final mode/time/baseline averaging remains equal as in the archived estimator'))
 dump(out/'summary.json',dict(status='COMPLETE_LOCAL_TABLE_NOISE',sefd_I_range_Jy=[float(sefd.min()),float(sefd.max())],variance_range_Jy2=[float(var.min()),float(var.max())]))
