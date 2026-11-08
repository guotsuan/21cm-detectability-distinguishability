#!/usr/bin/env python3
"""Fresh 21cmSense 3.0.0 weights at every band, gated by the saved 130 MHz forecast."""
import argparse,json,runpy,time
from pathlib import Path
import numpy as np
G=runpy.run_path(str(Path(__file__).resolve().parent/'common.py'));globals().update({k:v for k,v in G.items() if not k.startswith('__')})
def resolve_reference(arg,expected):
    if arg is not None:candidates=[arg]
    else:candidates=sorted({p.parent for p in Path.cwd().glob('**/reference_observation.json') if 'frequency_analysis' not in p.parts})
    matches=[]
    for p in candidates:
        f=p/'reference_observation.json'
        if f.exists() and json.loads(f.read_text())==expected and (p/'station_positions_native_xyz_m.npy').exists():matches.append(p)
    check(len(matches)==1,'Need one matching original reference directory. Pass --reference PATH containing reference_observation.json and station_positions_native_xyz_m.npy. Candidates: '+str(matches))
    return matches[0]
def compare_native(new,old,alpha):
    ix=np.lexsort((old['kpar_Mpc_inv'],old['kperp_Mpc_inv']))
    use=(new['k_Mpc_inv']>=.1)&(new['k_Mpc_inv']<=.3)
    jo=(old['k_Mpc_inv'][ix]>=.1)&(old['k_Mpc_inv'][ix]<=.3);ix=ix[jo]
    for name in ['k_Mpc_inv','kperp_Mpc_inv','kpar_Mpc_inv','thermal_sigma_mK2']:
        OLD['close'](new[name][use],old[name][ix],'130 MHz regression mismatch '+name,rtol=1e-6)
    OLD['close'](new['sample_coefficient'][use],alpha[ix],'130 MHz sample weights differ',rtol=1e-6)
    for case in CASES:check(np.array_equal(new[case][use],old[case][ix]),'130 MHz mask mismatch '+case)
    return dict(status='PASS',entries=int(use.sum()),thermal_and_sample_rtol=1e-6,masks='exact',extra_mode_factors_applied=False)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--forecast',type=Path,default=Path('forecast_130MHz_full'));p.add_argument('--reference',type=Path)
    p.add_argument('--frequencies',type=float,nargs='+',default=[130]);p.add_argument('--out',type=Path,default=Path('frequency_analysis/noise'))
    p.add_argument('--resume',action='store_true');p.add_argument('--preflight',action='store_true')
    a=p.parse_args();meta,old,alpha=OLD['saved_weights'](a.forecast);inst=meta['instrument'];cfg=inst['reference'];reference=resolve_reference(a.reference,cfg)
    check(len(set(a.frequencies))==len(a.frequencies),'Duplicate frequencies')
    B=cfg['bandwidth_MHz'];check(B==8 and cfg['n_channels']==128,'This package expects original 8 MHz / 128-channel reference')
    check(all(50+B/2<=nu<=200-B/2 for nu in a.frequencies),'8 MHz bands must be inside 50--200 MHz; use 55--195 MHz centres')
    from importlib.metadata import version
    check(version('21cmSense')=='3.0.0','Activate ska-forecast with 21cmSense==3.0.0')
    import astropy.units as u
    from astropy.cosmology import FlatLambdaCDM
    from astropy.cosmology.units import littleh,with_H0
    import py21cmsense as p21
    from py21cmsense import conversions as conv
    from py21cmsense.theory import TheoryModel
    c=cfg['cosmology'];cosmo=FlatLambdaCDM(H0=c['H0_km_s_Mpc'],Om0=c['Om0'],Ob0=c['Ob0'],Tcmb0=c['Tcmb0_K']*u.K);eq=with_H0(cosmo.H0)
    array=p21.Observatory.from_ska(cfg['array_template'],Trcv=cfg['receiver_temperature_K']*u.K)
    saved=np.load(reference/'station_positions_native_xyz_m.npy',allow_pickle=False)
    check(array.antpos.shape==saved.shape and np.allclose(array.antpos.to_value(u.m),saved,rtol=0,atol=1e-6),'Installed station layout differs from original')
    array=array.clone(baseline_filters=(p21.BaselineRange(bl_max=cfg['baseline_max_m']*u.m),))
    versions={v:version(v) for v in ['21cmSense','ska-ost-array-config','numpy','astropy']}
    job=dict(reference=cfg,station_sha256=sha(reference/'station_positions_native_xyz_m.npy'),time_step_seconds=inst['time_step_seconds'],buffer_h_Mpc_inv=inst['buffer_h_Mpc_inv'],versions=versions,code=code_hash(),original_forecast_sha256={f:sha(a.forecast/f) for f in ['forecast_metadata.json','forecast_native.npz','detection_snr.csv']})
    print('Reference:',reference.resolve());print('Versions:',versions);print('Band centres:',a.frequencies,flush=True)
    if a.preflight:return
    a.out.mkdir(parents=True,exist_ok=True);jp=a.out/'job.json'
    if jp.exists():check(a.resume and json.loads(jp.read_text())==job,'Noise config changed or --resume missing; choose new output')
    else:check(not any(a.out.iterdir()),'Nonempty output');json_write(jp,job)
    class PS(p21.PowerSpectrum):
        @property
        def X2Y(self):return conv.X2Y(self.observation.redshift,cosmo=self.cosmo,approximate=False)
    class Probe(TheoryModel):
        use_littleh=False
        def __init__(self,kind):self.kind=kind
        def delta_squared(self,z,k):
            k=np.asarray(k)
            return (np.ones_like(k)*self.kind if self.kind in [1,2] else 1+(k/.1)**2)*u.mK**2
    def extract(ps):
        r=ps.calculate_sensitivity_2d(compute_all=True);kz=ps.kparallel.to_value(1/u.Mpc,eq);rows=[]
        for key,tq in r['thermal'].items():
            kp=key.to_value(1/u.Mpc,eq);t=tq.to_value(u.mK**2);v=r['sample'][key].to_value(u.mK**2);b=r['both'][key].to_value(u.mK**2)
            ok=np.isfinite(t)&(t>0)&(kz>0)&(kp>0)
            if ok.any():rows.extend(zip(np.full(ok.sum(),kp),kz[ok],t[ok],v[ok],b[ok]))
        check(bool(rows),'No native samples');z=np.asarray(rows);z=z[np.lexsort((z[:,1],z[:,0]))]
        check(np.isfinite(z).all() and np.all(z[:,2:]>0),'Invalid native errors')
        OLD['close'](z[:,4],z[:,2]+z[:,3],'Public combined convention changed',rtol=1e-10)
        return z
    for nu in sorted(a.frequencies,key=lambda x:(x!=130,x)):
        folder=a.out/tag(nu);gate=a.out/'130MHz/COMPLETE.json'
        if nu!=130:check(gate.exists() and json.loads(gate.read_text())['regression_130']['status']=='PASS','First run and pass the 130 MHz regression')
        if (folder/'COMPLETE.json').exists():
            check(a.resume,'Output exists; use --resume');done=json.loads((folder/'COMPLETE.json').read_text());check(done['native_sha256']==sha(folder/'native.npz'),'Noise checkpoint corrupt');print('Verified resume:',nu,flush=True);continue
        folder.mkdir(parents=True,exist_ok=True);start=time.monotonic()
        obs=p21.Observation(observatory=array,frequency=nu*u.MHz,bandwidth=B*u.MHz,n_channels=cfg['n_channels'],n_days=cfg['n_days'],time_per_day=cfg['hours_per_day']*u.hour,lst_bin_size=cfg['lst_bin_hours']*u.hour,integration_time=inst['time_step_seconds']*u.s,phase_center_dec=cfg['phase_center_dec_deg']*u.deg,coherent=cfg['coherent'],cosmo=cosmo,use_approximate_cosmo=False,tsky_amplitude=260*u.K,tsky_ref_freq=150*u.MHz,spectral_index=2.6,max_chunk_mem_gb=.5)
        expected=cfg['receiver_temperature_K']+260*(nu/150)**-2.6
        check(np.isclose(obs.Tsys.to_value(u.K),expected,rtol=1e-6),'Tsys law mismatch')
        print(f'{nu:g} MHz: rebuilding UV coverage; Tsys={expected:.2f} K',flush=True);_=obs.uv_coverage
        base=PS(observation=obs,foreground_model='foreground_free',horizon_buffer=0*littleh/u.Mpc,theory_model=Probe(1));z=extract(base)
        for kind in [2,3]:
            zz=extract(base.clone(theory_model=Probe(kind)));OLD['close'](zz[:,:3],z[:,:3],'Probe altered coordinates/thermal errors',rtol=1e-10)
            kk=np.hypot(z[:,0],z[:,1]);scale=2 if kind==2 else 1+(kk/.1)**2
            OLD['close'](zz[:,3],z[:,3]*scale,'Sample coefficient not geometry-only',rtol=1e-6)
        beam=base.clone(foreground_model='optimistic');hor=base.clone(foreground_model='moderate',horizon_buffer=inst['buffer_h_Mpc_inv']*littleh/u.Mpc)
        dk=conv.dk_du(obs.redshift,cosmo=cosmo,approximate=False);umag=((z[:,0]/u.Mpc).to(dk.unit,eq)/dk).to_value(u.dimensionless_unscaled)
        bl=beam.horizon_limit(umag).to_value(1/u.Mpc,eq);hl=hor.horizon_limit(umag).to_value(1/u.Mpc,eq)
        native=dict(k_Mpc_inv=np.hypot(z[:,0],z[:,1]),kperp_Mpc_inv=z[:,0],kpar_Mpc_inv=z[:,1],thermal_sigma_mK2=z[:,2],sample_coefficient=z[:,3],no_wedge=np.ones(len(z),bool),beam_wedge=z[:,1]>=bl,horizon_buffer=z[:,1]>=hl)
        check(np.all(~native['horizon_buffer']|native['beam_wedge']),'Non-nested masks')
        regression=compare_native(native,old,alpha) if nu==130 else None
        np.savez_compressed(folder/'native.npz',**native)
        # Export only small UV bookkeeping, not a claim of independent 3D modes.
        pairs,counts=np.unique(np.column_stack([z[:,0],z[:,1]]),axis=0,return_counts=True)
        info=dict(frequency_MHz=nu,bandwidth_MHz=B,band_edges_MHz=[nu-B/2,nu+B/2],redshift=float(obs.redshift),Tsys_K=expected,native_entries=len(z),unique_coordinate_pairs=len(pairs),maximum_coordinate_multiplicity=int(counts.max()),probe_checks='1, 2 and scale-dependent positive spectra; public total=thermal+sample',regression_130=regression,native_sha256=sha(folder/'native.npz'),elapsed_seconds=time.monotonic()-start,observing_time='Original allocation per band; no assumption about simultaneous telescope bandwidth',note='Preserves duplicate native entries. No independent physical covariance, measured beam, foreground leakage or cross-frequency covariance validation.')
        json_write(folder/'COMPLETE.json',info);print('PASS',nu,'MHz',flush=True)
if __name__=='__main__':main()
