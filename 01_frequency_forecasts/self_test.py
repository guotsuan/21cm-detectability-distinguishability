#!/usr/bin/env python3
"""Synthetic numerical tests; no telescope calculation and no 21cmFAST needed."""
import json,runpy,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
G=runpy.run_path(str(HERE/'common.py'));CORE=G['CORE'];MODELS=G['MODELS'];CASES=G['CASES']
def run(*args):subprocess.run([sys.executable,'-I',str(HERE/args[0]),*map(str,args[1:])],check=True)
def main():
    CORE['self_test']()
    k=np.array([.1,.2,.3]);assert np.allclose(G['interp'](k,k,np.zeros(3)),0)
    assert np.isclose(G['interp'](np.array([.15]),k,np.array([0.,2.,3.]))[0],1)
    try:G['interp'](np.array([.4]),k,np.ones(3))
    except ValueError:pass
    else:raise AssertionError('Extrapolation allowed')
    R=runpy.run_path(str(HERE/'response_diagnostics.py'));A,w,g=R['operator'](32,[],'rectangular');assert np.allclose(A,np.eye(32),atol=1e-14)
    A,w,g=R['operator'](32,[3,4,15],'hann');rng=np.random.default_rng(13);x=rng.normal(size=32)+1j*rng.normal(size=32)
    assert np.allclose(A@np.fft.fft(x,norm='ortho'),np.fft.fft(w/g*x,norm='ortho'))
    O=R['overlap']([55,60,65],8,128);assert np.allclose(np.diag(O),1) and np.isclose(O[0,1],3/8) and O[0,2]==0
    J=runpy.run_path(str(HERE/'validate_joint_covariance.py'));D,c=J['joint'](np.array([3.,4.]),np.eye(2));assert np.isclose(D,5)
    try:J['joint'](np.ones(2),np.ones((2,2)))
    except (ValueError,np.linalg.LinAlgError):pass
    else:raise AssertionError('Singular joint covariance accepted')
    F=runpy.run_path(str(HERE/'forecast_bands.py'))
    native=dict(k_Mpc_inv=np.array([.15,.15,.26]),kperp_Mpc_inv=np.array([.1,.1,.1]),kpar_Mpc_inv=np.sqrt(np.array([.15,.15,.26])**2-.1**2),thermal_sigma_mK2=np.ones(3)*2,sample_coefficient=np.ones(3)*.5,no_wedge=np.ones(3,bool),beam_wedge=np.array([True,True,False]),horizon_buffer=np.array([False,False,True]))
    signals={n:{m:np.full(3,i+1.) for i,m in enumerate(MODELS)} for n in [96,128]}
    det,pairs,shifts,scale,err=F['evaluate'](130,native,signals)
    f=[r for r in det if r['model']=='F' and r['case']=='no_wedge' and r['kmax_Mpc_inv']==.3][0];assert np.isclose(f['snr_high'],np.sqrt(3)/2.5)
    pair=[r for r in pairs if r['model_a']=='F' and r['model_b']=='Hminus' and r['case']=='no_wedge' and r['kmax_Mpc_inv']==.3][0];assert np.isclose(pair['distance_high_fixed_F'],np.sqrt(3)/2.5)
    assert all(r['distance_using_low_covariance'] in (0,None) for r in shifts)
    assert any(r['snr_high'] is None for r in det)
    # Duplicate coordinates represent separate supplied native entries and are preserved.
    assert f['native_entries']==3
    signals[128]['F']=np.zeros(3);dd,*_=F['evaluate'](130,native,signals)
    assert all(r['snr_high']==0 for r in dd if r['model']=='F' and r['native_entries']>0)
    print('PASS: zero/mixed-zero handling, extrapolation rejection, duplicate-entry sums, analytic covariance metrics, exact response operator and overlap.')
    # Exercise actual command-line extraction and resumable checkpoints on finite cubes.
    with tempfile.TemporaryDirectory(prefix='eor-frequency-test-') as td:
        root=Path(td);ref=root/'forecast';ref.mkdir();common=dict(cosmo_params={'hlittle':.6766,'OMm':.30966,'OMb':.04897},flag_options={'APPLY_RSDS':False},random_seed=12345,user_params={'BOX_LEN':150.,'HII_DIM':64,'DIM':192,'N_THREADS':4})
        astro={m:{'test_amplitude':i+1} for i,m in enumerate(MODELS)};z=G['NU21']/130-1
        meta=dict(thermal_only=False,instrument={'reference':{'frequency_MHz':130}},signal_metadata={'common_configuration':common,'requested_astro_inputs':astro})
        (ref/'forecast_metadata.json').write_text(json.dumps(meta))
        for n in [96,128]:
            for i,m in enumerate(MODELS):
                f=root/f'n{n}'/m;f.mkdir(parents=True);cfg=dict(common,model=m,coordinate_space='real',required_version='3.4.0',redshifts=[z],astro_params_log10_inputs=astro[m]);cfg['user_params']=dict(common['user_params'],HII_DIM=n,DIM=384)
                (f/'requested_config.json').write_text(json.dumps(cfg));(f/'COMPLETE.json').write_text('{}')
                # All zero for Iplus tests explicit constant-field path.
                amplitude=0 if m=='Iplus' else (i+1)*(1 if n==96 else 1.01)
                cube=np.broadcast_to(amplitude*np.cos(2*np.pi*4*np.arange(n)/n)[:,None,None],(n,n,n)).copy()
                np.savez_compressed(f/f'coeval_z{z:.3f}.npz',brightness_temp=cube,xH_box=np.ones_like(cube)*.8,redshift=z,box_len_cMpc=150.,random_seed=12345,coordinate_space='real')
                (f/'history.csv').write_text('frequency_MHz,dTb_mean_mK,dTb_std_mK,xHI_volume\n130,'+str(cube.mean())+','+str(cube.std())+',0.8\n')
        args=['extract_spectra.py','--low-root',root/'n96','--high-root',root/'n128','--forecast',ref,'--frequencies','130','--out',root/'out']
        run(*args,'--preflight');run(*args);run(*args,'--resume')
        import csv
        with (root/'out/resolution_map.csv').open() as f:rows=list(csv.DictReader(f))
        r=[r for r in rows if r['model']=='F' and abs(float(r['k_low_Mpc_inv'])-.15)<1e-9][0]
        assert np.isclose(float(r['power_ratio']),1.01**2,rtol=1e-9)
        assert all(r['status']=='both_zero' and r['power_ratio']=='' for r in rows if r['model']=='Iplus')
        # Alter config; resume must reject rather than mixing different runs.
        cp=root/'n128/F/requested_config.json';conf=json.loads(cp.read_text());conf['random_seed']=999;cp.write_text(json.dumps(conf))
        fail=subprocess.run([sys.executable,'-I',str(HERE/args[0]),*map(str,args[1:]),'--resume'],capture_output=True,text=True);assert fail.returncode!=0
    print('PASS: command-line extraction, analytic inter-grid power ratio, zero-power snapshots, resume and changed-config rejection.')
    # End-to-end forecast CLI against independently calculated flat-spectrum values.
    with tempfile.TemporaryDirectory(prefix='eor-forecast-test-') as td:
        root=Path(td);spec=root/'spectra';noise=root/'noise';comp=root/'comparison'
        for f in [spec,noise,comp]:f.mkdir()
        for f in [spec/'job.json',spec/'COMPLETE.json',noise/'job.json']:f.write_text('{}')
        k=np.array([.15,.16,.18]);nn=dict(k_Mpc_inv=k,kperp_Mpc_inv=k*.6,kpar_Mpc_inv=k*.8,thermal_sigma_mK2=np.full(3,2.),sample_coefficient=np.full(3,.5),no_wedge=np.ones(3,bool),beam_wedge=np.array([True,True,False]),horizon_buffer=np.array([True,False,False]))
        for nu in [130,135]:
            folder=noise/G['tag'](nu);folder.mkdir();np.savez_compressed(folder/'native.npz',**nn)
            G['json_write'](folder/'COMPLETE.json',dict(frequency_MHz=nu,native_sha256=G['sha'](folder/'native.npz')))
            for n in [96,128]:
                for i,m in enumerate(MODELS):
                    pp=G['signal_file'](spec,nu,m,n);pp.parent.mkdir(parents=True,exist_ok=True)
                    val=(i+1)*(1 if n==96 else 1.02)
                    if nu==135 and m=='Iplus':val=0.
                    np.savez_compressed(pp,k=np.array([.05,.35]),delta2=np.full(2,val))
                    G['json_write'](pp.with_suffix('.json'),dict(spectrum_sha256=G['sha'](pp)))
        de=[];pa=[]
        for cut in [.2,.25,.3]:
            for case,num in zip(CASES,[3,2,1]):
                for i,m in enumerate(MODELS):
                    v=i+1;de.append(dict(case=case,kmax_Mpc_inv=cut,model=m,snr_low=np.sqrt(num)*v/(2+.5*v),snr_high=np.sqrt(num)*v*1.02/(2+.5*v*1.02)))
                for i,m in enumerate(MODELS):
                    for j,b in enumerate(MODELS[i+1:],i+1):
                        r=dict(case=case,kmax_Mpc_inv=cut,model_a=m,model_b=b)
                        for label,factor in [('low',1.),('high',1.02)]:
                            delta=(j-i)*factor
                            for lab,v in [('fixed_F',1),('if_a_true',i+1),('if_b_true',j+1)]:r['distance_'+label+'_'+lab]=np.sqrt(num)*delta/(2+.5*v*factor)
                            r['distance_'+label+'_common_low_F_covariance']=np.sqrt(num)*delta/2.5
                        pa.append(r)
        G['csv_write'](comp/'detection_comparison.csv',de);G['csv_write'](comp/'pairwise_comparison.csv',pa)
        args=['forecast_bands.py','--spectra',spec,'--noise',noise,'--comparison',comp,'--out',root/'out','--frequencies','130','135']
        run(*args);run(*args,'--resume')
        import csv
        with (root/'out/detection_all_bands.csv').open() as f:rr=list(csv.DictReader(f))
        assert len(rr)==126
        assert all(float(r['snr_high'])==0 for r in rr if float(r['frequency_MHz'])==135 and r['model']=='Iplus')
    print('PASS: end-to-end multi-band forecast CLI, independent analytic 130 MHz regression, zeros and resume.')
    print('ALL SYNTHETIC TESTS PASSED. Real SKA calculations must pass the separate 130 MHz regressions on your machine.')
if __name__=='__main__':main()
