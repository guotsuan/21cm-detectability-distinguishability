#!/usr/bin/env python3
"""Synthetic end-to-end analysis and failure gates; does not run 21cmFAST."""
import argparse,json,tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import analyze_scan as analysis
from common import sha,fit
from run_scan import jobs

def main():
    assert len(jobs([140],4))==3 and len(jobs([130,140],4))==6
    try: jobs([140],4,True)
    except ValueError: pass
    else: raise AssertionError('Refinement must be forbidden')
    a,d,ds=fit(np.array([2.,4.]),np.array([1.,2.]),np.ones(2))
    assert a==2 and abs(ds)<1e-12 and abs(d-np.sqrt(5))<1e-12
    with tempfile.TemporaryDirectory(prefix='eor_followup_test_') as tmp:
        root=Path(tmp);runs=root/'runs';forecast=root/'forecast'
        rng=np.random.default_rng(17);cube=rng.normal(size=(128,128,128)).astype('float32')
        density=rng.normal(size=cube.shape).astype('float32')
        k,p,*_=analysis.fft_power(cube,150);km,pm,de,count,sums=analysis.spherical(k,p,analysis.EDGES)
        sample=np.linspace(.105,.295,30)
        mapped=analysis.interp(pd.DataFrame(dict(k_Mpc_inv=km,Delta2_mK2=de)),sample)
        dest=forecast/'forecasts/130MHz';dest.mkdir(parents=True)
        sigma=np.ones_like(sample)*2
        np.savez(dest/'native_N128.npz',k_Mpc_inv=sample,F_signal_mK2=mapped,Mplus_signal_mK2=mapped*4,F_total_sigma_mK2=sigma,thermal_sigma_mK2=sigma,sample_coefficient=np.zeros_like(sample),no_wedge=np.ones(30,bool),beam_wedge=np.arange(30)%2==0,horizon_buffer=np.arange(30)%3==0)
        for (tag,c),scale in zip(jobs([130],4),[1,2,1]):
            folder=runs/'seed12345/production/130MHz'/tag/'real';folder.mkdir(parents=True)
            (folder/'requested.json').write_text(json.dumps(c))
            (folder/'resolved.json').write_text(json.dumps(dict(astro_params=c['astro_params'],flag_options=c['flag_options'])))
            (folder/'accepted_parameters.json').write_text(json.dumps({key:dict(log10_input=c['astro_params'][key],resolved=c['astro_params'][key]) for key in ['L_X','M_TURN']}))
            (folder/'environment.json').write_text(json.dumps(dict(packages={'21cmFAST':'3.4.0'},loaded_modules={'synthetic':True})))
            np.savez(folder/'cube.npz',brightness_temp=cube*scale,density=density,xH_box=np.ones(cube.shape,dtype='float32'),Tk_box=np.ones(cube.shape,dtype='float32'),Ts_box=np.ones(cube.shape,dtype='float32'),box_len_cMpc=150,random_seed=12345,redshift=c['redshift'],coordinate_space='real')
            (folder/'COMPLETE.json').write_text(json.dumps(dict(cube_sha256=sha(folder/'cube.npz'),elapsed_seconds=0)))
        args=SimpleNamespace(runs=runs,out=root/'success',forecast_root=forecast,frequencies=[130],seed=12345,regression_rtol=3e-5,refine=False)
        analysis.main(args)
        t=pd.read_csv(args.out/'scan_metrics.csv')
        assert len(t)==6 and np.all(t[t.L_X==40.25].D_fixed_F<1e-12)
        assert np.all(t[t.L_X==40].D_fixed_F>0)
        assert np.allclose(t.D_fixed_F,t.D_fresh_F)
        assert (args.out/'fixed_candidate_transfer.csv').exists()
        try:analysis.main(args)
        except ValueError as e:assert 'Output is not empty' in str(e)
        else:raise AssertionError('Must reject output overwrite')
        # A new seed must use matched fresh F, without forcing seed12345 baseline equality.
        import shutil
        other=runs/'seed23456'
        shutil.copytree(runs/'seed12345',other)
        for dest2 in (other/'production/130MHz').glob('*/real'):
            cfg=json.loads((dest2/'requested.json').read_text());cfg['random_seed']=23456
            (dest2/'requested.json').write_text(json.dumps(cfg))
            with np.load(dest2/'cube.npz') as zz:data2={k:zz[k] for k in zz.files}
            data2['random_seed']=23456
            # Reverse the preferred input between archived and fresh target means.
            scale={'F_control':2.0,'Mplus_LX40.00':0.95,'Mplus_LX40.25':1.0}[dest2.parent.name]
            data2['brightness_temp']*=scale
            np.savez(dest2/'cube.npz',**data2)
            (dest2/'COMPLETE.json').write_text(json.dumps(dict(cube_sha256=sha(dest2/'cube.npz'),elapsed_seconds=0)))
        args.seed=23456;args.out=root/'other_seed'
        analysis.main(args)
        tt=pd.read_csv(args.out/'scan_metrics.csv')
        assert np.all(tt[tt.L_X==40.25].D_fresh_F>0)
        assert np.all(tt[tt.L_X==40.25].D_fixed_F<1e-10)
        bb=pd.read_csv(args.out/'best_sampled.csv')
        assert np.all(bb.best_sampled_L_X==40.0)
        assert set(bb.selection_metric)=={'D_fresh_F'}
        assert np.allclose(bb.distance_ratio_fresh,1.0)
        transfer=pd.read_csv(args.out/'fixed_candidate_transfer.csv')
        assert np.all(transfer.distance_ratio_fresh>1)
        assert np.all(transfer.distance_ratio_archived_secondary<1e-10)
        assert np.allclose(tt.D_archived_F_secondary,tt.D_fixed_F)
        assert json.loads((args.out/'status.json').read_text())['archived_control_regression']=='not_applicable_different_realization'
        args.seed=12345
        # A changed native baseline must be reported and block metrics.
        nativefile=dest/'native_N128.npz'
        with np.load(nativefile) as z:native_data={k:z[k] for k in z.files}
        bad=dict(native_data);bad['F_signal_mK2']=native_data['F_signal_mK2']*1.1
        np.savez(nativefile,**bad);args.out=root/'regression_failure'
        try:analysis.main(args)
        except ValueError as e:assert 'Baseline regression failed' in str(e)
        else:raise AssertionError('Must reject baseline regression')
        assert (args.out/'control_regression.csv').exists()
        assert not (args.out/'scan_metrics.csv').exists()
        np.savez(nativefile,**native_data)
        # Mutate density while updating file checksum: density gate must still catch it.
        with np.load(folder/'cube.npz') as z:data={k:z[k] for k in z.files}
        data['density'][0,0,0]+=1
        np.savez(folder/'cube.npz',**data)
        (folder/'COMPLETE.json').write_text(json.dumps(dict(cube_sha256=sha(folder/'cube.npz'),elapsed_seconds=0)))
        args.out=root/'density_failure'
        try:analysis.main(args)
        except ValueError as e:assert 'Fresh density fields differ' in str(e)
        else:raise AssertionError('Must reject altered density')
        assert (args.out/'density_diagnostics.csv').exists()
        assert not (args.out/'scan_metrics.csv').exists()
        # A corrupted stored hash must fail independently of numerical comparisons.
        (folder/'COMPLETE.json').write_text(json.dumps(dict(cube_sha256='incorrect',elapsed_seconds=0)))
        args.out=root/'hash_failure'
        try:analysis.main(args)
        except ValueError as e:assert 'checksum' in str(e)
        else:raise AssertionError('Must reject corrupted cube')
    print('PASS: stage plans, amplitude metric, full synthetic analysis, transfer output, density diagnostics, checksum and output protection. No 21cmFAST simulation was executed.')
if __name__=='__main__':main()
