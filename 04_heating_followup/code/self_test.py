#!/usr/bin/env python3
"""Dependency-free CLI/plan/seed-isolation checks. Never imports py21cmfast."""
import hashlib, importlib.util, json, subprocess, sys, tempfile, py_compile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
import run_scan

def call(*args,ok=True):
    p=subprocess.run([sys.executable,str(ROOT/'run_scan.py'),*map(str,args)],text=True,capture_output=True)
    if ok and p.returncode:raise AssertionError(p.stdout+p.stderr)
    if not ok and not p.returncode:raise AssertionError('Expected rejection')
    return p

def main():
    for path in ROOT.rglob('*.py'):py_compile.compile(str(path),doraise=True)
    for seed in [12345,23456,34567]:
        plan=run_scan.jobs([130,140],4,False,seed)
        assert len(plan)==6
        for tag,c in plan:
            assert c['random_seed']==seed and c['required_version']=='3.4.0'
            assert c['user_params']['BOX_LEN']==150 and c['user_params']['HII_DIM']==128 and c['user_params']['DIM']==384
            assert not c['flag_options']['APPLY_RSDS'] and not c['flag_options']['SUBCELL_RSD']
            assert c['astro_params']['L_X'] in [40.,40.25]
            reference=run_scan.runner.config(c['model'],c['frequency_MHz'],'real',4)
            reference['random_seed']=seed
            if tag=='Mplus_LX40.25':reference['astro_params']['L_X']=40.25
            assert reference==c
    assert run_scan.jobs([130],4,False,12345)[0][1]['redshift']!=run_scan.jobs([140],4,False,12345)[0][1]['redshift']
    try:run_scan.jobs([140],4,True)
    except ValueError:pass
    else:raise AssertionError('Refinement forbidden')
    with tempfile.TemporaryDirectory(prefix='mnras_followup_plan_') as tmp:
        base=Path(tmp)/'runs'
        for seed,freq in [(12345,130),(23456,140),(34567,140)]:
            p=call('--out',base,'--seed',seed,'--frequencies',freq)
            assert 'Plan only' in p.stdout
            root=base/f'seed{seed}'
            plan=json.loads((root/f'plan_{freq}.json').read_text())
            assert len(plan)==3 and all(j['configuration']['random_seed']==seed for j in plan)
            assert not list(root.rglob('cube.npz')) and not (root/'cache').exists()
            for j in plan:
                dest=root/'production'/f'{freq}MHz'/j['tag']/'real'
                assert dest.parents[3]/'cache'/'shared_real'==root/'cache'/'shared_real'
        call('--out',base,'--seed',12345,'--frequencies',130,'--threads',5,ok=False)
        call('--out',base,'--seed',12345,'--frequencies',120,ok=False)
        call('--out',base,'--seed',11111,'--frequencies',130,ok=False)
        call('--out',base,'--seed',12345,'--refine',ok=False)
        # Packaging checks use labelled dummy bytes, never simulation outputs.
        import zipfile
        results=Path(tmp)/'synthetic_results';results.mkdir();(results/'status.json').write_text('{"synthetic_packaging_test": true}')
        cube=base/'seed12345/production/130MHz/F_control/real/cube.npz';cube.parent.mkdir(parents=True);cube.write_bytes(b'SYNTHETIC_PACKAGE_TEST_ONLY')
        cache=base/'seed12345/cache/shared_real/test.bin';cache.parent.mkdir(parents=True);cache.write_bytes(b'SYNTHETIC_CACHE_TEST_ONLY')
        for raw in [False,True]:
            dest=Path(tmp)/('raw.zip' if raw else 'small.zip')
            cmd=[sys.executable,str(ROOT/'package_results.py'),'--runs',str(base),'--results',str(results),'--out',str(dest)]
            if raw:cmd.append('--include-cubes')
            rr=subprocess.run(cmd,text=True,capture_output=True)
            assert rr.returncode==0,rr.stderr
            with zipfile.ZipFile(dest) as z:
                names=z.namelist();assert any(n.endswith('cube.npz') for n in names)==raw
                assert not any('/cache/' in n for n in names)
                assert 'manifest.json' in names

    summary={'CLI_and_seed_plan_checks':'PASS','packaging_checks':'PASS','production_simulations':'NOT_RUN','physics_validation':'NOT_RUN'}
    if all(importlib.util.find_spec(m) for m in ['numpy','pandas','matplotlib']):
        p=subprocess.run([sys.executable,str(ROOT/'synthetic_analysis_test.py')],text=True,capture_output=True)
        if p.returncode:raise AssertionError(p.stdout+p.stderr)
        summary['synthetic_analysis']='PASS';print(p.stdout)
    else:summary['synthetic_analysis']='NOT_RUN_missing_numpy_pandas_or_matplotlib'
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
