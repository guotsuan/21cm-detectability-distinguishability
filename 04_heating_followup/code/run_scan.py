#!/usr/bin/env python3
"""Seed-isolated fresh controls and fixed L_X=40.25 transfer; no parameter refinement."""
import argparse,sys,os,subprocess,shutil,inspect,importlib.metadata
from pathlib import Path
import json, hashlib
def require(ok,msg):
    if not ok: raise ValueError(msg)
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import run_physics as runner

def jobs(frequencies,threads,refine=False,seed=12345):
    require(not refine,"Parameter refinement is forbidden in this fixed-candidate protocol")
    out=[]
    for nu in frequencies:
        c=runner.config('F',nu,'real',threads);c['random_seed']=seed;out.append(('F_control',c))
        for lx in ([40.,40.20,40.25,40.30] if refine else [40.,40.25]):
            c=runner.config('Mplus',nu,'real',threads);c['random_seed']=seed;c['astro_params']['L_X']=lx
            out.append((f'Mplus_LX{lx:.2f}',c))
    return out

def check_environment(out):
    import py21cmfast as p21
    require(p21.__version__=='3.4.0',f'Requires py21cmfast 3.4.0, found {p21.__version__}; no automatic upgrade')
    source=inspect.getsource(p21.AstroParams)
    (out/'installed_AstroParams_source.txt').write_text(source)
    (out/'environment.json').write_text(json.dumps(dict(python=sys.version,py21cmfast=p21.__version__,packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions() if d.metadata.get('Name')},AstroParams_source_sha256=hashlib.sha256(source.encode()).hexdigest()),indent=2))
    print('Pinned version verified; installed AstroParams source archived for L_X convention audit.',flush=True)

def audit_accepted(c,dest):
    import numpy as np
    import py21cmfast as p21
    expected=dict(p21.AstroParams(c['astro_params']).defining_dict)
    resolved=json.loads((dest/'resolved.json').read_text())['astro_params']
    for key in ['L_X','M_TURN']:
        require(np.isclose(float(resolved[key]),float(expected[key]),rtol=1e-12),f'Resolved parameter differs: {key}')
    (dest/'accepted_parameters.json').write_text(json.dumps({k:dict(log10_input=c['astro_params'][k],resolved=resolved[k]) for k in ['L_X','M_TURN']},indent=2))

def main(a):
    require(a.threads>=1,'threads must be positive');require(a.seed in [12345,23456,34567],'Unsupported seed');a.out=a.out/f'seed{a.seed}';a.out.mkdir(parents=True,exist_ok=True)
    plan=jobs(a.frequencies,a.threads,False,a.seed)
    lock=a.out/"followup_root.json"
    identity=dict(protocol="fixed_candidate_seed_isolated_v2",threads=a.threads,seed=a.seed)
    if lock.exists():require(json.loads(lock.read_text())==identity,"Root configuration differs")
    else:
        require(not any(a.out.iterdir()),"Use a new empty followup root; old results cannot be reused")
        lock.write_text(json.dumps(identity,indent=2))
    print(f'{len(plan)} snapshot jobs; fixed L=150 cMpc, HII_DIM=128, DIM=384, seed={a.seed}, real space.',flush=True)
    for tag,c in plan:print(c['frequency_MHz'],tag,'M_TURN=',c['astro_params']['M_TURN'],'L_X=',c['astro_params']['L_X'],flush=True)
    (a.out/('plan_'+('refine_' if a.refine else '')+'_'.join(map(str,a.frequencies))+'.json')).write_text(json.dumps([dict(tag=tag,configuration=c) for tag,c in plan],indent=2))
    if not a.execute:
        print('Plan only. Add --execute for simulations; Only this new followup root is resumable.');return
    check_environment(a.out)
    for tag,c in plan:
        dest=a.out/'production'/f"{c['frequency_MHz']}MHz"/tag/'real'
        env=os.environ.copy();env.update(OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',json.dumps(c),'--out',str(dest)],check=True,env=env)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--frequencies',nargs='+',type=int,choices=[130,140],default=[140]);p.add_argument('--threads',type=int,default=4);p.add_argument('--execute',action='store_true');p.add_argument('--seed',type=int,choices=[12345,23456,34567],default=12345);p.set_defaults(refine=False);p.add_argument('--worker',help=argparse.SUPPRESS)
    a=p.parse_args()
    if a.worker:
        c=json.loads(a.worker);runner.run_one(c,a.out);audit_accepted(c,a.out)
    else:main(a)
