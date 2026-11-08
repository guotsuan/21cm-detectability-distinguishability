#!/usr/bin/env python3
"""Replay selected paper figures from compact inputs in a fresh output directory."""
import argparse, json, shutil, subprocess, sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();out=a.out.resolve()
    if out.exists():p.error('Output directory must not already exist')
    root=Path(__file__).resolve().parent
    shutil.copytree(root/'paper_snapshot',out/'work')
    work=out/'work';dest=out/'figures';dest.mkdir()
    jobs=[
      ('output/fig5_compact_20261104/plot_fig5.py', [], 'output/fig5_compact_20261104/main_F_Mplus.pdf','main_F_Mplus.pdf'),
      ('output/fig3_preview_20261104/plot_preview.py', [], 'output/fig3_preview_20261104/fig3_preview.pdf','seven_model_tracked.pdf'),
      ('output/submission_edit_20261104/figure/plot_heating.py', [], 'output/submission_edit_20261104/figure/physical_heating_comparison.pdf','physical_heating140.pdf'),
      ('matched_response_review/figure_source/plot_matched.py', ['--figure','3','--out',str(work/'matched_figures')], 'matched_figures/fig03.pdf','seven_model_tracked_all_pairs.pdf'),
      ('matched_response_review/figure_source/plot_matched.py', ['--figure','18','--out',str(work/'matched_figures')], 'matched_figures/fig18.pdf','matched_noise.pdf')]
    status=[]
    for script,args,product,name in jobs:
        proc=subprocess.run([sys.executable,str(work/script),*args],cwd=work,capture_output=True,text=True)
        (out/(name+'.log')).write_text(proc.stdout+proc.stderr)
        status.append({'figure':name,'returncode':proc.returncode})
        (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
        if proc.returncode:raise SystemExit(proc.returncode)
        shutil.copy2(work/product,dest/name)
    print(json.dumps(status,indent=2))
if __name__=='__main__':main()
