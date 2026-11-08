#!/usr/bin/env python3
"""Package derived results and required metadata; raw cubes are opt-in, caches excluded."""
import argparse,hashlib,json,zipfile
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--runs',type=Path,required=True);p.add_argument('--results',type=Path,nargs='+',required=True)
 p.add_argument('--out',type=Path,required=True);p.add_argument('--include-cubes',action='store_true');a=p.parse_args()
 if a.out.exists():p.error('Output exists; choose a new filename')
 if not a.runs.is_dir() or not all(x.is_dir() for x in a.results):p.error('Input directory missing')
 source=Path(__file__).resolve().parent;out=a.out.resolve()
 if out.is_relative_to(source) or out.is_relative_to(a.runs.resolve()) or any(out.is_relative_to(x.resolve()) for x in a.results):p.error('Put output ZIP outside source/runs/results')
 files=[]
 for i,root in enumerate(a.results):
  for f in sorted(root.rglob('*')):
   if f.is_file():files.append((f,f'results/{i}_{root.name}/'+f.relative_to(root).as_posix()))
 for f in sorted(a.runs.rglob('*')):
  rel=f.relative_to(a.runs)
  if not f.is_file() or 'cache' in rel.parts:continue
  if f.suffix in ['.json','.txt'] or (a.include_cubes and f.name=='cube.npz'):files.append((f,'run_metadata/'+rel.as_posix()))
 for f in sorted(source.rglob('*')):
  rel=f.relative_to(source)
  if f.is_file() and '__pycache__' not in rel.parts and f.suffix!='.zip':files.append((f,'code/'+rel.as_posix()))
 manifest={}
 with zipfile.ZipFile(a.out,'x',zipfile.ZIP_DEFLATED) as z:
  for f,name in files:
   z.write(f,name)
   h=hashlib.sha256()
   with f.open('rb') as stream:
    for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
   manifest[name]=h.hexdigest()
  z.writestr('manifest.json',json.dumps(manifest,indent=2))
 print('Saved',a.out,'; cubes included:',a.include_cubes,'; caches excluded')
if __name__=='__main__':main()
