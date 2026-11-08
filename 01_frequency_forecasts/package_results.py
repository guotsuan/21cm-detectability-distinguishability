#!/usr/bin/env python3
"""Package small derived results for review; never include simulation cubes."""
import argparse,json,zipfile
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path('frequency_analysis'));p.add_argument('--stage',choices=['pilot','full'],default='full');p.add_argument('--out',type=Path,default=Path('frequency_review.zip'));a=p.parse_args()
    if a.out.exists():raise FileExistsError('Choose a new archive filename')
    freqs=[130] if a.stage=='pilot' else list(range(55,196,5))
    for nu in freqs:
        for part in ['noise','forecasts']:
            f=a.root/part/f'{nu}MHz/COMPLETE.json'
            if not f.exists():raise FileNotFoundError(f)
    if not (a.root/'spectra/COMPLETE.json').exists():raise FileNotFoundError('Spectra incomplete')
    chosen=[]
    for pth in sorted(a.root.rglob('*')):
        if not pth.is_file():continue
        if pth.suffix in {'.csv','.json','.png'} and 'snapshots' not in pth.parts:chosen.append(pth)
        elif pth.suffix=='.npz' and ('130MHz' in pth.parts or pth.name=='spectral_response.npz'):chosen.append(pth)
    with zipfile.ZipFile(a.out,'w',zipfile.ZIP_DEFLATED) as z:
        for f in chosen:z.write(f,f.relative_to(a.root.parent))
    print('Created:',a.out.resolve(),'files:',len(chosen))
if __name__=='__main__':main()
