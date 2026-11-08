import csv,json,zipfile,itertools
from pathlib import Path
import numpy as np
from common import sha,dump,ROOT
MODELS=['F','Hminus','Hplus','Mminus','Mplus','Iminus','Iplus']
def package(out):
 rows=[];summary=[]
 for nu in [130,140]:
  folder=out/f'{nu}MHz'
  for seed in [12345,23456,34567]:
   result=folder/f'seed{seed}/results';done=json.loads((result/'complete.json').read_text())
   if done['results_sha256']!=sha(result/'results.npz'):raise ValueError('Results changed')
   cfg=json.loads((result/'configuration.json').read_text());P=np.load(result/'results.npz')['signal_power_Jy2'].reshape(3,7,2,3,4,7)
   for bi,beam in enumerate(cfg['beam_order']):
    for i,j in itertools.combinations(range(7),2):
     for k,kind in enumerate(['coeval','lightcone']):
      for op,pipeline in enumerate(cfg['pipeline_order']):
       for group in range(4):
        for b in range(7):
         x=P[bi,i,k,op,group,b];y=P[bi,j,k,op,group,b];total=x+y;delta=x-y
         rows.append([nu,seed,beam,MODELS[i],MODELS[j],kind,pipeline,group,b,b<6,x,y,delta,delta/total if total>0 else np.nan,abs(delta)>=1e-3*total if total>0 else False])
   summary.append(dict(frequency_MHz=nu,seed=seed,angular_span_deg=cfg['physical_metadata']['direction_cosine_span_deg'],frequency_range_MHz=cfg['physical_metadata']['target_frequency_range_MHz'],regression140=cfg['regression140']))
 with (out/'all_pair_contrasts.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['frequency_MHz','seed','beam','model_a','model_b','kind','pipeline','baseline_group','delay_bin','retained','power_a_Jy2','power_b_Jy2','difference_Jy2','normalized_contrast','above_reporting_floor']);w.writerows(rows)
 dump(out/'summary.json',dict(status='COMPLETE_MATCHED_FIXED_COMOVING_RESPONSE_REVIEW_REQUIRED',cases=summary,models=MODELS,pairs=21,seeds=3,frequencies=[130,140],interpretation='Same comoving window, different angular support. Equal delay support is not equal k_parallel. No independent-band combination or calibrated exclusion.',required_review=['130/140 SEFD selections and 140 noise/foreground/power regression','Thermal Monte Carlo covariance agreement and non-diagonal correlations','Foreground residual and signal cross-term biases','All21 pair contrasts on common reporting support','Finite-window and interpolation limits before a frequency-ranking claim']))
 target=out/'matched_response_review.zip'
 with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
  for name in ['all_pair_contrasts.csv','summary.json']:z.write(out/name,name)
  for nu in [130,140]:
   folder=out/f'{nu}MHz'
   files=[folder/'geometry.npz',folder/'prepared_complete.json']
   for d in ['noise','noise_validation']:files+=list((folder/d).glob('*'))
   files +=[folder/'beams/configuration.json',folder/'beams/complete.json']
   for seed in [12345,23456,34567]:
    files+=[folder/f'seed{seed}/prepared/metadata.json'];files+=list((folder/f'seed{seed}/results').glob('*'))
   for p in files:
    if p.is_file():z.write(p,p.relative_to(out))
  for p in ROOT.glob('*.py'):z.write(p,'code/'+p.name)
 print('Send:',target.resolve())
