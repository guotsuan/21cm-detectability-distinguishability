"""Heating comparison: exact matched shell ratios and archived same-seed distances."""
from pathlib import Path
import csv, os, json, hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/mnras_fig5_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent
sources=[]
def rows(p):
 sources.append(p)
 with p.open() as f:return list(csv.DictReader(f))
spectra={130:rows(R/'output/followup_review_20261104/returned/results/0_followup_results_130_s12345/spectra.csv'),140:rows(R/'physical_followup_review/spectra.csv')}
old=rows(R/'physical_followup_review/scan_metrics.csv')
sets=[(140,12345,old)]
for folder in ['output/followup_review_20261104/returned/results','output/heating_130_update_20261104/returned/results']:
 for path in sorted((R/folder).iterdir()):
  q=rows(path/'scan_metrics.csv');sets.append((int(q[0]['frequency_MHz']),int(q[0]['seed']),q))
sets=sorted(sets,key=lambda x:x[:2]);assert len(sets)==6
plt.rcParams.update({'font.size':8,'axes.titlesize':8.5,'axes.labelsize':8,'legend.fontsize':7,'pdf.fonttype':42,'lines.linewidth':1.2})
fig,ax=plt.subplots(1,3,figsize=(7.2,2.9),width_ratios=[1,1,1.55],layout='constrained')
colors=['#b55239','#137f8b'];markers=['o','s'];records=[]
for a,nu,panel in zip(ax[:2],[130,140],['(a)','(b)']):
 def spectrum(tag):
  q=sorted([r for r in spectra[nu] if r['tag']==tag],key=lambda r:float(r['k_Mpc_inv']))
  return np.array([float(r['k_Mpc_inv']) for r in q]),np.array([float(r['Delta2_mK2']) for r in q])
 k,f=spectrum('F_control');assert np.all(f>0)
 for lx,c,m in zip([40.,40.25],colors,markers):
  x,p=spectrum(f'Mplus_LX{lx:.2f}');assert np.array_equal(x,k),'Matched shell centres required; no interpolation used.'
  q=p/f;assert np.all(np.isfinite(q))
  a.plot(k,q,m+'-',ms=3,color=c)
  records.extend(dict(frequency_MHz=nu,seed=12345,L_X=lx,k_Mpc_inv=float(ki),ratio=float(qi)) for ki,qi in zip(k,q))
 a.axhline(1,color='.35',ls='--',lw=.8)
 a.set(xlim=(.1,.3),ylim=(0,3.3),xticks=[.1,.2,.3],xlabel=r'$k\ [{\rm Mpc}^{-1}]$',title=f'{panel} {nu} MHz')
 a.grid(alpha=.15);a.set_yticks([0,1,2,3])
ax[0].set_ylabel(r'Power ratio $\Delta^2_{M}/\Delta^2_F$')
ax[1].tick_params(labelleft=False)
distance_rows=[]
for i,(nu,seed,data) in enumerate(sets):
 q=sorted([r for r in data if r['case']=='horizon_buffer'],key=lambda r:float(r['L_X']));assert len(q)==2
 values=[float(r['D_fresh_F']) for r in q]
 ax[2].plot(values,[i,i],color='.75',lw=1.1,zorder=1)
 for v,c,m in zip(values,colors,markers):
  ax[2].scatter(v,i,c=c,marker=m,s=16,zorder=2)
  ax[2].annotate(f'{v:.2f}',(v,i),xytext=(0,5),textcoords='offset points',ha='center',fontsize=6.5)
 distance_rows.append(dict(frequency_MHz=nu,seed=seed,D_original=values[0],D_heating=values[1]))
ax[2].set(yticks=range(6),yticklabels=[f'{nu} / {seed}' for nu,seed,_ in sets],xscale='log',xlim=(.3,28),ylim=(5.5,-.7),xlabel='Horizon-buffer distance',title='(c) Frequency (MHz) / seed')
ax[2].tick_params(axis='y',labelsize=6.5);ax[2].grid(axis='x',alpha=.15)
fig.legend(handles=[Line2D([],[],color=c,marker=m,ms=4,label=l) for c,m,l in zip(colors,markers,[r'Original Mplus ($L_X=40.00$)',r'Heating +0.25 dex ($L_X=40.25$)'])],loc='outside lower center',ncol=2,frameon=False)
for ext in ['pdf','png']:fig.savefig(O/f'physical_heating_comparison.{ext}',dpi=240)
for filename,data in [('shell_ratios.csv',records),('distances.csv',distance_rows)]:
 with (O/filename).open('w') as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
(O/'provenance.json').write_text(json.dumps({'sources':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'shell_matching':'Exact equal shell-centre arrays within each frequency; no interpolation. Ratios calculated before plotting.','spectra_seed':12345,'distance_metric':'D_fresh_F','distance_rows':len(distance_rows),'shell_ratio_rows':len(records),'uncertainty':'None plotted; all points are computed cases.'},indent=2))
(O/'CAPTION.txt').write_text(r'''Fixed physical heating-candidate tests. (a,b) Mplus-to-F coeval power ratios at 130 and 140 MHz for seed 12345, evaluated at matched shell centres; the dashed line denotes equal power. (c) Horizon-buffer distances before and after the same 0.25-dex heating increase for three seeds at each frequency. Each distance uses the same-seed fresh F target and original F covariance, with no extra amplitude fit.''')
print('Saved three-panel figure and numerical provenance.')
