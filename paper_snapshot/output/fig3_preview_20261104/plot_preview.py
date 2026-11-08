from pathlib import Path
import os,csv,hashlib,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/mnras_fig3_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
source=ROOT/'matched_response_review/figure_source/physical_powers.npz'
with np.load(source) as z:
 P=z['power'];models=z['models'].astype(str)
assert models[0]=='F'
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'pdf.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(7.2,2.8),sharex=True,sharey=True,layout='constrained')
colors=['#0072B2','#D55E00'];records=[]
for ax,op,title in zip(axes,[0,2],['(a) Fixed catalogue','(b) Quadratic fitting']):
 for fi,nu in enumerate([130,140]):
  for row,j in enumerate(range(1,7)):
   vals=[]
   for si,seed in enumerate([12345,23456,34567]):
    a=P[fi,si,0,0,1,op,:,0];b=P[fi,si,0,j,1,op,:,0]
    assert np.all(a+b>0)
    vals.append(float(np.median(np.abs((a-b)/(a+b)))))
   lo,mid,hi=min(vals),float(np.median(vals)),max(vals)
   y=row+(fi-.5)*.28
   ax.errorbar(mid,y,xerr=[[mid-lo],[hi-mid]],fmt=['o','s'][fi],color=colors[fi],markersize=4,capsize=2,elinewidth=1.3,markeredgewidth=.8)
   records.append(dict(treatment=title[4:],frequency_MHz=nu,pair='F-'+models[j],seed12345=vals[0],seed23456=vals[1],seed34567=vals[2],minimum=lo,median=mid,maximum=hi))
 ax.set_xlim(-.02,1.02);ax.set_xticks([0,.2,.4,.6,.8,1]);ax.set_title(title,loc='left');ax.grid(axis='x',alpha=.17)
 ax.set_xlabel(r'Normalized contrast $|P_F-P_m|/(P_F+P_m)$')
 ax.spines[['top','right']].set_visible(False)
 for y in range(6):ax.axhline(y,color='.7',lw=.5,alpha=.2,zorder=0)
axes[0].set_yticks(range(6),['F–'+m for m in models[1:]])
axes[0].set_ylim(5.5,-.65)
fig.legend(handles=[Line2D([],[],color=c,marker=m,linewidth=1.3,label=f'{nu} MHz') for c,m,nu in zip(colors,['o','s'],[130,140])],loc='outside lower center',ncol=2,frameon=False)
for ext in ['png','pdf']:fig.savefig(OUT/f'fig3_preview.{ext}',dpi=220)
with (OUT/'summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
(OUT/'caption.txt').write_text('Normalized lightcone power contrasts between the fiducial model F and six parameter variations in the lowest absolute-delay interval, [0, 0.0625) microseconds, with dipole truth. Panels compare fixed-catalogue subtraction and quadratic fitting. For each seed, the contrast is first summarized by its median over four baseline groups; points and horizontal bars then give the median and full range across three seeds. Colours and symbols distinguish 130 and 140 MHz.\n')
(OUT/'provenance.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'power_shape':list(P.shape),'models':models.tolist(),'rows':len(records),'aggregation':'median over four baseline groups separately per seed, then median and min-max across three seeds'},indent=2))
print('Created preview, 24 summary points; original manuscript unchanged.')
