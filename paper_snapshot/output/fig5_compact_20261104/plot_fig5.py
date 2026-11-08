"""Compact Fig. 5; cosmetic x offsets only, original archived values unchanged."""
from pathlib import Path
import csv, hashlib, json, os
os.environ.setdefault('MPLCONFIGDIR','/tmp/mnras_fig5_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
E=ROOT/'focus_review/evidence'
for name, rec in json.loads((E/'provenance.json').read_text()).items():
    assert hashlib.sha256((E/name).read_bytes()).hexdigest()==rec['sha256']
def read(p):
    with p.open() as f: return list(csv.DictReader(f))
def num(r,k): return float(r[k])
cov=[r for r in read(E/'covariance_comparison.csv') if r['N']=='128' and num(r,'kmax_Mpc_inv')==.3 and r['case']=='horizon_buffer' and num(r,'frequency_MHz') in (120,130,140)]
seed=[r for r in read(E/'seed_metrics.csv') if num(r,'kmax_Mpc_inv')==.3 and r['case']=='horizon_buffer' and r['covariance']=='fixed_seed12345']
box=[r for r in read(E/'box_metrics.csv') if num(r,'kmin')==.1 and num(r,'kmax')==.3 and r['case']=='horizon_buffer' and r['covariance']=='fixed_original_F_or_model']
def choose(rows,stat):
    return [r for r in rows if (r['kind']=='pair' and r['model_a']=='F' and r['model_b']=='Mplus') if stat=='F-Mplus'] if stat=='F-Mplus' else [r for r in rows if r['kind']=='detection' and r['model_a']==stat]
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'legend.fontsize':6.5,'lines.linewidth':1.2,'lines.markersize':3.5,'pdf.fonttype':42})
fig,ax=plt.subplots(1,3,figsize=(7.2,2.55),layout='constrained')
colors={'F':'#2468a2','Mplus':'#d47719','F-Mplus':'#71459a'}
records=[]
def line(axis,rows,key,shift,label,color,style='o-',factor=1):
    rows=sorted(rows,key=lambda r:num(r,'frequency_MHz'))
    assert len(rows)==3
    x=[num(r,'frequency_MHz') for r in rows]; y=[factor*num(r,key) for r in rows]
    axis.plot([v+shift for v in x],y,style,color=color,label=label)
    records.extend(dict(panel=list(ax).index(axis),series=label,true_frequency_MHz=a,plotted_frequency_MHz=a+shift,value=b) for a,b in zip(x,y))
for stat in ['F','Mplus','F-Mplus']:
    q=[r for r in choose(cov,stat) if r['covariance_model']==('F' if stat=='F-Mplus' else stat)]
    pub=[r for r in q if r['convention']=='public']
    axis=ax[1] if stat=='F-Mplus' else ax[0]
    shift={'F':-.3,'Mplus':.3,'F-Mplus':-.6}[stat]
    line(axis,pub,'metric',shift,stat+' public',colors[stat])
    for nu in [120,130,140]:
        v=[num(r,'metric') for r in choose(seed,stat) if num(r,'frequency_MHz')==nu]
        assert len(v)==3
        axis.vlines(nu+shift,min(v),max(v),color=colors[stat],alpha=.45,linewidth=3)
    b=[r for r in box if r['statistic']==stat]
    axis.scatter([num(r,'frequency_MHz')+shift for r in b],[num(r,'large_box_metric') for r in b],marker='D',s=22,facecolors='none',edgecolors=colors[stat],linewidths=.8,zorder=5)
    if stat in ['F','F-Mplus']:
        line(ax[2],pub,'retained_information',-.6 if stat=='F' else 0,'F detection' if stat=='F' else 'F–Mplus public',colors[stat],factor=100)
    if stat=='F-Mplus':
        for convention,offset,style,label in [('joint_rho0',-.2,'^--',r'Joint, $\rho=0$'),('joint_rho1',.2,'v:',r'Joint, $\rho=1$')]:
            rows=[r for r in q if r['convention']==convention]
            line(ax[1],rows,'metric',offset,label,'#555555',style)
            if convention=='joint_rho1': line(ax[2],rows,'retained_information',.6,r'F–Mplus joint, $\rho=1$',colors[stat],'v:',100)
shape=[r for r in read(ROOT/'interpretation_review/results/shape/shape_summary.csv') if r['case']=='horizon_buffer']
line(ax[1],shape,'D_shape',.6,'Amplitude-profiled','#16805e','s-.')
for i,a in enumerate(ax):
    a.grid(alpha=.18);a.set_xlim(118,142);a.set_xticks([120,130,140]);a.set_xlabel('Band centre (MHz)')
    a.set_title(['(a) Detection','(b) F–Mplus separation','(c) Mode retention'][i],fontsize=9,loc='left')
    a.legend(loc='upper center',bbox_to_anchor=(.5,1),frameon=False,handlelength=1.8,labelspacing=.25)
ax[0].set_ylabel('Conditional SNR');ax[0].set_ylim(9,18)
ax[1].set_ylabel('Weighted distance');ax[1].set_ylim(0,28)
ax[2].set_ylabel('Retained squared metric (%)');ax[2].set_ylim(0,4.8)
fig.savefig(OUT/'main_F_Mplus.pdf');fig.savefig(OUT/'main_F_Mplus.png',dpi=220)
with (OUT/'plotted_values.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
print('Saved compact Fig. 5 with',len(records),'unchanged curve values.')
