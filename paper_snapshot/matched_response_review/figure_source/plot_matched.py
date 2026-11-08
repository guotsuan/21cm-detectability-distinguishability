#!/usr/bin/env python3
"""Figures3 and18: reviewed two-frequency, seven-model response data."""
from pathlib import Path
import argparse,itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]
COLORS=['#0072B2','#D55E00']
def plot(number,out):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if number==3:
  with np.load(ROOT/'figure_source/physical_powers.npz') as z:P=z['power'];models=z['models']
  pairs=list(itertools.combinations(range(7),2));fig,axes=plt.subplots(1,2,figsize=(8.4,6.8),sharey=True,layout='constrained')
  for ax,op,title in zip(axes,[0,2],['Fixed catalogue','Quadratic fitting']):
   for fi,nu in enumerate([130,140]):
    for si,seed in enumerate([12345,23456,34567]):
     vals=[]
     for i,j in pairs:
      a=P[fi,si,0,i,1,op,:,0];b=P[fi,si,0,j,1,op,:,0];vals.append(np.median(abs((a-b)/(a+b))))
     ax.scatter(vals,np.arange(21)+(fi-.5)*.28+(si-1)*.065,c=COLORS[fi],marker=['o','s','^'][si],s=16,alpha=.85)
   ax.set(xlim=(-.02,1.02),xlabel='Median absolute normalized contrast',title=title);ax.grid(axis='x',alpha=.2)
   for y in range(21):ax.axhline(y,c='grey',alpha=.09,lw=.6)
  axes[0].set_yticks(range(21),[models[i]+'–'+models[j] for i,j in pairs],fontsize=8);axes[0].invert_yaxis()
  handles=[Line2D([],[],color=c,marker='o',linestyle='',label=f'{nu} MHz') for c,nu in zip(COLORS,[130,140])]+[Line2D([],[],color='.4',marker=m,linestyle='',label=str(seed)) for m,seed in zip(['o','s','^'],[12345,23456,34567])]
  axes[1].legend(handles=handles,loc='lower right',fontsize=7)
 else:
  d=np.load(ROOT/'figure_source/fig18.npz');fig,axes=plt.subplots(2,3,figsize=(10,6),layout='constrained')
  for fi,nu in enumerate([130,140]):
   for t in range(3):axes[fi,0].plot(d['frequency_hz'][fi]/1e6,d['sefd'][fi,t],label=f'Time {t+1}')
   axes[fi,0].set(xlabel='Frequency (MHz)',ylabel='Station-equivalent SEFD I (Jy)',title=f'{nu} MHz');axes[fi,0].legend(fontsize=7)
   C=d['covariance'][fi,2,:7,:7];R=C/np.sqrt(np.outer(np.diag(C),np.diag(C)));im=axes[fi,1].imshow(R,vmin=-1,vmax=1,cmap='RdBu_r',origin='lower');axes[fi,1].set(xlabel='Delay-bin index',ylabel='Delay-bin index',title='Quadratic: noise correlation');fig.colorbar(im,ax=axes[fi,1],shrink=.8)
   for j,label in enumerate(['Fixed','Linear','Quadratic']):axes[fi,2].plot(np.arange(28),d['variance_ratio'][fi,j],'.-',ms=3,lw=.5,label=label)
   axes[fi,2].axhline(1,c='k',ls='--',lw=.6);axes[fi,2].set(xlabel='Group × delay output',ylabel='MC / analytic variance',ylim=(.89,1.11),title='2048 independent split pairs');axes[fi,2].legend(fontsize=7)
 for ext in ['pdf','png']:fig.savefig(out/f'fig{number:02d}.{ext}',dpi=180)
 plt.close(fig)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--figure',type=int,choices=[3,18],required=True);p.add_argument('--out',type=Path,default=ROOT/'regenerated');a=p.parse_args();plot(a.figure,a.out)
