#!/usr/bin/env python3
"""Conditional foreground/signal injection pilot; not a calibrated SKA forecast."""
import argparse,csv,hashlib,json,shutil,zipfile,platform
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
EDGES=np.array([0,.0625,.25,.5,1,2,4,8.01])*1e-6
MASKS=[(abs(np.fft.fftfreq(128,62500))>=a)&(abs(np.fft.fftfreq(128,62500))<b) for a,b in zip(EDGES[:-1],EDGES[1:])]
B=np.array([m/m.sum() for m in MASKS]);WINDOW=np.hanning(128);WINDOW/=np.sqrt(np.mean(WINDOW**2))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def base(v):return np.fft.fft((v-v.mean(-1,keepdims=True))*WINDOW,axis=-1,norm='ortho')
def grouped(x):
 q=np.stack([x[...,m].mean(-1) for m in MASKS],axis=-1)
 return q.reshape(*q.shape[:-3],3,4,4,7).mean(axis=(-4,-2)).reshape(*q.shape[:-3],28)
def point_response(uvw,nu,dirs):
 c=299792458.;dec=np.deg2rad(-26.82472208);l,m,n=dirs.T
 d=uvw@(dirs-np.array([0,0,1.])).T
 du=-np.sin(dec)*uvw[:,1]+np.cos(dec)*uvw[:,2]
 rate=7.272205217e-5*(du[:,None]*l+np.sin(dec)*uvw[:,0,None]*m-np.cos(dec)*uvw[:,0,None]*(n-1))
 beam=np.exp(-8*np.log(2)*(l*l+m*m)/np.sin(np.deg2rad(4)*140e6/nu)**2)
 return beam[None,:]*np.exp(2j*np.pi*nu*d/c)*np.sinc(62500*d/c)*np.sinc(10*nu*rate/c)
def make_catalogue(seed):
 rng=np.random.default_rng(seed);n=24
 # Finite illustrative catalogue; no measured source counts or diffuse emission.
 lm=rng.uniform(-np.deg2rad(4),np.deg2rad(4),(n,2))
 flux=10**rng.uniform(-1,1.5,n);alpha=rng.normal(-.8,.15,n);beta=rng.normal(0,.1,n)
 err=rng.normal(0,.01,n);da=rng.normal(0,.02,n)
 return np.c_[lm,flux,alpha,beta,err,da]
HEADER='l,m,flux140_Jy,alpha,curvature,model_fractional_flux_error,model_alpha_error'
def observe(cat,uvw,f):
 l,m=cat[:,:2].T;assert np.all(l*l+m*m<1) and np.all(cat[:,2]>0) and np.all(cat[:,5]>-1)
 dirs=np.c_[l,m,np.sqrt(1-l*l-m*m)];true=np.empty((3,16,128),complex);model=true.copy()
 for k,nu in enumerate(f):
  x=np.log(nu/140e6);K=point_response(uvw.reshape(-1,3),nu,dirs)
  sf=cat[:,2]*np.exp(cat[:,3]*x+cat[:,4]*x*x)
  sm=cat[:,2]*(1+cat[:,5])*np.exp((cat[:,3]+cat[:,6])*x)
  true[:,:,k]=(K@sf).reshape(3,16);model[:,:,k]=(K@sm).reshape(3,16)
 return true,model
class Operator:
 def __init__(self,var,filter_on):
  x=np.linspace(-1,1,128);self.X=np.polynomial.legendre.legvander(x,2);self.A=[]
  for v in var:
   w=1/v;self.A.append(np.linalg.solve(self.X.T@(w[:,None]*self.X),self.X.T*w))
  self.A=np.array(self.A);self.on=filter_on
  self.L=np.array([base((np.eye(128)-self.X@a).T).T if filter_on else base(np.eye(128)).T for a in self.A])
 def apply(self,v):
  if self.on:
   coeff=np.einsum('...tbf,tkf->...tbk',v,self.A,optimize=True)
   v=v-np.einsum('...tbk,fk->...tbf',coeff,self.X,optimize=True)
  return base(v)
def theory(var,L,signal):
 N=np.zeros((28,28));M=N.copy()
 for t in range(3):
  C=(L[t]*var[t])@L[t].conj().T
  for g in range(4):
   sl=slice(7*g,7*g+7);s=signal[t,g*4:(g+1)*4]
   N[sl,sl]+=4*(B@(.5*abs(C)**2)@B.T)/144
   M[sl,sl]+=B@np.real(C*(s.conj().T@s))@B.T/144
 return N,M

def run(a):
 a.out.mkdir(parents=True,exist_ok=False)
 inp=ROOT/'inputs/frozen_inputs.npz';pr=json.loads((ROOT/'inputs/provenance.json').read_text());assert sha(inp)==pr['frozen_input_sha256']
 r=np.load(inp);v=r['signal_vis_Jy'];var=r['complex_variance_I_Jy2'];f=r['frequencies_hz'];assert v.shape==(3,16,128) and var.shape==(3,128)
 assert np.isfinite(var).all() and (var>0).all()
 cat=make_catalogue(a.seed) if a.catalogue is None else np.atleast_2d(np.genfromtxt(a.catalogue,delimiter=',',skip_header=1))
 assert cat.ndim==2 and cat.shape[1]==7 and np.isfinite(cat).all()
 np.savetxt(a.out/'catalogue.csv',cat,delimiter=',',header=HEADER,comments='')
 fg,model=observe(cat,r['uvw_m'],f);residuals={'unsubtracted':fg,'exact_subtraction':np.zeros_like(fg),'imperfect_subtraction':fg-model}
 rows=[];reports={};arrays={'foreground_vis':fg,'model_vis':model,'signal_vis':v,'frequency_hz':f,'complex_variance_I':var};rawsignal=grouped(abs(base(v))**2)
 # Exclude very small baseline signal outputs from relative diagnostics; keep absolute values.
 supported=(rawsignal>rawsignal.max()*1e-10) & np.tile(np.arange(7)<6,4)
 for name,on in [('baseline',False),('quadratic',True)]:
  op=Operator(var,on);s=op.apply(v);ps=grouped(abs(s)**2)
  test=np.random.default_rng(1).normal(size=v.shape)+1j*np.random.default_rng(2).normal(size=v.shape)
  np.testing.assert_allclose(op.apply(test),np.einsum('tij,tbj->tbi',op.L,test),rtol=1e-11,atol=1e-12)
  N,MS=theory(var,op.L,s)
  if not on:
   np.testing.assert_allclose(N,r['reference_noise_covariance_I'],rtol=1e-11,atol=1e-20)
   np.testing.assert_allclose(MS,r['reference_mixed_covariance_I'],rtol=1e-11,atol=1e-25)
  arrays[name+'_operator']=op.L;arrays[name+'_noise_covariance']=N;arrays[name+'_signal_power']=ps
  signals={};means={};covs={};samples={key:[] for key in residuals}
  for key,res in residuals.items():
   h=op.apply(res);tot=s+h;mu=grouped(abs(tot)**2);prr=grouped(abs(h)**2);cross=grouped(2*np.real(s*h.conj()))
   np.testing.assert_allclose(mu,ps+prr+cross,rtol=1e-10,atol=1e-18)
   nn,mm=theory(var,op.L,tot);signals[key]=tot;means[key]=mu;covs[key]=nn+mm
   arrays[name+'_'+key+'_conditional_covariance']=nn+mm
   arrays[name+'_'+key+'_mean']=mu
   # Paired +/- injections remove the foreground cross term; finite-precision residual reported.
   recovery=[]
   for alpha in [1.,10.,100.]:
    rec=(grouped(abs(h+alpha*s)**2)+grouped(abs(h-alpha*s)**2)-2*prr)/(2*alpha**2)
    recovery.append(float(np.max(abs(rec[supported]-ps[supported])/rawsignal[supported])))
   for i in range(28):
    rows.append([name,key,i//7,i%7,bool(supported[i]),rawsignal[i],ps[i],ps[i]/rawsignal[i] if supported[i] else np.nan,prr[i],cross[i],mu[i]-ps[i],np.sqrt(covs[key][i,i]),(mu[i]-ps[i])/np.sqrt(covs[key][i,i])])
   reports[name+'_'+key]={'symmetric_injection_max_error_relative_to_baseline_signal_at_amplitudes_1_10_100':recovery}
  rng=np.random.default_rng(a.seed+1)
  for start in range(0,a.mc,16):
   n=min(16,a.mc-start);shape=(n,3,16,128);scale=np.sqrt(var[None,:,None,:]/2)
   aa=op.apply((rng.normal(size=shape)+1j*rng.normal(size=shape))*scale)
   bb=op.apply((rng.normal(size=shape)+1j*rng.normal(size=shape))*scale)
   zero=grouped(np.real(aa*bb.conj()))
   for key,t in signals.items():samples[key].append(zero+grouped(np.real(t[None]*bb.conj()+aa*t.conj()[None]))+means[key])
  for key in signals:
   q=np.concatenate(samples[key]);C=covs[key];E=np.cov(q,rowvar=False);rat=np.diag(E)/np.diag(C)
   z=(q.mean(0)-means[key])/np.sqrt(np.diag(C)/a.mc)
   reports[name+'_'+key].update({'variance_ratio_min_max':[float(rat.min()),float(rat.max())],'max_absolute_mean_z':float(abs(z).max()),'max_covariance_residual_in_correlation_units':float(np.max(abs(E-C)/np.sqrt(np.outer(np.diag(C),np.diag(C)))))})
   arrays[name+'_'+key+'_MC']=q
  print('Completed',name,flush=True)
 with (a.out/'diagnostics.csv').open('w') as fh:
  w=csv.writer(fh);w.writerow(['operator','foreground_case','group','delay_bin','relative_supported','baseline_signal_Jy2','processed_signal_Jy2','signal_power_ratio','residual_foreground_Jy2','signal_foreground_cross_Jy2','bias_vs_processed_signal_Jy2','conditional_noise_sd_Jy2','bias_over_conditional_noise_sd']);w.writerows(rows)
 np.savez_compressed(a.out/'results.npz',**arrays)
 dump(a.out/'summary.json',{'status':'COMPLETED_CONDITIONAL_FOREGROUND_PILOT_REVIEW_REQUIRED','mc':a.mc,'seed':a.seed,'source_count':len(cat),'cases':reports,'notes':['Fixed signal and foreground; covariance conditional on them, excludes foreground/sample-variance ensemble','Synthetic or supplied finite catalogue; Gaussian beam; no diffuse sky, gain solving, flagging or true station beam','Quadratic complex-visibility fit uses known inverse-noise weights; fixed linear operator, not data-adaptive basis selection','Delay power ratios can exceed one through redistribution; not information retention fractions','No F-Mplus distance, detection significance, 1000h scaling or covariance inversion']})
 dump(a.out/'configuration.json',{'source_input_sha256':sha(inp),'catalogue_sha256':sha(a.out/'catalogue.csv'),'mc':a.mc,'seed':a.seed,'python':platform.python_version(),'numpy':np.__version__,'fit':'Legendre degree 0..2 in linearly normalized frequency, weighted by inverse complex I variance at each time, before existing mean/Hann/FFT','split_seconds':10,'channels':128,'bandwidth_hz':62500,'mean_and_bin_definition':'same 4 groups x 4 baselines x 3 times; 7 absolute-delay bins; equal record/mode weights','catalogue_spectrum':'true S=S140 exp(alpha log(nu/140MHz)+curvature log(nu/140MHz)^2); model S=S140(1+flux_error) exp((alpha+alpha_error) log(nu/140MHz))','catalogue_origin':'synthetic seeded 24-source stress catalogue' if a.catalogue is None else str(a.catalogue),'limits':'All stations share Gaussian response and time-dependent uploaded-table sensitivity; model consistency and foreground field extent not validated'})
 shutil.copy2(__file__,a.out/'run_foregrounds.py');shutil.copy2(inp,a.out/'frozen_inputs.npz');shutil.copy2(ROOT/'inputs/provenance.json',a.out/'input_provenance.json')
 archive=a.out.with_name(a.out.name+'_review.zip')
 with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED) as z:
  for p in sorted(a.out.iterdir()):
   if p.is_file():z.write(p,p.name)
 print('Review archive:',archive.resolve())
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--mc',type=int,default=2048);p.add_argument('--seed',type=int,default=20261021);p.add_argument('--catalogue',type=Path);a=p.parse_args();assert a.mc>=128;run(a)
