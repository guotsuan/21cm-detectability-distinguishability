"""Project seven models and evaluate the frozen estimator with refreshed weights."""
import itertools,json,csv
from pathlib import Path
import numpy as np
from common import load,sha,dump,fringe,gaussian,C,ROOT
from cleaner import Cleaner
import estimator as est
MODELS=['F','Hminus','Hplus','Mminus','Mplus','Iminus','Iplus']
BEAMS=['dipole','isotropic','gaussian'];PIPE=['fixed_catalogue','fit_linear','fit_quadratic']

def project(folder,seed):
 pr,d,cat,_=load();prep=folder/f'seed{seed}/prepared';meta=json.loads((prep/'metadata.json').read_text())
 if sha(prep/'physical.npz')!=meta['physical_sha256']:raise ValueError('Prepared data changed')
 with np.load(prep/'physical.npz') as h:
  fields=h['temperatures_mK'].astype(float);dirs=h['directions'];solid=h['solid_angle_sr'];nu=h['frequencies_hz']
 np.testing.assert_array_equal(nu,d['frequencies_hz'])
 with np.load(folder/'geometry.npz') as g:np.testing.assert_array_equal(dirs,g['directions'])
 with np.load(folder/'beams/beam_power.npz') as h:beam=h['power']
 done=json.loads((folder/'beams/complete.json').read_text())
 if done['power_sha256']!=sha(folder/'beams/beam_power.npz') or done['geometry_sha256']!=sha(folder/'geometry.npz'):raise ValueError('Beam/geometry changed')
 v=np.empty((3,14,3,16,128),complex);fg=np.empty((3,3,16,128),complex)
 for ch,f in enumerate(nu):
  x=np.log(f/140e6);sf=cat[:,2]*np.exp(cat[:,3]*x+cat[:,4]*x*x)
  flux=fields[:,ch]*solid[None]*2*1.380649e-23*f*f/C**2*1e23
  for t in range(3):
   K=fringe(d['uvw_m'][t],f,dirs[1:],pr['phase_centre_dec_deg'])
   for mi,name in enumerate(BEAMS):
    response=gaussian(f,dirs) if name=='gaussian' else beam[mi,t,ch]
    fg[mi,t,:,ch]=K[:,:24]@(sf*response[1:25])
    v[mi,:,t,:,ch]=(K[:,24:]@(flux*response[25:]).T).T
  if (ch+1)%16==0:print(folder.name,'seed',seed,'projected',ch+1,'/128',flush=True)
 return v,fg,meta

def noise_check(var,ops,mc,seed):
 rng=np.random.default_rng(seed);samples=[[] for _ in ops]
 for start in range(0,mc,16):
  shape=(min(16,mc-start),3,16,128);scale=np.sqrt(var[None,:,None,:]/2)
  a=(rng.normal(size=shape)+1j*rng.normal(size=shape))*scale
  b=(rng.normal(size=shape)+1j*rng.normal(size=shape))*scale
  for j,op in enumerate(ops):samples[j].append(est.grouped(np.real(op.apply(a)*op.apply(b).conj())))
 arr=np.array([np.concatenate(s) for s in samples]);theory=[];rows=[]
 for j,op in enumerate(ops):
  N,_=est.theory(var,op.L,np.zeros((3,16,128),complex));E=np.cov(arr[j],rowvar=False);rat=np.diag(E)/np.diag(N)
  err=float(np.max(abs(E-N)/np.sqrt(np.outer(np.diag(N),np.diag(N)))))
  theory.append(N);rows.append(dict(pipeline=PIPE[j],variance_ratio_range=[float(rat.min()),float(rat.max())],max_covariance_residual_in_correlation_units=err,max_mean_z=float(np.max(abs(arr[j].mean(0))/np.sqrt(np.diag(N)/mc)))))
 return np.array(theory),arr,rows

def run(folder,mc):
 pr,d,_,_=load();var=d['complex_variance_I_Jy2'];nu=d['frequencies_hz'];centre=int(round(nu.mean()/1e6));ops=[Cleaner(var,nu,k) for k in [-1,1,2]]
 for op in ops:
  x=np.random.default_rng(11).normal(size=(3,16,128))+1j*np.random.default_rng(12).normal(size=(3,16,128))
  np.testing.assert_allclose(op.apply(x),np.einsum('tij,tbj->tbi',op.L,x),atol=2e-12,rtol=2e-12)
 noiseout=folder/'noise_validation'
 if not (noiseout/'complete.json').exists():
  noiseout.mkdir(exist_ok=True);N,draws,checks=noise_check(var,ops,mc,20261020+centre)
  np.savez_compressed(noiseout/'noise_validation.npz',analytic_covariance_Jy4=N,noise_crosspower_Jy2=draws,operators=np.array([o.L for o in ops]))
  dump(noiseout/'complete.json',dict(status='COMPLETE_REVIEW_REQUIRED',mc=mc,noise_sha256=sha(folder/'noise/noise.npz'),checks=checks))
 else:
  check=json.loads((noiseout/'complete.json').read_text())
  if check['mc']!=mc or check['noise_sha256']!=sha(folder/'noise/noise.npz'):raise ValueError('Noise MC resume differs')
  N=np.load(noiseout/'noise_validation.npz')['analytic_covariance_Jy4']
 for seed in [12345,23456,34567]:
  out=folder/f'seed{seed}/results';out.mkdir(exist_ok=True)
  dependencies={n:sha(p) for n,p in [('physical',folder/f'seed{seed}/prepared/physical.npz'),('noise',folder/'noise/noise.npz'),('beam',folder/'beams/beam_power.npz')]}
  if (out/'complete.json').exists():
   done=json.loads((out/'complete.json').read_text())
   if done['dependencies']!=dependencies or done['results_sha256']!=sha(out/'results.npz'):raise ValueError('Completed output changed')
   print('Reuse checked result',centre,seed);continue
  v,fg,meta=project(folder,seed)
  power=np.array([[op.power(v[mi]) for op in ops] for mi in range(3)]).transpose(0,2,1,3)
  if power.shape!=(3,14,3,28):raise ValueError('Power axes mismatch')
  regression={}
  if centre==140:
   refpath=ROOT/f'inputs/reference140_seed{seed}.npz'
   with np.load(refpath) as old:
    for key,now in [('signal_power_Jy2',power),('foreground_Jy',fg)]:
     error=float(np.linalg.norm(now-old[key])/max(np.linalg.norm(old[key]),1e-300));regression[key]=error
     if error>1e-8:raise ValueError('140 MHz reproduction failed: '+key+' '+str(error))
  bias=np.empty((14,3,2,28));mixed=np.empty((14,3,2,28,28));residual=np.empty((3,28));cross=np.empty((14,3,28));injection=[]
  for j,op in enumerate(ops):
   h=op.apply(fg[0]-fg[2]);residual[j]=est.grouped(abs(h)**2)
   for f in range(14):
    s=op.apply(v[0,f]);cross[f,j]=est.grouped(2*np.real(s*h.conj()))
    for c,r in enumerate([np.zeros_like(h),h]):
     bias[f,j,c]=est.grouped(abs(r)**2+2*np.real(s*r.conj()))
     nn,mm=est.theory(var,op.L,s+r);mixed[f,j,c]=nn+mm
    # Paired +/- injections measure the estimator signal response for a fixed linear fit.
    alpha=100.;rec=(est.grouped(abs(h+alpha*s)**2)+est.grouped(abs(h-alpha*s)**2)-2*residual[j])/(2*alpha**2)
    supported=power[0,f,j]>max(power[0,f,j].max()*1e-10,1e-300)
    injection.append(dict(field=meta['field_order'][f],pipeline=PIPE[j],amplitude=alpha,max_relative_error=float(np.max(abs(rec[supported]-power[0,f,j,supported])/power[0,f,j,supported])) if supported.any() else None))
  pairs=list(itertools.combinations(range(7),2));delta=np.array([power[:,2*i:2*i+2]-power[:,2*k:2*k+2] for i,k in pairs]).transpose(1,0,2,3,4)
  np.savez_compressed(out/'results.npz',visibilities_Jy=v,foreground_Jy=fg,signal_power_Jy2=power,pairwise_power_difference_Jy2=delta,conditional_covariance_Jy4=mixed,bias_Jy2=bias,foreground_residual_Jy2=residual,signal_foreground_cross_Jy2=cross)
  dump(out/'configuration.json',dict(frequency_MHz=centre,seed=seed,physical_metadata=meta,beam_order=BEAMS,field_order=meta['field_order'],pipeline_order=PIPE,power_axes=['beam','field','pipeline','group_times_delay28'],conditional_covariance_axes=['field','pipeline','matched_or_gaussian_residual','output','output'],estimator='slice transverse mean removal; fixed/weighted log-frequency degree1/2; unweighted frequency mean; RMS-normalised Hann; unitary FFT; equal average over 3 times x 4 baselines per group and Fourier modes per bin',angular_support='fixed 150 cMpc square; no common-angular-window claim',noise='independent equal-variance station model per time/channel; same 10-second split; thermal and fixed-sky mixed covariance only',limits='No sample variance ensemble covariance, diffuse foreground, calibration solver, station-specific SKALA or full survey significance',delay_edges_us=(est.EDGES*1e6).tolist(),regression140=regression,symmetric_injection=injection))
  dump(out/'complete.json',dict(status='COMPLETE_REVIEW_REQUIRED',dependencies=dependencies,results_sha256=sha(out/'results.npz')))
 print('Analysis complete',folder,flush=True)
