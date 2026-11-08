#!/usr/bin/env python3
"""Illustrative spectral-window response and overlap tests, NOT SKA joint covariance."""
import argparse,runpy
from pathlib import Path
import numpy as np
G=runpy.run_path(str(Path(__file__).resolve().parent/'common.py'));globals().update({k:v for k,v in G.items() if not k.startswith('__')})
def operator(n,flags,taper):
    w=np.ones(n) if taper=='rectangular' else np.hanning(n)
    w[np.asarray(flags,dtype=int)]=0
    gain=np.sqrt(np.mean(w*w));check(gain>0,'All channels removed')
    F=np.fft.fft(np.eye(n),axis=0,norm='ortho');A=(F*(w/gain)[None,:])@F.conj().T
    return A,w,gain

def overlap(centres,bandwidth,n):
    step=bandwidth/n;lo=min(centres)-bandwidth/2;hi=max(centres)+bandwidth/2
    grid=lo+(np.arange(round((hi-lo)/step))+.5)*step;weights=[]
    for c in centres:
        first=round((c-bandwidth/2-lo)/step);check(abs((c-bandwidth/2-lo)/step-first)<1e-8,'Bands must share the channel lattice')
        w=np.zeros(len(grid));w[first:first+n]=1;w/=np.linalg.norm(w);weights.append(w)
    W=np.array(weights);return W@W.T

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path('frequency_analysis/response_diagnostics'));p.add_argument('--draws',type=int,default=4096)
    a=p.parse_args();check(a.draws>=1024,'Use at least 1024 draws');check(not a.out.exists() or not any(a.out.iterdir()),'Choose new/empty output');a.out.mkdir(parents=True,exist_ok=True)
    n=128;nu=130+(np.arange(n)+.5-n/2)*8/n;rng=np.random.default_rng(24680);rows=[];responses={}
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for ax,(name,flags,taper) in zip(axs,[('clean_rectangular',[],'rectangular'),('clean_hann',[],'hann'),('flagged_hann',sorted(set(list(range(4))+list(range(124,128))+list(range(15,128,16)))),'hann')]):
        A,w,gain=operator(n,flags,taper);mix=abs(A)**2;cov=A@A.conj().T
        # Inject independent circular complex Gaussian Fourier amplitudes with covariance I.
        X=(rng.normal(size=(n,a.draws))+1j*rng.normal(size=(n,a.draws)))/np.sqrt(2);Y=A@X;emp=Y@Y.conj().T/a.draws
        rel=float(np.linalg.norm(emp-cov)/np.linalg.norm(cov))
        # A deterministic unit-normalized smooth spectral foreground; arbitrary amplitude.
        fg=(nu/130)**-2.6;fg=fg/np.linalg.norm(fg);out=np.fft.fft(w/gain*fg,norm='ortho')
        mode=np.rint(np.fft.fftfreq(n)*n).astype(int);outside=abs(mode)>5
        # Boundaries here are Fourier mode indices, not a horizon wedge or physical k cut.
        rows.append(dict(scenario=name,unflagged_fraction=float(1-len(flags)/n),nonzero_window_fraction=float(np.mean(w!=0)),window_rms_before_normalization=float(gain),covariance_mc_relative_frobenius_error=rel,expected_rms_relative_frobenius_error=float(np.trace(cov).real/(np.sqrt(a.draws)*np.linalg.norm(cov))),mean_same_mode_power_response=float(np.mean(np.diag(mix))),mean_power_mixed_elsewhere=float(np.mean(mix.sum(axis=0)-np.diag(mix))),smooth_fg_power_abs_mode_gt5=float(np.sum(abs(out[outside])**2))))
        responses.update({name+'_A':A,name+'_power_mixing':mix,name+'_white_noise_covariance':cov,name+'_window':w,name+'_foreground_output_power':abs(out)**2})
        ax.imshow(np.log10(np.maximum(np.fft.fftshift(mix),1e-12)),origin='lower',vmin=-8,vmax=0,extent=[-64,64,-64,64],cmap='magma');ax.set(title=name,xlabel='Input Fourier index',ylabel='Output Fourier index')
    fig.savefig(a.out/'spectral_power_mixing.png',dpi=170);plt.close(fig);np.savez_compressed(a.out/'spectral_response.npz',**responses);csv_write(a.out/'spectral_response_summary.csv',rows)
    centres=np.arange(55,196,5);corr=overlap(centres,8,128);ov=[]
    for i,x in enumerate(centres):
        for j,y in enumerate(centres):ov.append(dict(centre_a_MHz=int(x),centre_b_MHz=int(y),overlap_MHz=max(0.,8-abs(float(x-y))),rectangular_white_amplitude_correlation=float(corr[i,j]),proper_complex_white_power_correlation=float(abs(corr[i,j])**2)))
    csv_write(a.out/'band_overlap.csv',ov)
    fig,ax=plt.subplots(figsize=(7,6),layout='constrained');im=ax.imshow(corr,origin='lower',vmin=0,vmax=1);fig.colorbar(im,ax=ax,label='Illustrative white-amplitude correlation');ticks=np.arange(0,len(centres),4);ax.set_xticks(ticks,centres[ticks]);ax.set_yticks(ticks,centres[ticks]);ax.set(xlabel='Band centre [MHz]',ylabel='Band centre [MHz]');fig.savefig(a.out/'band_overlap.png',dpi=160);plt.close(fig)
    json_write(a.out/'interpretation.json',dict(status='ILLUSTRATIVE_ONLY',seed=24680,draws=a.draws,n_channels=n,bandwidth_MHz=8,frequency_MHz=130,normalization='Windows divided by RMS; matrices are uncorrected transfer response, not deconvolved bandpowers.',notes=['A = F diag(window/RMS) F^H; |A|^2 maps independent input Fourier powers.','A A^H is covariance only for unit circular complex white amplitudes.','Power-correlation = |amplitude-correlation|^2 applies only to the toy proper-complex white Gaussian case.','Smooth foreground has unit spectral L2 norm; leakage has no calibrated sky amplitude.','No baseline chromaticity, primary-beam sky coupling, wedge calibration, ionospheric response, or sample-variance evolution is modeled.','Do not insert these overlap coefficients into the cosmological forecast covariance.','Forecast_bands continues to use the untapered reference convention; this diagnostic does not alter its errors or masks.','A measured pipeline transfer matrix and full cross-band covariance are still needed before a joint significance claim.']))
if __name__=='__main__':main()
