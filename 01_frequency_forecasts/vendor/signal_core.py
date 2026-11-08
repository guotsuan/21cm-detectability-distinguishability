#!/usr/bin/env python3
"""Seven-model real-space spectra from existing coeval cubes; no simulation rerun.
Requires NumPy and Matplotlib only. Frequencies are checked from stored redshifts.
P = L^3/N^6 |FFT(T-mean(T))|^2; Delta2 = k_mean^3 * P_shell / (2*pi^2).
Full complex FFT counts BOTH conjugate cells. Counts are simulation cells,
NOT observational mode counts or covariance. No observational errors are inferred.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter

MODELS = ['F','Hminus','Hplus','Mminus','Mplus','Iminus','Iplus']
NU21 = 1420.40575177


def fft_power(cube, length):
    field = np.asarray(cube, dtype=np.float64)
    if field.ndim != 3 or len(set(field.shape)) != 1 or not np.isfinite(field).all():
        raise ValueError('brightness_temp must be a finite cubic 3D field')
    if not np.isfinite(length) or length <= 0:
        raise ValueError('Box length must be positive')
    n = field.shape[0]
    mean = field.mean()
    field = field - mean
    ft = np.fft.fftn(field)
    ft[0,0,0] = 0
    power = length**3 / field.size**2 * np.abs(ft)**2
    var_real = np.mean(field**2)
    var_fft = power.sum()/length**3
    # Absolute tolerance only matters for a field with essentially zero variance.
    if not np.isclose(var_real, var_fft, rtol=1e-10, atol=1e-20):
        raise AssertionError(f'Parseval failed: real={var_real}, FFT={var_fft}')
    axis = 2*np.pi*np.fft.fftfreq(n, d=length/n)
    k = np.sqrt(axis[:,None,None]**2+axis[None,:,None]**2+axis[None,None,:]**2)
    return k.ravel(), power.ravel(), float(mean), float(var_real), float(var_fft)


def spherical(k, power, edges):
    # Last bin includes its upper edge; DC is always excluded.
    sel = (k > 0) & (k >= edges[0]) & (k <= edges[-1])
    ki, pi = k[sel], power[sel]
    idx = np.searchsorted(edges, ki, side='right')-1
    idx[idx == len(edges)-1] -= 1
    nb = len(edges)-1
    count = np.bincount(idx, minlength=nb)
    sums = np.bincount(idx, weights=pi, minlength=nb)
    ksum = np.bincount(idx, weights=ki, minlength=nb)
    km = np.divide(ksum, count, out=np.full(nb,np.nan), where=count>0)
    pm = np.divide(sums, count, out=np.full(nb,np.nan), where=count>0)
    delta = km**3*pm/(2*np.pi**2)
    return km, pm, delta, count, sums


def load_model(root, model, frequency):
    folder = root/model
    config = json.loads((folder/'requested_config.json').read_text())
    if config['model'] != model or config['coordinate_space'] != 'real':
        raise ValueError(f'{model}: model label or real-space configuration mismatch')
    target_z = NU21/frequency - 1
    candidates = []
    for f in folder.glob('coeval_z*.npz'):
        try:
            if abs(float(f.stem.split('_z')[1])-target_z) <= .00051:
                candidates.append(f)
        except ValueError:
            continue
    if len(candidates) != 1:
        raise FileNotFoundError(f'{model}: expected one cube near z={target_z:.8f}, found {candidates}')
    path = candidates[0]
    with np.load(path, allow_pickle=False) as d:
        z = float(d['redshift'].item())
        nu = NU21/(1+z)
        if abs(nu-frequency) > 1e-4:
            raise ValueError(f'{model}: stored redshift gives {nu:.8f} MHz, not {frequency}')
        length = float(d['box_len_cMpc'].item())
        seed = int(d['random_seed'].item())
        if str(d['coordinate_space'].item()) != 'real':
            raise ValueError(f'{model}: cube is not real-space')
        cube = np.array(d['brightness_temp'],dtype=np.float64)
    n = config['user_params']['HII_DIM']
    if cube.shape != (n,n,n) or not np.isclose(length,config['user_params']['BOX_LEN']):
        raise ValueError(f'{model}: cube geometry disagrees with configuration')
    if seed != config['random_seed']:
        raise ValueError(f'{model}: seed disagrees with configuration')
    return cube,length,config,dict(model=model,cube_path=str(path.resolve()),redshift=z,
        frequency_MHz=nu,random_seed=seed,box_len_cMpc=length,N=n,
        cube_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def self_test():
    n=16; length=80.; amp=3.
    cube=np.broadcast_to(amp*np.cos(2*np.pi*2*np.arange(n)/n)[:,None,None],(n,n,n))+7.
    k,p,mean,vr,vf=fft_power(cube,length)
    assert np.isclose(mean,7) and np.isclose(vr,amp**2/2)
    assert np.isclose(p.max(),length**3*amp**2/4)
    edges=np.array([.01,.14,.18,.7])
    km,pm,de,cnt,sums=spherical(k,p,edges)
    assert np.isclose(sums.sum()/length**3,amp**2/2)
    assert np.isclose(sums[1]/length**3,amp**2/2)
    _,p0,_,v0,_=fft_power(np.ones((8,8,8))*2,40)
    assert v0==0 and np.all(p0==0)
    print('PASS: analytic sinusoid normalization, binning, Parseval, and constant-field checks.')

