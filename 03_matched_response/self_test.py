#!/usr/bin/env python3
"""No simulation: estimator identity, covariance scaling, and geometry conventions."""
import numpy as np
from common import fringe,gaussian,load_reference,C
from cleaner import Cleaner
import estimator as est

def main():
 rng=np.random.default_rng(91);v=rng.normal(size=(3,16,128))+1j*rng.normal(size=(3,16,128));nu=126e6+(np.arange(128)+.5)*62500
 var=np.broadcast_to(np.linspace(2,4,128),(3,128)).copy();ops=[Cleaner(var,nu,k) for k in [-1,1,2]]
 for op in ops:
  np.testing.assert_allclose(op.apply(v),np.einsum('tij,tbj->tbi',op.L,v),rtol=1e-12,atol=2e-12)
  N,M=est.theory(var,op.L,op.apply(v));N2,M2=est.theory(var*2,op.L,op.apply(v))
  np.testing.assert_allclose(N2,4*N,atol=1e-12);np.testing.assert_allclose(M2,2*M,atol=1e-12)
  assert np.linalg.eigvalsh(N).min()>-1e-12
  # Cross-power estimator averages 12 records/group and modes/bin.
  expected=[]
  f=abs(op.apply(v))**2
  for g in range(4):
   for mask in est.MASKS:expected.append(f[:,g*4:(g+1)*4][:,:,mask].mean())
  np.testing.assert_allclose(op.power(v),expected,rtol=1e-14)
  # Changing reference inside log does not change scaled Legendre coordinates.
  x=np.log(nu/140e6);y=np.log(nu/130e6)
  np.testing.assert_allclose(2*(x-x.min())/np.ptp(x)-1,2*(y-y.min())/np.ptp(y)-1,atol=1e-13)
 a=np.array([[0.,0.,1.],[.01,0,np.sqrt(1-.01**2)]])
 pr,d,_,_=load_reference();k1=fringe(d['uvw_m'][0],130e6,a,pr['phase_centre_dec_deg']);k2=fringe(d['uvw_m'][0],140e6,a,pr['phase_centre_dec_deg'])
 np.testing.assert_allclose(k1[:,0],1);assert not np.allclose(k1[:,1],k2[:,1]);assert gaussian(130e6,a)[1]>gaussian(140e6,a)[1]
 assert sum(int(m.sum()) for m in est.MASKS)==128
 assert len(est.MASKS)==7
 print('PASS: operator identity, equal-record/bin averaging, noise/mixed covariance scaling, PSD, frequency-dependent fringes and Gaussian scaling; all 128 delay modes assigned')
if __name__=='__main__':main()
