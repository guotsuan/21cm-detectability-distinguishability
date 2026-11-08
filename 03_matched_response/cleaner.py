import numpy as np
import estimator as est
class Cleaner:
 def __init__(self,var,nu,degree):
  self.degree=degree;self.var=var;self.Q=np.tile(np.eye(128),(3,1,1))
  if degree>=0:
   x=np.log(nu/140e6);x=2*(x-x.min())/(x.max()-x.min())-1
   X=np.polynomial.legendre.legvander(x,degree)
   for t in range(3):
    A=X/np.sqrt(var[t])[:,None];U,s,Vh=np.linalg.svd(A,full_matrices=False)
    if s[-1]<s[0]*1e-10:raise ValueError('Ill-conditioned fit basis')
    pinv=(Vh.T/s)@U.T
    self.Q[t]-=X@pinv/np.sqrt(var[t])[None,:]
    np.testing.assert_allclose(self.Q[t]@X,0,atol=1e-12)
    np.testing.assert_allclose(self.Q[t]@self.Q[t],self.Q[t],atol=1e-12)
  L0=est.Operator(var,False).L;self.L=L0@self.Q
 def residual(self,v):return np.einsum('tij,...tbj->...tbi',self.Q,v,optimize=True)
 def apply(self,v):return est.base(self.residual(v))
 def power(self,v):return est.grouped(abs(self.apply(v))**2)

