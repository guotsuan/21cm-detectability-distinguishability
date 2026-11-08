"""Gaussian likelihood utilities; algebra only, not a calibrated EoR likelihood."""
import numpy as np

def loglike(data,mean,cov):
    data=np.asarray(data,float);mean=np.asarray(mean,float);cov=np.asarray(cov,float)
    if cov.shape!=(data.size,data.size) or mean.shape!=data.shape:raise ValueError('Incompatible dimensions')
    if not np.isfinite(cov).all() or not np.allclose(cov,cov.T,rtol=1e-10,atol=1e-12):raise ValueError('Invalid covariance')
    L=np.linalg.cholesky(cov);r=np.linalg.solve(L,data-mean)
    return -.5*(r@r+2*np.log(np.diag(L)).sum()+data.size*np.log(2*np.pi))

def expected_2logLR_diagonal(truth,candidate,var_truth,var_candidate):
    """2 KL[N(truth,Ctruth)||N(candidate,Ccandidate)], not sigma squared."""
    f,m,vf,vm=map(lambda x:np.asarray(x,float),(truth,candidate,var_truth,var_candidate))
    if not np.isfinite([f,m,vf,vm]).all() or not np.all(vf>0) or not np.all(vm>0):raise ValueError('Invalid mean/variance')
    ratio=vf/vm
    return float(np.sum(ratio-1-np.log(ratio)+(f-m)**2/vm))
