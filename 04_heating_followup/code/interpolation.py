from common import *
def interp(t,k):
    t=t.sort_values('k_Mpc_inv');x=t.k_Mpc_inv.to_numpy();v=t.Delta2_mK2.to_numpy()
    require(len(x)>1 and np.all(v>0) and np.all(np.diff(x)>0),'Invalid/duplicate spectral support')
    y=np.exp(np.interp(np.log(k),np.log(x),np.log(v)));y[(k<x.min())|(k>x.max())]=np.nan;return y
