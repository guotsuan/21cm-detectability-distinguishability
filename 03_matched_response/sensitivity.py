import numpy as np
CASES={'baseline128':('linear',128),'cubic128':('cubic4',128),'crop96':('linear',96),'crop64':('linear',64)}
def cubic4(cube,dist,target):
 dist=np.asarray(dist);target=np.asarray(target)
 if len(dist)<4 or not np.all(np.diff(dist)>0):raise ValueError('Need four increasing radial nodes')
 if target.min()<dist[0] or target.max()>dist[-1]:raise ValueError('No extrapolation allowed')
 start=np.clip(np.searchsorted(dist,target,side='right')-2,0,len(dist)-4)
 idx=start[:,None]+np.arange(4);nodes=dist[idx];weights=np.ones_like(nodes)
 for i in range(4):
  for j in range(4):
   if i!=j:weights[:,i]*=(target-nodes[:,j])/(nodes[:,i]-nodes[:,j])
 result=sum(cube[...,idx[:,i]]*weights[:,i] for i in range(4))
 return np.moveaxis(result,-1,0)
def window_field(f,size):
 if f.shape[1:]!=(128,128):raise ValueError('Expected fixed128x128 directions')
 if size not in (128,96,64):raise ValueError('Unsupported patch')
 if size==128:
  mean=f.mean((1,2));return f-mean[:,None,None],mean
 start=(128-size)//2;sl=slice(start,start+size)
 mean=f[:,sl,sl].mean((1,2));out=np.zeros_like(f)
 out[:,sl,sl]=f[:,sl,sl]-mean[:,None,None]
 return out,mean
