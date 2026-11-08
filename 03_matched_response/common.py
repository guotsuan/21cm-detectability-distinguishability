import csv,hashlib,json,zipfile,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
C=299792458.
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_reference():
 pr=json.loads((ROOT/'inputs/provenance.json').read_text())
 for n,h in pr['input_sha256'].items():
  if sha(ROOT/'inputs'/n)!=h:raise ValueError('Input hash mismatch: '+n)
 d=np.load(ROOT/'inputs/frozen.npz');cat=np.loadtxt(ROOT/'inputs/catalogue.csv',delimiter=',',skiprows=1)
 fg=np.c_[cat[:,:2],np.sqrt(1-np.sum(cat[:,:2]**2,axis=1))]
 dirs=np.vstack([[0,0,1],fg,d['signal_directions_lmn']])
 return pr,d,cat,dirs

def load():
 pr,old,cat,dirs=load_reference()
 context=Path(os.environ['EOR_MATCHED_CONTEXT'])
 with np.load(context/'noise.npz') as z:
  d={k:old[k] for k in old.files}
  d['frequencies_hz']=z['frequencies_hz'].copy()
  d['complex_variance_I_Jy2']=z['complex_variance_I_Jy2'].copy()
 pr['matched_noise_sha256']=sha(context/'noise.npz')
 return pr,d,cat,dirs

def radec(dirs,ra0,dec0):
 a,d=np.deg2rad([ra0,dec0]);u=np.array([-np.sin(a),np.cos(a),0.]);v=np.array([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)]);w=np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
 basis=np.stack([u,v,w]);xyz=dirs@basis
 ra=np.arctan2(xyz[:,1],xyz[:,0]);dec=np.arcsin(np.clip(xyz[:,2],-1,1))
 back=np.c_[np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)]@basis.T
 np.testing.assert_allclose(back,dirs,atol=2e-14)
 return np.c_[np.rad2deg(ra)%360,np.rad2deg(dec)]
def fringe(uvw,nu,dirs,dec_deg):
 l,m,n=dirs.T;distance=uvw@(dirs-[0,0,1]).T;dec=np.deg2rad(dec_deg)
 du=-np.sin(dec)*uvw[:,1]+np.cos(dec)*uvw[:,2]
 rate=7.272205217e-5*(du[:,None]*l+np.sin(dec)*uvw[:,0,None]*m-np.cos(dec)*uvw[:,0,None]*(n-1))
 return np.exp(2j*np.pi*nu*distance/C)*np.sinc(62500*distance/C)*np.sinc(10*nu*rate/C)
def gaussian(nu,dirs):return np.exp(-8*np.log(2)*np.sum(dirs[:,:2]**2,axis=1)/np.sin(np.deg2rad(4)*140e6/nu)**2)
def archive(out):
 dest=out.with_name(out.name+'_review.zip')
 with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
  for f in sorted(out.rglob('*')):
   if f.is_file():z.write(f,f.relative_to(out))
 print('Review archive:',dest,flush=True)
