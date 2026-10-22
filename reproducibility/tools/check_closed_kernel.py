"""Independent algebraic simplification, positivity and endpoint-law checks."""
import sys, json, time, hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.physics import bound_energy,dissociation_kernel,born_rate,cutoff_boundary

def closed_kernel(q,K,U,J=1.):
 E=bound_energy(K,U,J)
 t=2*J*np.cos((K+q)/2)
 D2=4*t*t-E*E
 if D2<=0:return 0.
 D=np.sqrt(D2)
 c=np.cos(K/2)
 return float(8*(U/E)**3*D/(D2+U*U)/(1-4*c*c*t*t/(E*E))**2)

def edge_amplitude(K,U,J=1.):
 E=bound_energy(K,U,J);Q=2*np.arccos(E/(4*J))
 return 8*U/(E**3*np.sin(K/2)**4)*np.sqrt(8*J*J*np.sin(Q))

def main():
 start=time.perf_counter();cpu=time.process_time()
 rng=np.random.default_rng(513784);rows=[]
 for _ in range(1500):
  K=rng.uniform(.01,np.pi);U=rng.uniform(.05,3.99);q=rng.uniform(-np.pi,np.pi)
  a=dissociation_kernel(q,K,U);b=closed_kernel(q,K,U)
  if max(a,b)>1e-12:rows.append(abs(a-b)/max(a,b))
 K=.75*np.pi;threshold=4*np.sqrt(-np.cos(K));edges=[]
 for U,mode in [(2.,'wrap_open'),(threshold,'critical'),(3.5,'internal')]:
  bc=cutoff_boundary(K,U);qc=bc['qmax'];A=edge_amplitude(K,U)
  vals=[]
  for d in [1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5]:
   q0=qc-d
   c,_=born_rate(K,U,q0,epsabs=1e-14,epsrel=1e-9)
   cp=c*(1-q0/np.pi)
   expected=A/(3*np.pi)*(d**1.5)/(1-qc/np.pi) if mode=='internal' else (A/3*np.sqrt(d) if mode=='critical' else closed_kernel(np.pi,K,U))
   vals.append(dict(delta=d,c=c,passband=cp,asymptotic=expected,ratio=c/expected))
  slope=np.polyfit(np.log([r['delta'] for r in vals[-3:]]),np.log([r['c'] for r in vals[-3:]]),1)[0]
  edges.append(dict(U=U,mode=mode,qc_over_pi=qc/np.pi,A=A,slope_last3=float(slope),rows=vals))
 out=dict(max_relative_kernel_error=max(rows),positive_points=len(rows),edges=edges,wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu)
 target=Path(__file__).resolve().parents[1]/'results'/'theory_audit'
 target.mkdir(parents=True,exist_ok=True)
 dest=target/'closed_kernel_initial.json'
 if dest.exists():raise FileExistsError(dest)
 dest.write_text(json.dumps(out,indent=2),encoding='utf-8')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
