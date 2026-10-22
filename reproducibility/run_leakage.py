"""One permitted robustness extension: finite leakage in the spectral gap."""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,hashlib,time,sys
import numpy as np
from scipy.sparse.linalg import expm_multiply
from run_confirmation import basis
from src.physics import initial_bound_state,two_body_hamiltonian,bound_population,born_rate

ROOT=Path(__file__).resolve().parent
def job(seed,epsilon):
 start=time.perf_counter();cpu=time.process_time();n=256;U=3.5;W=.02;f=.97
 qs=2*np.pi*np.fft.fftfreq(n);shape=np.where(np.abs(qs)>f*np.pi+1e-12,1.,epsilon);shape[0]=0.
 raw=np.random.default_rng(seed).standard_normal(n);v=W*np.fft.ifft(np.fft.fft(raw)*np.sqrt(shape/shape.mean())).real
 b=basis(n,U);psi=initial_bound_state(b,3*n//8);h=two_body_hamiltonian(v,U);ts=np.arange(0.,33.,2.)
 states=expm_multiply(-1j*h,psi,start=0.,stop=32.,num=len(ts),traceA=-1j*h.diagonal().sum())
 norm=float(np.max(abs(np.linalg.norm(states,axis=1)-1)));sym=float(np.linalg.norm(states[-1].reshape(n,n)-states[-1].reshape(n,n).T))
 modespath=sorted((ROOT/'results/sample_diagnostic').glob('*/modes_N256_U3.5.npz'))[-1];modes=np.load(modespath)['mode_kernels']
 ell=(abs(np.fft.fft(v/W))**2/n**2)@modes
 loss=np.array([1-bound_population(s,b) for s in states])
 return dict(n=n,U=U,W=W,f=f,seed=seed,epsilon=epsilon,times=ts.tolist(),loss=loss.tolist(),predicted_loss=(W*W*ell).tolist(),
   spectral_normalizer=float(shape.mean()),sample_variance=float(v@v/n),max_norm_error=norm,exchange_error=sym,
   potential_sha256=hashlib.sha256(v.tobytes()).hexdigest(),relative_curve_error=float(np.linalg.norm(loss-W*W*ell)/np.linalg.norm(W*W*ell)),
   wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu,valid=bool(norm<1e-9 and sym<1e-9))

def main():
 run=ROOT/'results/leakage'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');run.mkdir(parents=True)
 start=time.perf_counter();cpu=time.process_time();rows=[]
 srcs=[Path(__file__),ROOT/'run_confirmation.py',ROOT/'src/physics.py',ROOT/'notes/G3_LEAKAGE_PROTOCOL.md']
 meta=dict(status='RUNNING',source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in srcs},command=sys.argv)
 (run/'manifest.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
 with ProcessPoolExecutor(max_workers=4) as pool:
  fs=[pool.submit(job,s,e) for e in [0.,.003,.01] for s in range(200000,200016)]
  for f in as_completed(fs):
   row=f.result();assert row['valid'];rows.append(row)
   (run/f"epsilon{row['epsilon']}_seed{row['seed']}.json").write_text(json.dumps(row,indent=2),encoding='utf-8')
   if len(rows)%16==0:print('completed',len(rows),flush=True)
 cwhite=born_rate(.75*np.pi,3.5,0.)[0]
 analytic=[dict(epsilon=e,c=e*cwhite/(.03+.97*e),relative_to_white=e/(.03+.97*e)) for e in [0.,.001,.003,.01,.03]]
 meta.update(status='PASS',trajectories=len(rows),wall_seconds=time.perf_counter()-start,parent_cpu_seconds=time.process_time()-cpu,
   worker_cpu_seconds=sum(r['process_cpu_seconds'] for r in rows),analytic=analytic,
   diagnostics=[dict(epsilon=e,mean_curve_error=float(np.mean([r['relative_curve_error'] for r in rows if r['epsilon']==e]))) for e in [0.,.003,.01]],
   output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in run.glob('epsilon*.json')})
 (run/'manifest.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
 print('Saved',run,meta['diagnostics'],flush=True)
if __name__=='__main__':main()
