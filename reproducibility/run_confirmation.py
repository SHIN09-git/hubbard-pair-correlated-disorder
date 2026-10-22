"""Frozen paired-seed full N^2 dynamics, auditable per-worker CPU costs."""
from pathlib import Path
from datetime import datetime,timezone
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,sys,time,platform,os
import numpy as np
import scipy
from scipy.sparse.linalg import expm_multiply
from src.physics import clean_bound_basis,initial_bound_state,two_body_hamiltonian,bound_population,filtered_noise

ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(maxsize=16)
def basis(n,U):return clean_bound_basis(n,U)

def evolve(case):
    start=time.perf_counter();cpu=time.process_time()
    n=case['n'];U=case['U'];W=case['W'];f=case['f'];seed=case['seed'];norm=case['normalization']
    b=basis(n,U);psi=initial_bound_state(b,3*n//8)
    v=filtered_noise(n,f*np.pi,W,seed)
    modes=int(np.sum(np.abs(2*np.pi*np.fft.fftfreq(n))>f*np.pi+1e-12))
    if norm=='passband':v*=np.sqrt(modes/n)
    ts=np.arange(0.,case.get('tmax',32.)+1.,2.)
    h=two_body_hamiltonian(v,U)
    states=expm_multiply(-1j*h,psi,start=0.,stop=float(ts[-1]),num=len(ts),traceA=-1j*h.diagonal().sum())
    pops=[bound_population(p,b) for p in states]
    norms=np.linalg.norm(states,axis=1)
    sym=float(np.linalg.norm(states[-1].reshape(n,n)-states[-1].reshape(n,n).T))
    row=dict(case,times=ts.tolist(),bound_populations=pops,allowed_modes=modes,
             expected_variance=W*W*(modes/n if norm=='passband' else 1.),
             sample_variance=float(v@v/n),potential_sha256=hashlib.sha256(v.tobytes()).hexdigest(),
             max_norm_error=float(np.max(np.abs(norms-1))),exchange_symmetry_error=sym,
             wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu)
    row['valid']=bool(row['max_norm_error']<1e-9 and sym<1e-9 and np.isfinite(pops).all())
    return row

def main():
    cfg_path=Path(sys.argv[1]);cfg=json.loads(cfg_path.read_text(encoding='utf-8'))
    run=ROOT/'results'/'confirmation'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run.mkdir(parents=True);(run/'config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    sources=list(ROOT.glob('*.py'))+list((ROOT/'src').glob('*.py'))
    protocol=ROOT/'CONFIRMATION_PROTOCOL.md'
    manifest=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),config_sha256=sha(cfg_path),
      protocol_sha256=sha(protocol),source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},
      command=sys.argv,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
      thread_environment={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']})
    (run/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    cases=[dict(c,seed=s) for c in cfg['cases'] for s in cfg['seeds']]
    start=time.perf_counter();cpu_parent=time.process_time();rows=[]
    with ProcessPoolExecutor(max_workers=cfg['workers']) as pool:
        future_to_case={pool.submit(evolve,c):c for c in cases}
        for future in as_completed(future_to_case):
            c=future_to_case[future];row=future.result();rows.append(row)
            key=f"N{c['n']}_U{c['U']}_W{c['W']}_f{c['f']}_{c['normalization']}_s{c['seed']}"
            (run/(key+'.json')).write_text(json.dumps(row,indent=2),encoding='utf-8')
            if len(rows)%16==0:print(len(rows),'/',len(cases),'worker CPU',round(sum(r['process_cpu_seconds'] for r in rows),1),flush=True)
            if not row['valid']:raise RuntimeError('Numerical invariant failed: '+key)
    manifest.update(status='PASS',finished_utc=datetime.now(timezone.utc).isoformat(),trajectories=len(rows),
       wall_seconds=time.perf_counter()-start,parent_cpu_seconds=time.process_time()-cpu_parent,
       worker_cpu_seconds=sum(r['process_cpu_seconds'] for r in rows),
       max_norm_error=max(r['max_norm_error'] for r in rows),max_symmetry_error=max(r['exchange_symmetry_error'] for r in rows),
       output_sha256={p.name:sha(p) for p in run.glob('*.json') if p.name!='manifest.json'})
    (run/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('PASS',run,'wall',manifest['wall_seconds'],'worker CPU',manifest['worker_cpu_seconds'],flush=True)
if __name__=='__main__':main()
