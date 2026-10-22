"""Post-confirmation diagnostic: separate sample powers from finite-W corrections.

This does not change the frozen primary endpoints or exclusions.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import json,hashlib,time
import numpy as np
from src.finite_time import finite_time_loss
from src.physics import filtered_noise
from analyze_confirmation import interval,slope
ROOT=Path(__file__).resolve().parent

def main():
 start=time.perf_counter();cpu=time.process_time()
 run=ROOT/'results/sample_diagnostic'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');run.mkdir(parents=True)
 rs=[];group=defaultdict(list)
 for p in (ROOT/'results/confirmation').glob('*/N*.json'):rs.append(json.loads(p.read_text()))
 ts=np.array(rs[0]['times']);cache={};rows=[]
 for n,U in sorted(set((r['n'],r['U']) for r in rs)):
  mode=np.array(finite_time_loss(n,.75*np.pi,U,ts,return_modes=True)['mode_kernels'])
  np.savez_compressed(run/f'modes_N{n}_U{U}.npz',times=ts,mode_kernels=mode)
  cache[n,U]=mode
 for r in rs:
  n,U,W,f,norm=(r[x] for x in ['n','U','W','f','normalization']);v=filtered_noise(n,f*np.pi,W,r['seed'])
  if norm=='passband':v*=np.sqrt(r['allowed_modes']/n)
  assert hashlib.sha256(v.tobytes()).hexdigest()==r['potential_sha256']
  ell=(abs(np.fft.fft(v/W))**2/n**2)@cache[n,U]
  actual=(1-np.array(r['bound_populations']))/W**2
  w=(12,24) if U==2 else (16,28);mask=(ts>=w[0])&(ts<=w[1]);pred_s=float(np.polyfit(ts[mask],ell[mask],1)[0]);act_s=slope(r,w)
  rec=dict(n=n,U=U,W=W,f=f,normalization=norm,seed=r['seed'],predicted_slope=pred_s,actual_slope=act_s,
   relative_curve_error=float(np.linalg.norm(actual[1:]-ell[1:])/np.linalg.norm(ell[1:])),
   predicted_loss=(W*W*ell).tolist(),actual_loss=(W*W*actual).tolist())
  rows.append(rec);group[n,U,W,f,norm].append(rec)
 summary=[]
 for (n,U,W,f,norm),rr in sorted(group.items()):
  pred=np.array([r['predicted_slope'] for r in rr]);act=np.array([r['actual_slope'] for r in rr])
  summary.append(dict(n=n,U=U,W=W,f=f,normalization=norm,n_seeds=len(rr),
    predicted_mean=float(pred.mean()),actual_mean=float(act.mean()),
    ratio=float(act.mean()/pred.mean()) if abs(pred.mean())>1e-8 else None,
    residual=interval(act-pred),relative_curve_error=interval([r['relative_curve_error'] for r in rr])))
 out=dict(status='POST_CONFIRMATION_DIAGNOSTIC_NOT_NEW_PRIMARY_TEST',rows=rows,summary=summary,
    wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu,
    source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'src/finite_time.py']})
 (run/'summary.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
 for r in summary:
  if r['n']==256:print({k:v for k,v in r.items() if k not in ['residual','relative_curve_error']},'curve_error',r['relative_curve_error']['mean'])
 print('Saved',run)
if __name__=='__main__':main()
