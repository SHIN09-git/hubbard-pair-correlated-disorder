"""Apply pre-frozen windows, paired comparisons and descriptive intervals."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import numpy as np
from scipy.stats import t as student
import json,hashlib,time
from src.physics import born_rate

ROOT=Path(__file__).resolve().parent
def interval(x):
 x=np.asarray(x);n=len(x);mean=float(x.mean());sd=float(x.std(ddof=1));half=float(student.ppf(.975,n-1)*sd/np.sqrt(n))
 return dict(n=n,mean=mean,sd=sd,ci95=[mean-half,mean+half])
def slope(r,window):
 ts=np.array(r['times']);loss=1-np.array(r['bound_populations']);m=(ts>=window[0])&(ts<=window[1])
 return float(np.polyfit(ts[m],loss[m],1)[0]/r['W']**2)
def main():
 start=time.perf_counter();cpu=time.process_time();groups=defaultdict(list);inputs={}
 runs=[]
 for run in sorted((ROOT/'results/confirmation').iterdir()):
  if not run.is_dir():continue
  manifest=json.loads((run/'manifest.json').read_text())
  if manifest['status']!='PASS':continue
  runs.append(str(run))
  for p in run.glob('N*.json'):
   inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
   r=json.loads(p.read_text());key=(r['n'],r['U'],r['W'],r['f'],r['normalization']);groups[key].append(r)
 summary=[];lookup={};rowlookup={}
 pilot=sorted((ROOT/'results/window_pilot').glob('*/manifest.json'))[-1].parent
 for key,rs in sorted(groups.items()):
  rs.sort(key=lambda r:r['seed']);assert len(set(r['seed'] for r in rs))==len(rs)
  n,U,W,f,norm=key;windows=[(12,24),(8,20),(16,28)] if U==2 else [(16,28),(12,24),(20,32)]
  rate=born_rate(.75*np.pi,U,f*np.pi,fixed_variance=norm=='variance')[0]
  ts=np.array(rs[0]['times']);curves=np.array([1-np.array(r['bound_populations']) for r in rs]);se=curves.std(axis=0,ddof=1)/np.sqrt(len(rs))
  pr=json.loads((pilot/f'N{n}_U{U}.json').read_text());pred=next(r for r in pr['rows'] if r['q0_over_pi']==f)
  pt=np.array(pr['times']);pl=np.array(pred['fixed_variance' if norm=='variance' else 'fixed_passband'])
  sp={f'{a}-{b}':float(np.polyfit(pt[(pt>=a)&(pt<=b)],pl[(pt>=a)&(pt<=b)],1)[0]) for a,b in windows}
  row=dict(n=n,U=U,W=W,f=f,normalization=norm,seeds=[r['seed'] for r in rs],
    coefficient=rate,main_window=list(windows[0]),main=interval([slope(r,windows[0]) for r in rs]),
    sensitivity={f'{a}-{b}':interval([slope(r,(a,b)) for r in rs]) for a,b in windows},
    finite_time_prediction_slopes=sp,times=ts.tolist(),mean_loss=curves.mean(axis=0).tolist(),sem_loss=se.tolist(),
    finite_time_prediction_loss=(W*W*np.interp(ts,pt,pl)).tolist(),
    max_mean_loss=float(curves.mean(axis=0).max()),allowed_modes=rs[0]['allowed_modes'],
    expected_variance=rs[0]['expected_variance'],sample_variance=interval([r['sample_variance'] for r in rs]),
    mean_loss_at_24=interval(curves[:,list(ts).index(24.)]))
  summary.append(row);lookup[key]=rs;rowlookup[key]=row
 comparisons=[];rng=np.random.default_rng(746205)
 for n in [128,256]:
  for W in [.02,.04,.08]:
   k0=(n,2.,W,0.,'variance')
   if k0 not in lookup:continue
   for norm in ['variance','passband']:
    k1=(n,2.,W,.8,norm)
    if k1 not in lookup:continue
    base={r['seed']:r for r in lookup[k0]};hi={r['seed']:r for r in lookup[k1]};seeds=sorted(set(base)&set(hi))
    x=np.array([slope(base[s],(12,24)) for s in seeds]);y=np.array([slope(hi[s],(12,24)) for s in seeds])
    ix=rng.integers(0,len(x),size=(10000,len(x)));rr=y[ix].mean(axis=1)/x[ix].mean(axis=1)
    comparisons.append(dict(n=n,W=W,normalization=norm,seeds=seeds,difference=interval(y-x),ratio=float(y.mean()/x.mean()),
      bootstrap_ratio_ci95=np.quantile(rr,[.025,.975]).tolist(),theory_ratio=born_rate(.75*np.pi,2,.8*np.pi,fixed_variance=norm=='variance')[0]/born_rate(.75*np.pi,2,0.)[0]))
 run=ROOT/'results/analysis'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');run.mkdir(parents=True)
 out=dict(status='APPLIED_FROZEN_PROTOCOL',input_runs=runs,groups=summary,comparisons=comparisons,
  input_sha256=inputs,analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu)
 (run/'summary.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
 for c in comparisons:print('PAIR',c['n'],c['W'],c['normalization'],'diff',c['difference'],'ratio',c['ratio'],c['bootstrap_ratio_ci95'])
 for r in summary:
  if r['n']==256:print('GROUP',r['U'],r['W'],r['f'],r['normalization'],'slope',r['main'],'c',r['coefficient'],'loss24',r['mean_loss_at_24']['mean'])
 print('Saved',run)
if __name__=='__main__':main()
