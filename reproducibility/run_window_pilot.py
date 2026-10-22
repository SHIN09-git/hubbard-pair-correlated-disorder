"""Deterministic time-window design before confirmation seeds are run."""
from pathlib import Path
from datetime import datetime,timezone
import json,time,hashlib,sys,platform
import numpy as np
from src.finite_time import finite_time_loss
from src.physics import born_rate

def main():
    start=time.perf_counter();cpu=time.process_time()
    root=Path(__file__).resolve().parent
    run=root/'results'/'window_pilot'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run.mkdir(parents=True)
    rows=[];summary=[]
    ts=np.arange(0.,65.,2.)
    for n in (64,128,256):
        for U in (2.,3.5,5.):
            tic=time.perf_counter()
            row=finite_time_loss(n,.75*np.pi,U,ts)
            rows.append(row)
            (run/f'N{n}_U{U}.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
            for r in row['rows']:
                cc=born_rate(.75*np.pi,U,r['q0_over_pi']*np.pi)[0]
                ss={f'{lo}-{hi}':float(np.polyfit(ts[(ts>=lo)&(ts<=hi)],np.array(r['fixed_variance'])[(ts>=lo)&(ts<=hi)],1)[0]) for lo,hi in [(8,20),(12,24),(16,28),(24,40),(32,48)]}
                summary.append(dict(n=n,U=U,f=r['q0_over_pi'],c=cc,slopes=ss))
            print(n,U,'seconds',time.perf_counter()-tic,flush=True)
    out=dict(status='DESIGN_ONLY_NO_CONFIRMATION_SEEDS',summary=summary,wall_seconds=time.perf_counter()-start,
             process_cpu_seconds=time.process_time()-cpu,python=sys.version,platform=platform.platform(),
             source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),root/'src/finite_time.py',root/'src/physics.py']})
    (run/'manifest.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print('Saved',run,flush=True)
if __name__=='__main__':main()
