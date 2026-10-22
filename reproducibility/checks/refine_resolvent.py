"""Post hoc eta/chain-size refinement of the three worst saved resolvent points.

This checks numerical resolution of the existing diagnostic, not new dynamics.
It deliberately reuses the existing banded-solve implementation. The original
25-point source and run are read-only; selection is by the saved absolute error.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, platform, sys, time, traceback
import numpy as np
import scipy
from verify_manuscript import kernel_resolvent

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    saved = ROOT/'checks/results/20261014_v1/summary.json'
    src = Path(__file__).resolve()
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(),
                    role='exploratory_resolution_check',seed=None,
                    selection='three largest absolute relative errors in the original 25-point run',
                    argv=sys.argv,cwd=str(Path.cwd()),
                    inputs={str(p.relative_to(ROOT)):sha(p) for p in (saved,src,src.with_name('verify_manuscript.py'))},
                    config=dict(etas=[2e-3,1e-3,5e-4,2.5e-4],half_chain_lengths=[150000,300000]),
                    environment=dict(executable=sys.executable,python=sys.version,platform=platform.platform(),
                                     numpy=np.__version__,scipy=scipy.__version__),status='RUNNING')
    wall,cpu = time.perf_counter(),time.process_time()
    code = 0
    try:
        old = json.loads(saved.read_text())['resolvent_check']['points']
        worst = sorted(old,key=lambda x:abs(x['rel']),reverse=True)[:3]
        rows=[]
        for point in worst:
            values={}
            for L in manifest['config']['half_chain_lengths']:
                values[L]=[float(kernel_resolvent(point['K'],point['U'],point['q'],L,eta))
                           for eta in manifest['config']['etas']]
            v=values[300000]
            extrap=[2*v[i+1]-v[i] for i in range(3)]
            relative=[(x-point['formula'])/point['formula'] for x in extrap]
            size_error=max(abs(a-b)/point['formula'] for a,b in zip(values[150000],v))
            rows.append(dict(point=point,raw_by_half_chain=values,extrapolated=extrap,
                             relative_errors=relative,max_chain_size_relative_difference=size_error,
                             error_reduction_factors=[abs(relative[i]/relative[i+1]) for i in range(2)]))
        max_final=max(abs(x['relative_errors'][-1]) for x in rows)
        max_size=max(x['max_chain_size_relative_difference'] for x in rows)
        thresholds=dict(final_extrapolation_max_relative=5e-6,chain_size_max_relative=2e-8)
        if not (max_final<thresholds['final_extrapolation_max_relative'] and
                max_size<thresholds['chain_size_max_relative'] and
                all(3<r<5 for row in rows for r in row['error_reduction_factors'])):
            raise AssertionError('The planned resolution criteria are not met; inspect saved rows.')
        summary=dict(rows=rows,max_final_relative_error=max_final,max_chain_size_relative_difference=max_size,
                     thresholds=thresholds,scope='Post hoc numerical resolution, three selected points only')
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
        manifest['status']='PASS'
        print(json.dumps(summary,indent=2))
    except Exception:
        code=1
        manifest['status']='FAIL'
        manifest['traceback']=traceback.format_exc()
        print(manifest['traceback'],file=sys.stderr)
    finally:
        manifest.update(finished_utc=datetime.now(timezone.utc).isoformat(),exit_code=code,
                        wall_seconds=time.perf_counter()-wall,process_cpu_seconds=time.process_time()-cpu)
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    sys.exit(code)

if __name__=='__main__':
    main()
