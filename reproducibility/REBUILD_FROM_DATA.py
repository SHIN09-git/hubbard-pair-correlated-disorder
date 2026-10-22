"""Portable default: regenerate tables/figures and structure checks, no dynamics."""
from pathlib import Path
from datetime import datetime,timezone
from uuid import uuid4
import argparse, shutil, subprocess, sys, os, json, hashlib, time, math, platform

def differences(a,b,path=''):
    if isinstance(a,dict) and isinstance(b,dict):
        if a.keys()!=b.keys(): return [path+': keys differ']
        return sum((differences(a[k],b[k],path+'/'+str(k)) for k in a),[])
    if isinstance(a,list) and isinstance(b,list):
        if len(a)!=len(b):return [path+': length differs']
        return sum((differences(x,y,path+'/'+str(i)) for i,(x,y) in enumerate(zip(a,b))),[])
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        return [] if math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12) else [path+': numeric mismatch']
    return [] if a==b else [path+': mismatch']

def main():
    root=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    out=(args.output or root/'rebuilds'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid4().hex[:8])).resolve()
    if out.exists():raise FileExistsError('Output must be a new directory: '+str(out))
    snap=out/'source';snap.mkdir(parents=True)
    files=[root/n for n in ['analyze_confirmation.py','make_figures.py','tools/check_closed_kernel.py',
               'tools/audit_core_structure.py','tools/export_figure_data.py']]
    files+=list((root/'src').glob('*.py'))
    for folder in ['confirmation','window_pilot','sample_diagnostic','leakage']:
        files += [p for p in (root/'results'/folder).rglob('*') if p.is_file()]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    inputs={}
    for p in files:
        target=snap/p.relative_to(root);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,target);inputs[p.relative_to(root).as_posix()]=sha(p)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONIOENCODING='utf-8')
    commands=['tools/check_closed_kernel.py','analyze_confirmation.py','make_figures.py',
              'tools/export_figure_data.py','tools/audit_core_structure.py']
    manifest=dict(status='RUNNING',python=sys.version,platform=platform.platform(),
             input_sha256=inputs,commands=[],note='Saved-data regeneration only. No confirmation or leakage dynamics invoked.')
    start=time.perf_counter()
    for cmd in commands:
        t=time.perf_counter()
        proc=subprocess.run([sys.executable,cmd],cwd=snap,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8')
        (out/(Path(cmd).stem+'.log')).write_text(proc.stdout,encoding='utf-8')
        manifest['commands'].append(dict(script=cmd,returncode=proc.returncode,wall_seconds=time.perf_counter()-t))
        print(cmd,'PASS' if proc.returncode==0 else 'FAIL',flush=True)
        if proc.returncode:
            manifest['status']='FAIL';break
    else:
        original=json.loads(sorted(root.glob('results/analysis/*/summary.json'))[-1].read_text())
        rebuilt=json.loads(sorted(snap.glob('results/analysis/*/summary.json'))[-1].read_text())
        delta=[]
        for key in ['groups','comparisons']:
            delta += differences(original[key],rebuilt[key],key)
        oldedge=json.loads((root/'results/theory_audit/closed_kernel_initial.json').read_text())['edges']
        newedge=json.loads((snap/'results/theory_audit/closed_kernel_initial.json').read_text())['edges']
        delta += differences(oldedge,newedge,'edges')
        manifest['numeric_comparison']=dict(groups=len(original['groups']),paired_comparisons=len(original['comparisons']),
            endpoint_branches=len(oldedge),rtol=1e-10,atol=1e-12,mismatches=delta)
        audit=json.loads(sorted(snap.glob('results/core_review/*/checks.json'))[-1].read_text())
        manifest['structure_check_status']=audit['status']
        manifest['status']='PASS' if not delta and audit['status']=='PASS' else 'FAIL'
    manifest['wall_seconds']=time.perf_counter()-start
    manifest['output_sha256']={p.relative_to(snap).as_posix():sha(p) for folder in ['docs','data','results/analysis','results/core_review','results/theory_audit'] for p in (snap/folder).rglob('*') if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(manifest['status'],out,flush=True)
    return 0 if manifest['status']=='PASS' else 1
if __name__=='__main__':sys.exit(main())
