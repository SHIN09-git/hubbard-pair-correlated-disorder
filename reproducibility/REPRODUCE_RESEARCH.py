"""Clean-source replay of theory, frozen confirmation, leakage and figures.

Creates an isolated source snapshot and all results underneath a new directory.
No downloads, no overwrite of prior runs, no public uploads.
"""
from pathlib import Path
from datetime import datetime,timezone
from uuid import uuid4
import subprocess,sys,os,time,json,hashlib,platform,shutil

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 root=Path(__file__).resolve().parent
 run=root/'results/full_reproductions'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid4().hex[:8])
 snap=run/'source';snap.mkdir(parents=True)
 files=list(root.glob('*.py'))+list(root.glob('requirements*.txt'))+[root/'CONFIRMATION_PROTOCOL.md',root/'notes/G3_LEAKAGE_PROTOCOL.md']
 for folder in ['src','tools','tests','configs']:
  files += [p for p in (root/folder).glob('*') if p.is_file() and p.suffix in ['.py','.json']]
 sources={}
 for p in files:
  dest=snap/p.relative_to(root);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest);sources[str(p.relative_to(root))]=sha(p)
 env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONIOENCODING='utf-8')
 commands=[('tests',['-m','unittest','discover','-s','tests','-v']),('theory',['tools/check_closed_kernel.py']),
  ('windows',['run_window_pilot.py']),('confirmation_batch1',['run_confirmation.py','configs/confirmation_v1.json']),
  ('confirmation_batch2',['run_confirmation.py','configs/confirmation_v1_batch2.json']),('analysis',['analyze_confirmation.py']),
  ('sample_diagnostic',['analyze_sample_response.py']),('leakage',['run_leakage.py']),('figures',['make_figures.py'])]
 manifest=dict(status='RUNNING',source_sha256=sources,python=sys.version,platform=platform.platform(),commands=[],
   note='Replay uses the same seeds; it adds reproducibility evidence, not independent samples. Child CPU costs are in each scientific step manifest.')
 (run/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 start=time.perf_counter()
 for label,args in commands:
  t=time.perf_counter();cmd=[sys.executable]+args
  p=subprocess.run(cmd,cwd=snap,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8')
  (run/(label+'.log')).write_text(p.stdout,encoding='utf-8')
  manifest['commands'].append(dict(label=label,command=cmd,returncode=p.returncode,wall_seconds=time.perf_counter()-t))
  print(label,'PASS' if p.returncode==0 else 'FAIL',round(time.perf_counter()-t,2),flush=True)
  if p.returncode:
   manifest['status']='FAIL';break
 else:manifest['status']='PASS'
 manifest['wall_seconds']=time.perf_counter()-start
 manifest['output_sha256']={str(p.relative_to(run)):sha(p) for p in (snap/'results').rglob('*') if p.is_file() and p.suffix in ['.json','.npz','.log']}
 (run/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 print(manifest['status'],run,flush=True)
 return 0 if manifest['status']=='PASS' else 1
if __name__=='__main__':sys.exit(main())
