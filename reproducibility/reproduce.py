"""Run the selected saved-data reproduction in a new directory, without dynamics."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ctypes
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def child_cpu(process):
    if os.name != 'nt':
        return None
    from ctypes import wintypes
    values = [wintypes.FILETIME() for _ in range(4)]
    fn = ctypes.WinDLL('kernel32', use_last_error=True).GetProcessTimes
    fn.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
    fn.restype = wintypes.BOOL
    if not fn(int(process._handle), *[ctypes.byref(x) for x in values]):
        return None
    return sum((x.dwHighDateTime << 32) + x.dwLowDateTime for x in values[2:]) / 1e7


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', choices=['verify', 'tables', 'posthoc', 'checks', 'figures', 'all'], default='all')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    wall, cpu = time.perf_counter(), time.process_time()
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(), argv=sys.argv,
                    target=args.target, seed=None, source_sha256=sha(Path(__file__)),
                    python=sys.version, executable=sys.executable, platform=platform.platform(),
                    scope='Saved-data reconstruction and bounded structural checks; no new trajectories.',
                    seeds='Original seeds/configurations and deterministic check seeds remain in input files.',
                    steps=[], status='RUNNING')
    record = out / 'manifest.json'
    code = 0
    try:
        saved = json.loads((ROOT / 'MANIFEST.json').read_text(encoding='utf-8'))
        mismatches = [name for name, digest in saved['files'].items()
                      if not (ROOT / name).is_file() or sha(ROOT / name) != digest]
        manifest['input_manifest_sha256'] = sha(ROOT / 'MANIFEST.json')
        manifest['files_verified'] = len(saved['files'])
        manifest['input_mismatches'] = mismatches
        if mismatches:
            raise ValueError('Input checksum mismatch: ' + repr(mismatches[:5]))
        env = os.environ.copy()
        env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONIOENCODING='utf-8')

        def run(name, command, cwd=ROOT):
            start = time.perf_counter()
            with (out / (name + '.log')).open('wb') as log:
                process = subprocess.Popen(command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
                status = process.wait()
                process_cpu = child_cpu(process)
            manifest['steps'].append(dict(name=name, argv=command, exit_code=status,
                                          wall_seconds=time.perf_counter() - start,
                                          process_cpu_seconds=process_cpu,
                                          cpu_scope='Direct child only; not a total including descendants.',
                                          log=name + '.log'))
            record.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
            print(name, 'PASS' if status == 0 else 'FAIL', flush=True)
            if status:
                raise RuntimeError(name + ' exited ' + str(status))

        targets = {'tables', 'posthoc', 'checks', 'figures'} if args.target == 'all' else {args.target}
        if 'tables' in targets:
            run('tables', [sys.executable, 'REBUILD_FROM_DATA.py', '--output', str(out / 'tables')])
            table_record = json.loads((out / 'tables/manifest.json').read_text(encoding='utf-8'))
            if table_record['status'] != 'PASS':
                raise AssertionError('Saved table/endpoint reconstruction did not match')
            manifest['tables_comparison'] = table_record['numeric_comparison']
        if 'posthoc' in targets:
            run('posthoc', [sys.executable, 'tools/reanalyze_existing_data.py', '--input-root', str(ROOT),
                            '--output', str(out / 'posthoc')])
            from REBUILD_FROM_DATA import differences
            errors = differences(json.loads((ROOT / 'results/reanalysis/summary.json').read_text('utf-8')),
                                 json.loads((out / 'posthoc/summary.json').read_text('utf-8')))
            manifest['posthoc_comparison'] = dict(rtol=1e-10, atol=1e-12, mismatches=errors)
            if errors:
                raise AssertionError('Post hoc results differ: ' + repr(errors[:5]))
        if 'checks' in targets:
            for name, command in [
                ('tests', ['-m', 'unittest', 'discover', '-s', 'tests', '-v']),
                ('kernel_25', ['checks/verify_manuscript.py', '--output', str(out / 'kernel_25')]),
                ('refined_resolvent', ['checks/refine_resolvent.py', '--output', str(out / 'refined_resolvent')]),
                ('finite_time', ['tools/infinite_chain_finite_time.py', '--output', str(out / 'finite_time')]),
                ('higher_order_support', ['tools/check_momentum_support.py', '--output', str(out / 'higher_order_support')]),
            ]:
                run(name, [sys.executable] + command)
            # The original 25-point program's PASS only records successful execution.
            # Retain the actual maximum error and compare against the stated 7e-5 scale.
            kernel = json.loads((out / 'kernel_25/summary.json').read_text('utf-8'))['resolvent_check']['points']
            error = max(abs(row['rel']) for row in kernel)
            manifest['kernel_25_error'] = dict(points=len(kernel), maximum_relative=error, threshold=7e-5)
            if len(kernel) != 25 or error >= 7e-5:
                raise AssertionError('25-point resolvent result exceeds the reported scale')
        if 'figures' in targets:
            dest = out / 'figures'
            shutil.copytree(ROOT / 'figures', dest)
            run('figures', [sys.executable, str(dest / 'make_manuscript_figures.py')], cwd=out)
        manifest['status'] = 'PASS'
    except Exception as error:
        code = 1
        manifest.update(status='FAIL', error=repr(error))
    finally:
        manifest.update(finished_utc=datetime.now(timezone.utc).isoformat(), exit_code=code,
                        wall_seconds=time.perf_counter() - wall,
                        parent_process_cpu_seconds=time.process_time() - cpu,
                        output_sha256={p.relative_to(out).as_posix(): sha(p)
                                       for p in out.rglob('*') if p.is_file() and p != record})
        record.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': manifest['status'], 'record': str(record)}, indent=2))
    return code


if __name__ == '__main__':
    sys.exit(main())
