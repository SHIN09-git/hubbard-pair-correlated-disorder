"""Kinematic support only. Does not compute a T-matrix or a dissociation rate."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,math,platform,sys,time
import numpy as np


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter();cpu=time.process_time();pi=math.pi
    J=1.;U=3.5;k=.75*pi;E=math.sqrt(U*U+16*J*J*math.cos(k/2)**2);Q=2*math.acos(E/(4*J))
    odd_gap=pi-k-Q;even_gap=k-Q
    cases=[]
    for f in [.95,.97,.99]:
        delta=(1-f)*pi
        possible=[n for n in range(1,101) if n*delta>(odd_gap if n%2 else even_gap)]
        first=possible[0]
        cases.append(dict(q0_over_pi=f,first_kinematically_allowed_amplitude_order=first,
            conditional_candidate_rate_power=2*first,nonzero_amplitude_proven=False))
    assert [r['first_kinematically_allowed_amplitude_order'] for r in cases]==[3,3,7]
    N=65536;q=2*pi*np.fft.fftfreq(N);wrap=lambda z:(z+pi)%(2*pi)-pi
    a=(np.abs(q)>.97*pi).astype(float);a/=a.sum();fa=np.fft.fft(a)
    omega=abs(wrap(k+q))<Q
    convolution=[]
    for n in [1,2,3]:
        dist=np.fft.ifft(fa**n).real;mass=float(dist[omega].sum())
        convolution.append(dict(insertions=n,mass_in_clean_open_arc=mass,minimum_roundoff=float(dist.min())))
        assert (abs(mass)<1e-12) if n<3 else (mass>0.01)
    data=dict(status='PASS',scope='Support geometry and cyclic convolution; no amplitude or rate calculation.',
        J=J,U=U,K_over_pi=k/pi,E=E,Q_over_pi=Q/pi,
        Omega_over_pi=[(-k-Q)/pi,(-k+Q)/pi],
        odd_positive_measure_threshold_over_pi=odd_gap/pi,
        even_positive_measure_threshold_over_pi=even_gap/pi,
        cases=cases,grid=N,convolution=convolution,
        equality_rule='Strict > for positive-measure overlap; equality is endpoint contact, established analytically, not by the floating-point grid.')
    result=args.output/'support_check.json';result.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    metadata=dict(status='PASS',created_utc=datetime.now(timezone.utc).isoformat(),python=sys.version,numpy=np.__version__,
        platform=platform.platform(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        output_sha256=hashlib.sha256(result.read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-start,
        process_cpu_seconds=time.process_time()-cpu,exit_code=0)
    (args.output/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False))


if __name__=='__main__':main()
