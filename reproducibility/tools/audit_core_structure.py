"""New structural checks, independent occupation-number construction; no dynamics rerun."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,time,hashlib
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.physics import relative_block,two_body_hamiltonian

def fock(n,U,J,v):
    pairs=[(x,y) for x in range(n) for y in range(x,n)]
    occ=[]
    for x,y in pairs:
        row=np.zeros(n,dtype=int);row[x]+=1;row[y]+=1;occ.append(row)
    occ=np.array(occ);index={tuple(row):i for i,row in enumerate(occ)}
    h=np.diag(occ@v+U/2*np.sum(occ*(occ-1),axis=1)).astype(complex)
    T=np.zeros_like(h);W=np.zeros((n*n,len(pairs)),complex)
    for j,row in enumerate(occ):
        T[index[tuple(np.roll(row,1))],j]=1
        for source in range(n):
            if not row[source]:continue
            for target in [(source-1)%n,(source+1)%n]:
                changed=row.copy();changed[source]-=1;changed[target]+=1
                h[index[tuple(changed)],j]+=-J*np.sqrt(row[source]*(row[target]+1))
        x,y=pairs[j]
        if x==y:W[x*n+y,j]=1
        else:W[x*n+y,j]=W[y*n+x,j]=1/np.sqrt(2)
    return h,T,W,occ

def bloch(n,K):
    B=np.zeros((n*n,n),complex)
    for x in range(n):
        for r in range(n):B[x*n+(x+r)%n,r]=np.exp(1j*K*x)/np.sqrt(n)
    return B

def closed(q,K,U,J=1):
    E=np.sqrt(U*U+16*J*J*np.cos(K/2)**2);t=2*J*np.cos((K+q)/2)
    d2=4*t*t-E*E
    if d2<=0:return 0.
    B=1-4*np.cos(K/2)**2*t*t/(E*E)
    return 8*(U/E)**3*np.sqrt(d2)/((d2+U*U)*B*B)

def bound_parameters(K,U,J=1):
    t=2*J*np.cos(K/2);E=np.sqrt(U*U+4*t*t)
    a=-2*t/(E+U)
    return E,a,np.sqrt(U/E)

def elastic_amplitude(q,K,U):
    E,a,C=bound_parameters(K,U);Ep,ap,Cp=bound_parameters(K+q,U)
    rho=a*ap
    return 2*C*Cp*(1-rho*rho)/(1-2*rho*np.cos(q/2)+rho*rho)

def main():
    start=time.perf_counter();cpu=time.process_time()
    run=ROOT/'results/core_review'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    run.mkdir(parents=True)
    rng=np.random.default_rng(841725)
    records=[]
    for n in [7,8,9,10]:
        U,J=1.7,1.2;v=rng.normal(size=n)
        h,T,W,occ=fock(n,U,J,v)
        mapping=np.max(np.abs(W.conj().T@two_body_hamiltonian(v,U,J)@W-h))
        h0,_,_,_=fock(n,U,J,np.zeros(n))
        dim=[];source_error=0.;energy_error=0.;symmetry_error=0.
        for m in range(n):
            K=2*np.pi*m/n
            P=sum(np.exp(1j*K*l)*np.linalg.matrix_power(T,l) for l in range(n))/n
            vals,vec=eigh((P+P.conj().T)/2);Q=vec[:,vals>.5]
            expected=(n+1)//2 if n%2 else n//2+(m%2==0)
            assert Q.shape[1]==expected
            dim.append(expected)
            ef=eigh(Q.conj().T@h0@Q,eigvals_only=True)
            er,vr=eigh(relative_block(n,K,U,J));b=vr[:,-1];B=bloch(n,K)
            energy_error=max(energy_error,abs(ef[-1]-er[-1]))
            state=(B@b).reshape(n,n)
            symmetry_error=max(symmetry_error,np.linalg.norm(state-state.T))
            ff=W.conj().T@(B@b)
            assert abs(np.linalg.norm(ff)-1)<1e-11
            for iq in [0,1,n//2,n-1]:
                q=2*np.pi*iq/n
                dq=occ@np.exp(1j*q*np.arange(n))
                actual=bloch(n,K+q).conj().T@W@(dq*ff)
                wanted=(1+np.exp(1j*q*np.arange(n)))*b
                source_error=max(source_error,np.max(np.abs(actual-wanted)))
        assert sum(dim)==n*(n+1)//2
        assert max(mapping,energy_error,source_error,symmetry_error)<1e-11
        records.append(dict(n=n,sector_dimensions=dim,mapping_error=float(mapping),
            bound_energy_error=float(energy_error),source_error=float(source_error),
            exchange_error=float(symmetry_error)))
    # Spectral resolution of the identity: an absolute normalization check.
    nodes,weights=leggauss(800);p=(nodes+1)*np.pi/2;weights*=np.pi/2
    r=np.arange(-3,4);identity_even=(np.eye(7)+np.fliplr(np.eye(7)))/2
    spectral=[]
    for t in [1.3,-1.3]:
        for U in [.7,3.,7.]:
            E=np.sqrt(U*U+4*t*t);a=-2*t/(E+U);b=np.sqrt(U/E)*a**abs(r)
            D=2*abs(t)*np.sin(p)
            phi=(D[:,None]*np.cos(p[:,None]*r)+U*np.sign(t)*np.sin(p[:,None]*abs(r)))/np.sqrt(np.pi*(D[:,None]**2+U*U))
            residual=(phi.T*weights)@phi+np.outer(b,b)-identity_even
            err=float(np.max(abs(residual)));assert err<2e-11
            spectral.append(dict(t=t,U=U,completeness_error=err))
    # Folding covariance in the source and in the full retarded spectral measure.
    folding=[]
    for K,U,q in [(2.4,.9,2.8),(-2.2,1.3,-2.9),(1.6,2.,-2.2)]:
        E,a,C=bound_parameters(K,U);r=np.arange(-70,71)
        t=2*np.cos((K+q)/2);D=np.sqrt(4*t*t-E*E)
        assert D>0
        z=(-E+1j*D)/(2*t);g=z**abs(r)/(1j*D)
        gmat=z**abs(r[:,None]-r[None,:])/(1j*D)
        s=2*C*a**abs(r)*np.cos(q*r/2)
        # Direct convolution of the free resolvent and a rank-one Dyson update.
        value=-2*np.imag(s@gmat@s+U/(1-U/(1j*D))*(s@g)**2)
        expected=closed(q,K,U)
        err=float(abs(value-expected)/expected)
        period_error=float(abs(closed(q+2*np.pi,K,U)-expected)/expected)
        assert max(err,period_error)<1e-9
        folding.append(dict(K=K,U=U,q=q,resolvent=value,closed=expected,relative_error=err,fold_error=period_error))
    # A hard spectral jump is tested by the finite-time approximate identity.
    K=.75*np.pi;U=2.;qe=2*np.pi-2*K;E,_,_=bound_parameters(K,U)
    velocity=4*abs(np.sin(K))/E;target=elastic_amplitude(qe,K,U)**2/velocity
    edge=[]
    for T in [100,500,2000,10000]:
        def integrand(q):
            d=bound_parameters(K+q,U)[0]-E
            return elastic_amplitude(q,K,U)**2*T*np.sinc(T*d/(2*np.pi))**2/(2*np.pi)
        left=quad(integrand,qe-.2,qe,epsabs=1e-9,limit=4000)[0]
        right=quad(integrand,qe,qe+.2,epsabs=1e-9,limit=4000)[0]
        edge.append(dict(time=T,full_over_on_shell=(left+right)/target,
                         cutoff_equals_qB_over_on_shell=right/target))
    assert abs(edge[-1]['cutoff_equals_qB_over_on_shell']-.5)<.002
    assert abs(edge[-1]['full_over_on_shell']-1)<.002
    # Narrow U peak: integrate in scaled outgoing momentum, not on a coarse q grid.
    weak=[];k=.75*np.pi;s=np.sin(k/2);c=np.cos(k/2)
    for U in [.08,.04,.02]:
        E=np.sqrt(U*U+16*c*c);Q=2*np.arccos(E/4)
        value=quad(lambda y:closed(U*y/(2*c)-k,k,U)*U/(2*c)/(2*np.pi),
                   -2*c*Q/U,2*c*Q/U,epsabs=1e-10,limit=500)[0]
        weak.append(dict(U=U,white_coefficient=value,limit=1/s,ratio=value*s))
    assert abs(weak[-1]['ratio']-1)<.002
    result=dict(status='PASS',description='Independent structure checks; no confirmation trajectories rerun.',
        occupation_number_checks=records,continuum_completeness=spectral,folding_resolvent=folding,
        elastic_hard_edge=edge,weak_U_limit=weak,process_cpu_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-start,python=sys.version,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__),ROOT/'src/physics.py']})
    (run/'checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2));print('OUTPUT',run)

if __name__=='__main__':main()
