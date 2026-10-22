"""O(W^2) finite-time loss, including off-shell continuum transitions."""
import numpy as np
from scipy.linalg import eigh
from .physics import relative_block

def finite_time_loss(n,K,U,ts,fractions=(0.,.8,.9,.97),return_modes=False):
    e0,v0=eigh(relative_block(n,K,U));Eb=e0[-1];b=v0[:,-1]
    qs=2*np.pi*np.fft.fftfreq(n);r=np.arange(n)
    kernels=np.zeros((n,len(ts)))
    for iq,q in enumerate(qs):
        e,v=eigh(relative_block(n,K+q,U))
        weights=np.abs(v[:,:-1].conj().T@((1+np.exp(1j*q*r))*b))**2
        de=Eb-e[:-1]
        transition=np.asarray(ts)[None,:]**2*np.sinc(de[:,None]*np.asarray(ts)[None,:]/(2*np.pi))**2
        kernels[iq]=weights@transition
    rows=[]
    for f in fractions:
        mask=np.abs(qs)>f*np.pi+1e-12
        ell=kernels[mask].mean(axis=0)
        rows.append(dict(q0_over_pi=f,allowed_modes=int(mask.sum()),finite_fraction=float(mask.mean()),
                         fixed_variance=ell.tolist(),fixed_passband=(ell*mask.mean()).tolist()))
    out=dict(n=n,U=U,K_over_pi=K/np.pi,times=np.asarray(ts).tolist(),rows=rows)
    if return_modes:out['mode_kernels']=kernels.tolist()
    return out
