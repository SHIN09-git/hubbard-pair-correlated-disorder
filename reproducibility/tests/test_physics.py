import unittest
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
from src.physics import (bound_energy,cutoff_boundary,elastic_backscatter_transfer,
    dissociation_kernel,born_rate,relative_block,clean_bound_basis,initial_bound_state,
    bound_population,two_body_hamiltonian,filtered_noise,finite_born_coefficients)


def independent_resolvent_kernel(q,K,U,J=1.,R=60):
    E=bound_energy(K,U,J); t=2*J*np.cos((K+q)/2)
    if E>=2*abs(t): return 0.
    D=np.sqrt(4*t*t-E*E); lam=(-E+1j*D)/(2*t)
    r=np.arange(-R,R+1); a=-4*J*np.cos(K/2)/(E+U)
    C=np.sqrt((1-a*a)/(1+a*a)); s=2*C*a**np.abs(r)*np.cos(q*r/2)
    g=lam**np.abs(r)/(1j*D)
    G=lam**np.abs(r[:,None]-r[None,:])/(1j*D)
    z=s@G@s+U/(1-U/(1j*D))*(s@g)**2
    return float(-2*z.imag)

class PhysicsTests(unittest.TestCase):
    def test_01_relative_hermiticity(self):
        h=relative_block(32,.75*np.pi,2.)
        np.testing.assert_allclose(h,h.conj().T,atol=1e-14)
    def test_02_bound_energy(self):
        for U in (2.,3.5,5.):
            self.assertAlmostEqual(eigh(relative_block(64,.75*np.pi,U),eigvals_only=True)[-1],bound_energy(.75*np.pi,U),places=11)
    def test_03_full_realspace_eigenstate(self):
        n=32; U=2.; basis=clean_bound_basis(n,U);p=initial_bound_state(basis,12)
        h=two_body_hamiltonian(np.zeros(n),U)
        self.assertLess(np.linalg.norm(h@p-bound_energy(.75*np.pi,U)*p),1e-11)
    def test_04_bosonic_symmetry(self):
        n=32;p=initial_bound_state(clean_bound_basis(n,2.),12).reshape(n,n)
        np.testing.assert_allclose(p,p.T,atol=1e-12)
    def test_05_bound_band_projection(self):
        b=clean_bound_basis(32,2.)
        for m in (0,5,12,16,31):
            self.assertAlmostEqual(bound_population(initial_bound_state(b,m),b),1.,places=12)
    def test_06_projection_not_double_occupation(self):
        n=32;p=initial_bound_state(clean_bound_basis(n,2.),12).reshape(n,n)
        onsite=np.sum(abs(np.diag(p))**2)
        self.assertAlmostEqual(onsite,2./bound_energy(.75*np.pi,2.),places=10)
        self.assertLess(onsite,.9)
    def test_07_kernel_resolvent_agreement(self):
        for U,K,q in ((2.,.75*np.pi,.85*np.pi),(2.,.75*np.pi,-.75*np.pi),(3.5,.75*np.pi,-.82*np.pi)):
            a=dissociation_kernel(q,K,U);b=independent_resolvent_kernel(q,K,U)
            self.assertGreater(a,0.)
            np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-12)
    def test_08_kernel_nonnegative(self):
        for q in np.linspace(-np.pi,np.pi,151):
            self.assertGreaterEqual(dissociation_kernel(q,.75*np.pi,2.),0.)
    def test_09_energy_closed_channel(self):
        self.assertTrue(cutoff_boundary(.75*np.pi,5.)['energy_closed'])
        self.assertEqual(born_rate(.75*np.pi,5.,.8*np.pi),(0.,0.))
    def test_10_interior_cutoff_boundary(self):
        b=cutoff_boundary(.75*np.pi,3.5)
        self.assertAlmostEqual(b['qmax']/np.pi,.9416570958029887,places=12)
        self.assertEqual(born_rate(.75*np.pi,3.5,.97*np.pi),(0.,0.))
        self.assertGreater(born_rate(.75*np.pi,3.5,.9*np.pi)[0],0.)
    def test_11_no_highpass_closure_U2(self):
        self.assertEqual(cutoff_boundary(.75*np.pi,2.)['qmax'],np.pi)
        self.assertGreater(born_rate(.75*np.pi,2.,.99*np.pi)[0],0.)
    def test_12_backscatter_q(self):
        self.assertAlmostEqual(elastic_backscatter_transfer(.75*np.pi),.5*np.pi,places=12)
    def test_13_normalization_control(self):
        c=born_rate(.75*np.pi,2.,.8*np.pi)[0]
        raw=born_rate(.75*np.pi,2.,.8*np.pi,fixed_variance=False)[0]
        self.assertAlmostEqual(raw/c,.2,places=12)
        white=born_rate(.75*np.pi,2.,0.)[0]
        self.assertGreater(c,white);self.assertLess(raw,white)
    def test_14_time_reversal(self):
        self.assertAlmostEqual(born_rate(.75*np.pi,2.,.8*np.pi)[0],born_rate(-.75*np.pi,2.,.8*np.pi)[0],places=10)
    def test_15_noise_spectral_gap(self):
        n=128;v=filtered_noise(n,.8*np.pi,.04,17)
        f=np.fft.fft(v);q=2*np.pi*np.fft.fftfreq(n)
        self.assertLess(np.max(abs(f[abs(q)<=.8*np.pi])),1e-12)
        self.assertAlmostEqual(v.mean(),0.,places=14)
    def test_16_noise_ensemble_variance(self):
        for frac in (0.,.8):
            var=np.mean([np.mean(filtered_noise(128,frac*np.pi,.04,s)**2) for s in range(300)])
            self.assertLess(abs(var/.04**2-1),.08)
    def test_17_unitarity(self):
        n=16;U=2.;b=clean_bound_basis(n,U);p=initial_bound_state(b,6)
        h=two_body_hamiltonian(filtered_noise(n,.5*np.pi,.04,7),U)
        pt=expm_multiply(-3j*h,p,traceA=-3j*h.diagonal().sum())
        self.assertAlmostEqual(np.linalg.norm(pt),1.,places=11)
        self.assertTrue(-1e-12<=bound_population(pt,b)<=1+1e-12)
    def test_18_invalid_inputs(self):
        for fn in (lambda:bound_energy(4.,2.),lambda:bound_energy(0.,0.),lambda:born_rate(0.,2.,np.pi),lambda:filtered_noise(2,0.,.1,0)):
            with self.assertRaises(ValueError):fn()
    def test_19_finite_eigensum_check(self):
        out=finite_born_coefficients(128,.75*np.pi,2.,etas=(.08,))
        for row in out['rows']:
            c=born_rate(.75*np.pi,2.,row['q0_over_pi']*np.pi)[0]
            self.assertLess(abs(row['coefficients'][0]/c-1),.03)
    def test_20_clean_dynamics_preserves_band(self):
        n=16;b=clean_bound_basis(n,2.);p=initial_bound_state(b,6)
        h=two_body_hamiltonian(np.zeros(n),2.)
        pt=expm_multiply(-10j*h,p,traceA=-10j*h.diagonal().sum())
        self.assertAlmostEqual(bound_population(pt,b),1.,places=11)

if __name__=='__main__':unittest.main(verbosity=2)
