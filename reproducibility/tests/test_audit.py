"""Differential checks for the new scientific derivation and finite-time route."""
import unittest
import numpy as np
from scipy.sparse.linalg import expm_multiply
from src.physics import dissociation_kernel,born_rate,clean_bound_basis,initial_bound_state,bound_population,two_body_hamiltonian
from src.finite_time import finite_time_loss
from tools.check_closed_kernel import closed_kernel

class AuditTests(unittest.TestCase):
 def test_closed_kernel_differential_general_J_and_negative_K(self):
  rng=np.random.default_rng(871325)
  count=0
  for _ in range(500):
   J=rng.uniform(.4,2.);K=rng.uniform(-np.pi,np.pi);U=rng.uniform(.05,3.9)*J;q=rng.uniform(-np.pi,np.pi)
   x=dissociation_kernel(q,K,U,J);y=closed_kernel(q,K,U,J)
   np.testing.assert_allclose(x,y,rtol=2e-9,atol=1e-12)
   if x>0:count+=1
  self.assertGreater(count,100)
 def test_all_open_outgoing_momenta_have_positive_matrix_elements(self):
  for K,U in [(.75*np.pi,2.),(.75*np.pi,3.5),(-.75*np.pi,3.5),(np.pi,2.)]:
   E=np.sqrt(U*U+16*np.cos(K/2)**2);Q=2*np.arccos(E/4)
   for kp in np.linspace(-.999*Q,.999*Q,29):
    q=(kp-K+np.pi)%(2*np.pi)-np.pi
    self.assertGreater(closed_kernel(q,K,U),0)
 def test_finite_time_spectral_prediction_against_full_realspace(self):
  n=16;K=.75*np.pi;U=2.;times=np.array([0.,2.,4.]);eps=3e-4
  b=clean_bound_basis(n,U);psi=initial_bound_state(b,3*n//8)
  pred=finite_time_loss(n,K,U,times,fractions=(0.,.8))
  for row in pred['rows']:
   f=row['q0_over_pi'];mask=np.abs(2*np.pi*np.fft.fftfreq(n))>f*np.pi+1e-12
   raw=np.zeros(n);raw[0]=np.sqrt(n)
   v=np.fft.ifft(np.fft.fft(raw)*mask).real/np.sqrt(mask.mean())
   losses=[]
   for sign in [-1,1]:
    h=two_body_hamiltonian(sign*eps*v,U)
    states=expm_multiply(-1j*h,psi,start=0.,stop=4.,num=3,traceA=-1j*h.diagonal().sum())
    losses.append(np.array([1-bound_population(p,b) for p in states])/eps**2)
   # Equal Fourier powers give the exact ensemble coefficient at second order.
   np.testing.assert_allclose(np.mean(losses,axis=0)[1:],np.array(row['fixed_variance'])[1:],rtol=8e-5,atol=1e-6)
 def test_critical_endpoint_is_half_not_three_halves_at_fixed_variance(self):
  K=.75*np.pi;U=4*np.sqrt(-np.cos(K));ds=[1e-4,1e-5]
  rs=[born_rate(K,U,np.pi-d,epsabs=1e-14)[0] for d in ds]
  slope=np.log(rs[0]/rs[1])/np.log(ds[0]/ds[1])
  self.assertAlmostEqual(slope,.5,places=3)

if __name__=='__main__':unittest.main()
