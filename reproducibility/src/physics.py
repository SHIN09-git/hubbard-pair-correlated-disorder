"""Two distinguishable particles with common on-site disorder; symmetric sector.

Units: hbar = lattice spacing = 1; J is the hopping energy. The infinite-chain
routine returns the leading Born/Fermi-golden-rule coefficient gamma/W**2,
NOT a finite-W lifetime or an all-orders protection theorem.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy.sparse import diags, eye, kron, csr_matrix


def _check(K: float, U: float, J: float = 1.0) -> None:
    if not np.isfinite([K, U, J]).all() or U <= 0 or J <= 0:
        raise ValueError("K must be finite, U>0 and J>0 are required.")
    if not -np.pi <= K <= np.pi:
        raise ValueError("Use K in [-pi, pi].")


def bound_energy(K: float, U: float, J: float = 1.0) -> float:
    _check(K, U, J)
    return float(np.sqrt(U * U + 16 * J * J * np.cos(K / 2) ** 2))


def cutoff_boundary(K: float, U: float, J: float = 1.0) -> dict:
    """Kinematic support condition only; equality is a zero-measure edge.

    qmax=pi means no nonempty flat high-pass spectrum can close the leading
    dissociation channel. E>=4J closes the clean continuum channel, but does
    not imply P_bound(t)=1 after switching on disorder (virtual dressing).
    """
    E = bound_energy(K, U, J)
    if E >= 4 * J:
        return {"energy_closed": True, "Q": None, "qmin": None, "qmax": None}
    Q = float(2 * np.arccos(E / (4 * J)))
    return {"energy_closed": False, "Q": Q,
            "qmin": max(0.0, abs(K) - Q), "qmax": min(np.pi, abs(K) + Q)}


def elastic_backscatter_transfer(K: float) -> float:
    return float(abs((2 * K + np.pi) % (2 * np.pi) - np.pi))


def dissociation_kernel(q: float, K: float, U: float, J: float = 1.0) -> float:
    """Infinite relative-chain continuum spectral kernel, -2 Im <s|G|s>.

    The source s_r = 2 C a**abs(r) cos(q r / 2). The interacting continuum
    wave, not a free plane wave, is used. Bound-to-bound poles are excluded
    by restricting to the continuum |E|<2|t|. No artificial linewidth eta.
    See notes/THEORY.md for normalization and independent resolvent route.
    """
    E = bound_energy(K, U, J)
    t = 2 * J * np.cos((K + q) / 2)
    if E >= 2 * abs(t):
        return 0.0
    D = float(np.sqrt(4 * t * t - E * E))
    a = -4 * J * np.cos(K / 2) / (E + U)
    C = np.sqrt((1 - a * a) / (1 + a * a))
    p = np.arccos(np.clip(-E / (2 * t), -1.0, 1.0))
    z1, z2 = a * np.exp(1j * (p + q / 2)), a * np.exp(1j * (p - q / 2))
    R = z1 / (1 - z1) + z2 / (1 - z2)
    # Use a normalized scattering phase to avoid A=U/D overflow at edges.
    den = np.hypot(D, U)
    amp = 2 * C * ((1 + R.real) * D / den + R.imag * U * np.sign(t) / den)
    return float(2 * amp * amp / D)


def born_rate(K: float, U: float, q0: float, J: float = 1.0,
              fixed_variance: bool = True, epsabs: float = 1e-10,
              epsrel: float = 1e-9) -> tuple[float, float]:
    """Return (Gamma/W^2, numerical quadrature error estimate).

    S(q)=1/(1-q0/pi) on q0<|q|<=pi for fixed_variance=True; otherwise
    S(q)=1 there (fixed pass-band height). Error is not a confidence interval
    and excludes physical approximation error.
    """
    if not np.isfinite(q0) or not 0 <= q0 < np.pi:
        raise ValueError("Require 0 <= q0 < pi; pi leaves no disorder modes.")
    boundary = cutoff_boundary(K, U, J)
    if boundary["energy_closed"] or q0 >= boundary["qmax"]:
        return 0.0, 0.0
    Q = boundary["Q"]
    wrap = lambda x: (x + np.pi) % (2 * np.pi) - np.pi
    points = sorted(set(abs(wrap(x - K)) for x in (Q, -Q)))
    points = [p for p in points if q0 < p < np.pi]
    value, err = quad(lambda q: dissociation_kernel(q, K, U, J)
                     + dissociation_kernel(-q, K, U, J),
                     q0, np.pi, points=points, epsabs=epsabs,
                     epsrel=epsrel, limit=200)
    divisor = 2 * np.pi * ((1 - q0 / np.pi) if fixed_variance else 1)
    return float(value / divisor), float(err / divisor)


def relative_block(n: int, K: float, U: float, J: float = 1.0) -> np.ndarray:
    """Finite periodic block in |K,r>=sum_x exp(iKx)|x,x+r>/sqrt(n).

    K must be a ring momentum for comparison with the full n*n problem.
    This matrix routine itself also accepts other K for numerical checks.
    """
    if n < 3:
        raise ValueError("Require n>=3.")
    h = np.zeros((n, n), dtype=complex)
    r = np.arange(n)
    h[(r + 1) % n, r] += -J * (1 + np.exp(1j * K))
    h[(r - 1) % n, r] += -J * (1 + np.exp(-1j * K))
    h[0, 0] = U
    return h


def finite_born_coefficients(n: int, K: float, U: float,
                             fractions=(0.0, 0.8), etas=(0.08, 0.04),
                             J: float = 1.0) -> dict:
    """Independent finite-block eigensum using Gaussian-broadened delta.

    Finite n and eta are numerical regulators, not physical dissipation.
    Small eta at fixed n can be less accurate. No artificial delta at zero.
    """
    _check(K, U, J)
    if abs(K * n / (2 * np.pi) - round(K * n / (2 * np.pi))) > 1e-9:
        raise ValueError("K must be a permitted ring momentum.")
    if any(e <= 0 for e in etas):
        raise ValueError("eta>0 required.")
    if any(not 0 <= f < 1 for f in fractions):
        raise ValueError("cutoff fractions must lie in [0,1).")
    e0, v0 = eigh(relative_block(n, K, U, J))
    Eb, b = float(e0[-1]), v0[:, -1]
    qs, r = 2 * np.pi * np.fft.fftfreq(n), np.arange(n)
    kernels = np.zeros((len(etas), n))
    for iq, q in enumerate(qs):
        e, v = eigh(relative_block(n, K + q, U, J))
        f = v[:, :-1].conj().T @ ((1 + np.exp(1j * q * r)) * b)
        for ie, eta in enumerate(etas):
            delta = np.exp(-0.5 * ((Eb - e[:-1]) / eta) ** 2) / (np.sqrt(2*np.pi)*eta)
            kernels[ie, iq] = 2 * np.pi * np.dot(np.abs(f) ** 2, delta)
    rows = []
    for fraction in fractions:
        mask = np.abs(qs) > fraction * np.pi + 1e-12
        if not mask.any():
            raise ValueError("No allowed finite-ring Fourier modes.")
        S = mask.astype(float) / mask.mean()
        rows.append({"q0_over_pi": float(fraction), "allowed_modes": int(mask.sum()),
                     "coefficients": (kernels @ S / n).tolist()})
    return {"n": n, "K_over_pi": K/np.pi, "U": U, "J": J, "etas": list(etas),
            "bound_energy": Eb, "bound_energy_exact_infinite": bound_energy(K,U,J),
            "rows": rows}


def clean_bound_basis(n: int, U: float, J: float = 1.0) -> np.ndarray:
    return np.array([eigh(relative_block(n, 2*np.pi*m/n, U, J),
                          subset_by_index=(n-1,n-1))[1][:,0] for m in range(n)])


def initial_bound_state(basis: np.ndarray, momentum_index: int) -> np.ndarray:
    n = basis.shape[0]
    x, y = np.indices((n, n))
    m = momentum_index % n
    return (np.exp(2j*np.pi*m*x/n) * basis[m, (y-x)%n] / np.sqrt(n)).reshape(-1)


def bound_population(psi: np.ndarray, basis: np.ndarray) -> float:
    n = basis.shape[0]
    if psi.size != n*n:
        raise ValueError("State size must be n*n.")
    p = np.asarray(psi).reshape(n,n)
    x, r = np.indices((n,n))
    amplitudes = np.fft.fft(p[x,(x+r)%n], axis=0, norm='ortho')
    overlap = np.sum(basis.conj()*amplitudes,axis=1)
    return float(np.sum(np.abs(overlap)**2))


def two_body_hamiltonian(v: np.ndarray, U: float, J: float = 1.0) -> csr_matrix:
    v = np.asarray(v, dtype=float)
    n = len(v)
    if n < 3 or not np.isfinite(v).all():
        raise ValueError("Finite real potential with n>=3 required.")
    h = diags([-J*np.ones(n-1), -J*np.ones(n-1)],[-1,1],format='lil')
    h[0,n-1] = h[n-1,0] = -J
    h = h.tocsr()
    V = v[:,None] + v[None,:] + U*np.eye(n)
    return (kron(h,eye(n),format='csr') + kron(eye(n),h,format='csr')
            + diags(V.reshape(-1),format='csr')).tocsr()


def filtered_noise(n: int, q0: float, W: float, seed: int) -> np.ndarray:
    """Gaussian ensemble at fixed expected on-site variance W^2.

    Do not renormalize each sample: that would change the Gaussian ensemble.
    Using the same seed for white/high-pass gives a paired comparison.
    The mean Fourier mode is removed in both cases; q0=0 is white up to DC.
    """
    if n < 3 or not 0 <= q0 < np.pi or W < 0 or not np.isfinite(W):
        raise ValueError("Invalid n, q0 or W.")
    qs = 2*np.pi*np.fft.fftfreq(n)
    mask = (np.abs(qs) > q0+1e-12)
    if not mask.any():
        raise ValueError("No allowed modes.")
    raw = np.random.default_rng(seed).standard_normal(n)
    return W*np.fft.ifft(np.fft.fft(raw)*mask).real/np.sqrt(mask.mean())
