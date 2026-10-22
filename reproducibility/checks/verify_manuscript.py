"""Independent checks of the analytic statements and numbers quoted in main.tex.

Exploratory verification for the English draft; no dynamics, no saved data are modified.
Run from the manuscript directory:  python checks/verify_manuscript.py --output checks/results/<new_dir>
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, os, platform, sys, time, traceback
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import solve_banded

J = 1.0
SEED = 20261014
FROZEN_C = {"ungapped": 1.0611341268352634, "fixed_variance": 2.0133460220473625, "fixed_height": 0.4026692044094724}


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def bound(K, U):
    t0 = 2 * J * np.cos(K / 2)
    E = np.sqrt(U * U + 4 * t0 * t0)
    return t0, E, -2 * t0 / (E + U), np.sqrt(U / E)


def kernel(K, U, q):
    """Equation (5)."""
    t0, E, a, C = bound(K, U)
    t = 2 * J * np.cos((K + q) / 2)
    D2 = 4 * t * t - E * E
    B = 1 - 4 * np.cos(K / 2) ** 2 * t * t / E ** 2
    D = np.sqrt(np.clip(D2, 0, None))
    return np.where(D2 > 0, 8 * (U / E) ** 3 * D / ((D2 + U * U) * B * B), 0.0)


def kernel_resolvent(K, U, q, L, eta):
    """-2 Im <s|(E+i eta-h)^{-1}|s> on the relative chain r=-L..L (independent of eq. 5)."""
    t0, E, a, C = bound(K, U)
    t = 2 * J * np.cos((K + q) / 2)
    r = np.arange(-L, L + 1)
    s = 2 * C * a ** np.abs(r) * np.cos(q * r / 2)
    ab = np.zeros((3, len(r)), complex)
    ab[0, 1:] = t
    ab[2, :-1] = t
    ab[1, :] = E + 1j * eta
    ab[1, L] -= U
    x = solve_banded((1, 1), ab, s.astype(complex))
    return -2 * np.imag(np.vdot(s, x))


def kinematics(K, U):
    k = abs(K)
    _, E, _, _ = bound(K, U)
    if E >= 4 * J:
        return dict(E=E, Q=None, qc=None)
    Q = 2 * np.arccos(E / (4 * J))
    return dict(E=E, Q=Q, qc=min(np.pi, k + Q))


def coefficient(K, U, q0, family):
    _, E, _, _ = bound(K, U)
    pts = [-np.pi, -q0, q0, np.pi]
    kin = kinematics(K, U)
    if kin["Q"] is not None:
        pts += [wrap(-K - kin["Q"]), wrap(-K + kin["Q"])]
    pts = sorted(set(p for p in pts if -np.pi <= p <= np.pi))
    total = 0.0
    for lo, hi in zip(pts[:-1], pts[1:]):
        if abs((lo + hi) / 2) <= q0:
            continue
        total += quad(lambda q: float(kernel(K, U, q)), lo, hi, epsabs=1e-13, epsrel=1e-12, limit=400)[0]
    c = total / (2 * np.pi)
    return c / (1 - q0 / np.pi) if family == "V" else c


def run():
    out = {}
    rng = np.random.default_rng(SEED)

    # 1. Algebraic identities of appendix A.3 at random points.
    err = dict(one_plus_ReR=0.0, ImR=0.0, identity1=0.0, identity2_rel=0.0, F_rel=0.0, overlap=0.0)
    for _ in range(4000):
        K = rng.uniform(-np.pi, np.pi); U = rng.uniform(0.05, 5); q = rng.uniform(-np.pi, np.pi); p = rng.uniform(0.01, np.pi - 0.01)
        t0, E, a, C = bound(K, U); h = q / 2; t = 2 * J * np.cos(K / 2 + h)
        R = sum(a * np.exp(1j * (p + s * h)) / (1 - a * np.exp(1j * (p + s * h))) for s in (1, -1))
        dp = 1 - 2 * a * np.cos(p + h) + a * a; dm = 1 - 2 * a * np.cos(p - h) + a * a
        err["one_plus_ReR"] = max(err["one_plus_ReR"], abs(1 + R.real - (1 - a * a) * (1 + a * a - 2 * a * np.cos(p) * np.cos(h)) / (dp * dm)))
        err["ImR"] = max(err["ImR"], abs(R.imag - 2 * a * np.sin(p) * ((1 + a * a) * np.cos(h) - 2 * a * np.cos(p)) / (dp * dm)))
        err["identity1"] = max(err["identity1"], abs(t * t + t0 * t0 - 2 * t * t0 * np.cos(h) - 4 * J * J * np.sin(h) ** 2))
        # off-shell overlap <phi|s> = 2 C D_p F / sqrt(pi (D_p^2+U^2)) against a real-space sum
        r = np.arange(-120, 121); Dp = 2 * abs(t) * np.sin(p)
        phi = (Dp * np.cos(p * r) + U * np.sign(t) * np.sin(p * np.abs(r))) / np.sqrt(np.pi * (Dp * Dp + U * U))
        F = 1 + R.real + U / (2 * t * np.sin(p)) * R.imag
        if abs(a) < 0.85:
            err["overlap"] = max(err["overlap"], abs(phi @ (2 * C * a ** np.abs(r) * np.cos(q * r / 2)) - 2 * C * Dp * F / np.sqrt(np.pi * (Dp * Dp + U * U))))
        if 4 * t * t > E * E:
            ps = np.arccos(-E / (2 * t))
            dp = 1 - 2 * a * np.cos(ps + h) + a * a; dm = 1 - 2 * a * np.cos(ps - h) + a * a
            lhs = t * t * dp * dm; rhs = 4 * np.sin(h) ** 2 * (J * J * (1 + a * a) ** 2 - a * a * t * t)
            err["identity2_rel"] = max(err["identity2_rel"], abs(lhs - rhs) / (1 + abs(lhs)))
            R = sum(a * np.exp(1j * (ps + s * h)) / (1 - a * np.exp(1j * (ps + s * h))) for s in (1, -1))
            F = 1 + R.real + U / (2 * t * np.sin(ps)) * R.imag
            B = 1 - 4 * np.cos(K / 2) ** 2 * t * t / E ** 2
            err["F_rel"] = max(err["F_rel"], abs(F - (U / E) / B) / abs(F))
    out["identities_max_error"] = err

    # 2. Kernel versus resolvent at 25 random open points (eta -> 0 by Richardson).
    rows = []
    while len(rows) < 25:
        K = rng.uniform(-np.pi, np.pi); U = rng.uniform(0.1, 5); q = rng.uniform(-np.pi, np.pi)
        f = float(kernel(K, U, q))
        if f == 0:
            continue
        b = 2 * kernel_resolvent(K, U, q, 300000, 1e-3) - kernel_resolvent(K, U, q, 300000, 2e-3)
        rows.append(dict(K=K, U=U, q=q, formula=f, resolvent=b, rel=(b - f) / f))
    out["resolvent_check"] = dict(points=rows, max_abs_rel=max(abs(r["rel"]) for r in rows),
                                  chain_sites=600001, etas=[1e-3, 2e-3], extrapolation="2*K(1e-3)-K(2e-3)")

    # 3. Numbers at K = 0.75 pi.
    K = 0.75 * np.pi; k = K
    nums = {}
    for U in (2.0, 3.5, 5.0):
        kin = kinematics(K, U)
        nums[f"U={U}"] = dict(E=kin["E"], Q_over_pi=None if kin["Q"] is None else kin["Q"] / np.pi,
                             qc_over_pi=None if kin["qc"] is None else kin["qc"] / np.pi)
    nums["Ustar"] = 4 * J * np.sqrt(-np.cos(k)); nums["U_E"] = 4 * J * np.sin(k / 2)
    nums["qB_over_pi"] = abs(wrap(2 * K)) / np.pi
    c0 = coefficient(K, 2.0, 0.0, "H"); cV = coefficient(K, 2.0, 0.8 * np.pi, "V"); cH = coefficient(K, 2.0, 0.8 * np.pi, "H")
    nums["c_U2"] = dict(ungapped=c0, fixed_variance=cV, fixed_height=cH,
                        rel_to_frozen=dict(ungapped=c0 / FROZEN_C["ungapped"] - 1, fixed_variance=cV / FROZEN_C["fixed_variance"] - 1,
                                           fixed_height=cH / FROZEN_C["fixed_height"] - 1),
                        ratio_V=cV / c0, ratio_H=cH / c0)
    out["numbers_K075"] = nums

    # 4. Corollary 1 and window (8): random scan against direct kinematics.
    mism_star = mism_window = 0; n_star = n_window = 0
    for _ in range(20000):
        kk = rng.uniform(np.pi / 2 + 1e-3, np.pi - 1e-3); U = rng.uniform(0.01, 4.0)
        kin = kinematics(kk, U)
        if kin["Q"] is None:
            continue
        n_star += 1
        if (kin["qc"] >= np.pi - 1e-12) != (U <= 4 * J * np.sqrt(-np.cos(kk))):
            mism_star += 1
        if kk < 2 * np.pi / 3:
            n_window += 1
            qB = abs(wrap(2 * kk))
            if (qB < kin["qc"]) != (U < 4 * J * np.sqrt(-np.sin(2 * kk) * np.sin(kk))):
                mism_window += 1
    out["classification_scan"] = dict(points_Ustar=n_star, mismatches_Ustar=mism_star, points_window=n_window, mismatches_window=mism_window)

    # 5. Backscattering overlap and group velocity.
    Kb, Ub = 0.6 * np.pi, 1.3
    t0, E, a, C = bound(Kb, Ub); r = np.arange(-400, 401)
    FB_sum = 2 * C * C * np.sum(a ** (2 * np.abs(r)) * np.cos(Kb * r)); FB = 2 * (1 - a * a) ** 2 / (1 - 2 * a * a * np.cos(Kb) + a ** 4)
    Eb = lambda x: np.sqrt(Ub * Ub + 16 * J * J * np.cos(x / 2) ** 2); d = 1e-6
    out["backscattering"] = dict(FB_sum=FB_sum, FB_formula=FB, vB_numeric=abs((Eb(Kb + d) - Eb(Kb - d)) / (2 * d)), vB_formula=4 * J * J * abs(np.sin(Kb)) / E)

    # 6. Endpoint asymptotics (proposition 3) at K = 0.75 pi.
    ends = {}
    Ustar = nums["Ustar"]
    for label, U in (("interior_U3.5", 3.5), ("critical_Ustar", Ustar), ("wrapped_U2", 2.0)):
        kin = kinematics(K, U); Q = kin["Q"]; qc = kin["qc"]; _, E, _, _ = bound(K, U)
        A = 8 * U * np.sqrt(8 * J * J * np.sin(Q)) / (E ** 3 * np.sin(k / 2) ** 4)
        res = []
        for Delta in (1e-3, 1e-4):
            q0 = qc - Delta
            cH = coefficient(K, U, q0, "H"); cV = coefficient(K, U, q0, "V")
            if label.startswith("interior"):
                pH = A * Delta ** 1.5 / (3 * np.pi); pV = pH / (1 - qc / np.pi)
            elif label.startswith("critical"):
                pH = A * Delta ** 1.5 / (3 * np.pi); pV = A * np.sqrt(Delta) / 3
            else:
                Kpi = float(kernel(K, U, np.pi)); pH = Kpi * Delta / np.pi; pV = Kpi
            res.append(dict(Delta=Delta, cH=cH, cH_pred=pH, cH_ratio=cH / pH, cV=cV, cV_pred=pV, cV_ratio=cV / pV))
        ends[label] = dict(A=A, qc_over_pi=qc / np.pi, rows=res)
    out["endpoint_asymptotics"] = ends

    # 7. Small-U limit (6).
    smallU = []
    for U in (0.02, 0.005, 0.002):
        q = np.linspace(-np.pi, np.pi, 4000001)
        smallU.append(dict(U=U, c=float(np.trapezoid(kernel(K, U, q), q) / (2 * np.pi))))
    out["small_U"] = dict(rows=smallU, limit=1 / (J * np.sin(k / 2)))

    # 8. Lemma 1: first kinematically allowed order at U = 3.5.
    kin = kinematics(K, 3.5); Q, qc = kin["Q"], kin["qc"]
    nk = {}
    for f in (0.95, 0.97, 0.99):
        delta = np.pi * (1 - f)
        n = 1
        while not ((n % 2 == 1 and n * delta > np.pi - qc) or (n % 2 == 0 and n * delta > k - Q)):
            n += 1
        nk[str(f)] = n
    out["n_kin_U3.5"] = dict(n_kin=nk, Omega_over_pi=[(-k - Q) / np.pi, (-k + Q) / np.pi])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    outdir = args.output.resolve(); outdir.mkdir(parents=True, exist_ok=False)
    me = Path(__file__).resolve()
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(), role="exploratory_manuscript_verification",
                    argv=sys.argv, cwd=str(Path.cwd()), seed=SEED,
                    source=dict(file="checks/verify_manuscript.py", sha256=hashlib.sha256(me.read_bytes()).hexdigest()),
                    environment=dict(python=sys.version, executable=sys.executable, numpy=np.__version__, scipy=scipy.__version__,
                                     platform=platform.platform(), cpu_count=os.cpu_count()),
                    status="RUNNING")
    wall, cpu = time.perf_counter(), time.process_time(); code = 0
    try:
        res = run()
        (outdir / "summary.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")
        manifest["status"] = "PASS"
        print(json.dumps({k: v for k, v in res.items() if k != "resolvent_check"}, indent=1, default=float))
        print("resolvent max |rel| =", res["resolvent_check"]["max_abs_rel"])
    except Exception:
        code = 1; manifest["status"] = "FAIL"; manifest["traceback"] = traceback.format_exc(); print(manifest["traceback"], file=sys.stderr)
    finally:
        manifest.update(finished_utc=datetime.now(timezone.utc).isoformat(), wall_seconds=time.perf_counter() - wall,
                        process_cpu_seconds=time.process_time() - cpu, exit_code=code)
        (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    sys.exit(code)


if __name__ == "__main__":
    main()
