"""Exploratory infinite-chain off-shell O(W^2) sinc^2 response; no dynamics."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, csv, hashlib, json, os, platform, sys, time, traceback
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import roots_legendre

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def save(p, value):
    p.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

def gauss(n, lo, hi):
    x, w = roots_legendre(n)
    return (lo + hi) / 2 + x * (hi - lo) / 2, w * (hi - lo) / 2

def cuts(K, E, J, fraction):
    q0 = np.pi * fraction
    points = [-np.pi, -q0, 0., q0, np.pi]
    # Split the signed relative hopping and the on-shell edges for quadrature.
    for n in range(-2, 3):
        zero = np.pi - K + 2 * np.pi * n
        if -np.pi < zero < np.pi:
            points.append(zero)
    if E < 4 * J:
        Q = 2 * np.arccos(E / (4 * J))
        for edge in [-K - Q, -K + Q]:
            points.append((edge + np.pi) % (2 * np.pi) - np.pi)
    return sorted(set(float(x) for x in points))

def matrix_element(q, p, K, U, J):
    E = np.sqrt(U * U + 16 * J * J * np.cos(K / 2)**2)
    a = -4 * J * np.cos(K / 2) / (E + U)
    C = np.sqrt(U / E)
    t = 2 * J * np.cos((K + q) / 2)
    d = 2 * np.abs(t) * np.sin(p)
    z1 = a * np.exp(1j * (p + q / 2))
    z2 = a * np.exp(1j * (p - q / 2))
    R = z1 / (1 - z1) + z2 / (1 - z2)
    amp = 2 * C / np.sqrt(np.pi) * (d * (1 + R.real) + U * np.sign(t) * R.imag) / np.hypot(d, U)
    return amp, -2 * t * np.cos(p)

def on_shell(q, K, U, J):
    E2 = U * U + 16 * J * J * np.cos(K / 2)**2
    t2 = 4 * J * J * np.cos((K + q) / 2)**2
    d2 = 4 * t2 - E2
    if d2 <= 0:
        return 0.
    B = 1 - 4 * np.cos(K / 2)**2 * t2 / E2
    return 8 * (U / np.sqrt(E2))**3 * np.sqrt(d2) / ((d2 + U * U) * B * B)

def compute(config):
    J, U = config["J"], config["U"]
    K = np.pi * config["K_over_pi"]
    E = np.sqrt(U * U + 16 * J * J * np.cos(K / 2)**2)
    a = -4 * J * np.cos(K / 2) / (E + U)
    C = np.sqrt(U / E)
    ts = np.arange(config["window"][0], config["window"][1] + config["time_step"] / 2, config["time_step"])
    centered = ts - ts.mean()
    fit_weights = centered / (centered @ centered)
    families = [(0., "ungapped", 1.), (config["gap_fraction"], "fixed_height", 1.),
                (config["gap_fraction"], "fixed_variance", 1 / (1 - config["gap_fraction"]))]
    coefficients = {}
    for fraction, label, scale in families:
        pts = cuts(K, E, J, fraction)
        value, error = 0., 0.
        for lo, hi in zip(pts[:-1], pts[1:]):
            if abs((lo + hi) / 2) < fraction * np.pi:
                continue
            v, e = quad(on_shell, lo, hi, args=(K, U, J), epsabs=1e-12, epsrel=1e-11, limit=200)
            value += v
            error += e
        coefficients[label] = {"c": value * scale / (2 * np.pi), "quad_error_estimate": error * scale / (2 * np.pi)}
    runs = []
    closure_max = 0.
    for order in config["quadrature_orders"]:
        p, wp = gauss(order, 0., np.pi)
        curves = {}
        for fraction in [0., config["gap_fraction"]]:
            loss = np.zeros(len(ts))
            pts = cuts(K, E, J, fraction)
            for lo, hi in zip(pts[:-1], pts[1:]):
                if abs((lo + hi) / 2) < fraction * np.pi:
                    continue
                q, wq = gauss(order, lo, hi)
                amp, energies = matrix_element(q[:, None], p[None, :], K, U, J)
                weights = amp**2 * wq[:, None] * wp[None, :] / (2 * np.pi)
                de = E - energies
                for it, T in enumerate(ts):
                    loss[it] += np.sum(weights * T**2 * np.sinc(de * T / (2 * np.pi))**2)
                if order == config["quadrature_orders"][-1]:
                    # Independent completeness check including the outgoing bound state.
                    t = 2 * J * np.cos((K + q) / 2)
                    Eb = np.sqrt(U * U + 4 * t * t)
                    ab = -2 * t / (Eb + U)
                    Cb = np.sqrt(U / Eb)
                    z = a * ab * np.exp(1j * q / 2)
                    bound_weight = (2 * C * Cb * (1 + 2 * (z / (1 - z)).real))**2
                    z = a * a * np.exp(1j * q)
                    norm_source = 4 * C * C * (1 + a * a / (1 - a * a) + (z / (1 - z)).real)
                    cont_weight = amp**2 @ wp
                    closure_max = max(closure_max, float(np.max(abs(cont_weight + bound_weight - norm_source))))
            curves[fraction] = loss
        rows = []
        for fraction, label, scale in families:
            loss = curves[fraction] * scale
            slope = float(loss @ fit_weights)
            intercept = float(loss.mean() - slope * ts.mean())
            c = coefficients[label]["c"]
            rows.append({"family": label, "gap_fraction": fraction, "slope": slope, "intercept": intercept,
                         "c": c, "slope_over_c": slope / c, "relative_difference": slope / c - 1,
                         "loss_over_W2": loss.tolist()})
        runs.append({"q_order_per_segment": order, "p_order": order, "rows": rows})
    # Real-space inner product checks avoid reusing the geometric-series formula.
    radius = 80
    r = np.arange(-radius, radius + 1)
    direct_error = 0.
    for q in [-2.9, -1.2, 0., 0.9, 2.8]:
        for p in [0.17, 1.05, 2.65]:
            t = 2 * J * np.cos((K + q) / 2)
            d = 2 * abs(t) * np.sin(p)
            phi = (d * np.cos(p * r) + U * np.sign(t) * np.sin(p * abs(r))) / np.sqrt(np.pi * (d * d + U * U))
            source = 2 * C * a**abs(r) * np.cos(q * r / 2)
            expected = float(phi @ source)
            actual = float(matrix_element(q, p, K, U, J)[0])
            direct_error = max(direct_error, abs(expected - actual))
    qzero = float(np.max(abs(matrix_element(0., np.linspace(.01, np.pi - .01, 201), K, U, J)[0])))
    difference = max(abs(a["slope"] - b["slope"]) for a, b in zip(runs[-1]["rows"], runs[-2]["rows"]))
    checks = {"last_two_grid_max_slope_difference": difference, "completeness_max_absolute_error": closure_max,
              "real_space_matrix_element_max_error": direct_error, "uniform_q_zero_max_amplitude": qzero,
              "real_space_radius": radius}
    assert difference < config["convergence_absolute_tolerance"], checks
    assert max(closure_max, direct_error, qzero) < 1e-10, checks
    return {"role": config["role"], "times": ts.tolist(), "quadrature_runs": runs, "rows": runs[-1]["rows"],
            "coefficient_integrals": coefficients, "checks": checks,
            "scope": "Infinite-chain ensemble-averaged O(W^2) response, J=1 U=2 K=0.75pi, original T=12:2:24 OLS with intercept. Does not exclude finite-N/time coupling, finite-W effects, other windows or endpoint regimes."}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--config", default="configs/infinite_chain_finite_time_v1.json")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = args.root.resolve()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    start_utc = datetime.now(timezone.utc).isoformat()
    wall, cpu = time.perf_counter(), time.process_time()
    config_path = root / args.config
    config = json.loads(config_path.read_text("utf-8"))
    source_time_file = root / config["reference_time_file"]
    old_times = json.loads(source_time_file.read_text("utf-8"))["times"]
    inputs = {args.config: sha(config_path), config["reference_time_file"]: sha(source_time_file)}
    manifest = {"started_utc": start_utc, "role": config["role"], "argv": sys.argv, "cwd": str(Path.cwd()),
                "config": config, "seed": None, "input_hashes": inputs,
                "source": {"file": "tools/infinite_chain_finite_time.py", "sha256": sha(Path(__file__))},
                "environment": {"python": sys.version, "executable": sys.executable, "numpy": np.__version__,
                                "scipy": scipy.__version__, "platform": platform.platform(), "cpu_count": os.cpu_count()},
                "status": "RUNNING"}
    save(out / "manifest.json", manifest)
    exit_code = 0
    try:
        result = compute(config)
        assert result["times"] == [t for t in old_times if config["window"][0] <= t <= config["window"][1]]
        result["checks"]["original_time_grid_identical"] = True
        save(out / "summary.json", result)
        with (out / "curves.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["family", "T", "L2_over_W2"])
            for row in result["rows"]:
                writer.writerows(zip([row["family"]] * len(result["times"]), result["times"], row["loss_over_W2"]))
        manifest["status"] = "PASS"
        print(json.dumps({"rows": [{k: v for k, v in row.items() if k != "loss_over_W2"} for row in result["rows"]],
                          "checks": result["checks"]}, ensure_ascii=False, indent=2))
    except Exception:
        exit_code = 1
        manifest["status"] = "FAIL"
        manifest["traceback"] = traceback.format_exc()
        print(manifest["traceback"], file=sys.stderr)
    finally:
        manifest.update(finished_utc=datetime.now(timezone.utc).isoformat(), wall_seconds=time.perf_counter() - wall,
                        process_cpu_seconds=time.process_time() - cpu, exit_code=exit_code)
        save(out / "manifest.json", manifest)
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
