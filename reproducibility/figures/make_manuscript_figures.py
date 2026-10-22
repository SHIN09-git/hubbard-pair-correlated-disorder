"""Vector versions of the two manuscript figures, sized for a single A4 text width.

Reads only the saved figure data of the v1.5 deliverable (no recomputation of physics);
the plotted content is identical to docs/figures/01_support_and_normalization and 03_endpoint_laws,
with larger fonts. Run from the manuscript directory: python figures/make_manuscript_figures.py
"""
from pathlib import Path
import csv, hashlib, json, sys, time, platform
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"  # Byte-identical copies of the frozen v1.5 figure data.
COL = ["#2166ac", "#b35806", "#5e3c99"]
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"], "mathtext.fontset": "stix",
                     "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
                     "legend.fontsize": 8.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})


def rows(name):
    with open(DATA / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fig1():
    g = np.load(DATA / "01_phase_grid.npz")
    b = rows("01_boundaries.csv")
    r = rows("01_rates.csv")
    fig, ax = plt.subplots(1, 3, figsize=(6.6, 2.65), layout="constrained")
    ax[0].pcolormesh(g["K_over_pi"], g["U_over_J"], g["phase"], cmap=ListedColormap(["#e1e1e1", "#efbd95", "#8dbadd"]),
                     vmin=-.5, vmax=2.5, rasterized=True)
    k = np.array([float(x["K_over_pi"]) for x in b])
    ax[0].plot(k, [float(x["U_energy"]) for x in b], "k-", lw=.9)
    us = np.array([float(x["U_star"]) if x["U_star"] else np.nan for x in b])
    ax[0].plot(k, us, "k--", lw=.9)
    ax[0].text(.32, 4.75, "Energy closed", ha="center", fontsize=8.5)
    ax[0].text(.34, .5, "Finite gap\ncloses", ha="center", fontsize=8.5)
    ax[0].text(.82, 1.0, "No high-pass\nclosure", ha="center", fontsize=8.5)
    for u in (2, 3.5, 5):
        ax[0].plot(.75, u, "ko", ms=2.5)
    ax[0].set(xlabel=r"$|K|/\pi$", ylabel=r"$U/J$", title="(a)", xlim=(0, 1), ylim=(0, 5.5))
    for U, c in ((2.0, COL[0]), (3.5, COL[1])):
        sel = [x for x in r if float(x["U"]) == U]
        f = np.array([float(x["f"]) for x in sel])
        ax[1].plot(f, [float(x["c_variance"]) for x in sel], color=c, lw=1.2, label=rf"$U/J={U:g}$")
        ax[2].plot(f, [float(x["c_passband"]) for x in sel], color=c, lw=1.2, label=rf"$U/J={U:g}$")
    for a, t in ((ax[1], "(b) fixed variance"), (ax[2], "(c) fixed height")):
        a.axvline(.5, color="gray", ls=":", lw=.8, label=r"$q_B=0.5\pi$")
        a.axvline(.9416570958, color=COL[1], ls="--", lw=.8)
        a.set(xlabel=r"$q_0/\pi$", ylabel=r"$c=\Gamma^{(2)}/W^2$", title=t, xlim=(0, 1), ylim=(0, None))
        a.legend(frameon=False, loc="best")
    fig.savefig(HERE / "fig1_classification.pdf", dpi=600)
    plt.close(fig)


def fig2():
    r = rows("03_endpoints.csv")
    fig, ax = plt.subplots(1, 2, figsize=(6.0, 2.65), layout="constrained")
    lv = {"wrap_open": r"$U/J=2$: plateau", "critical": r"$U=U_*$: $\Delta^{1/2}$", "internal": r"$U/J=3.5$: $\Delta^{3/2}$"}
    lh = {"wrap_open": r"$U/J=2$: $\Delta$", "critical": r"$U=U_*$: $\Delta^{3/2}$", "internal": r"$U/J=3.5$: $\Delta^{3/2}$"}
    for i, mode in enumerate(("wrap_open", "critical", "internal")):
        sel = [x for x in r if x["mode"] == mode]
        d = np.array([float(x["delta"]) for x in sel])
        ax[0].loglog(d, [float(x["c"]) for x in sel], "o-", color=COL[i], ms=2.5, lw=1, label=lv[mode])
        ax[1].loglog(d, [float(x["passband"]) for x in sel], "o-", color=COL[i], ms=2.5, lw=1, label=lh[mode])
    for a, t in zip(ax, ("(a) fixed variance", "(b) fixed height")):
        a.set(xlabel=r"$\Delta=q_c-q_0$", ylabel=r"$c$", title=t)
        a.legend(frameon=False)
        a.grid(alpha=.15, which="major")
    fig.savefig(HERE / "fig2_endpoints.pdf", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    started = datetime.now(timezone.utc).isoformat()
    t0, c0 = time.perf_counter(), time.process_time()
    fig1(); fig2()
    inputs = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DATA.glob("0[13]_*"))}
    outputs = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob("*.pdf"))}
    (HERE / "figures_manifest.json").write_text(json.dumps(dict(
        source=dict(file="figures/make_manuscript_figures.py", sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        inputs=inputs, outputs=outputs, python=sys.version, matplotlib=matplotlib.__version__,
        started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
        executable=sys.executable, platform=platform.platform(), argv=sys.argv, seed=None,
        role="figure_formatting_only_saved_data_unchanged",
        wall_seconds=time.perf_counter() - t0, process_cpu_seconds=time.process_time() - c0, exit_code=0), indent=2) + "\n", encoding="utf-8")
    print("ok", list(outputs))
