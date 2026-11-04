#!/usr/bin/env python3
"""Fig. 5: simulated time-to-peak (TTP) vs. empirical F-TRACT peak delay.

Reproduces the MATLAB figure made by ``yihang_ftract/RC_delay_scatter_plot_Amp.m``
(all pairs) and ``RC_delay_scatter_plot_Prob.m`` (unconnected pairs), whose
simulated TTP came from ``RC_delay_scatter_refine.m`` (SC_data_neuron, edges
with F < 0.25 removed, C = 1, V0(stim) = 1, 200 ms / 2,000 points, TTP >= 200
ms set to NaN), and recomputes the same statistics with the X_struct.v8.W run
behind Fig. 3 (``rc_peak_xstruct_w.npz``).

Masks (MATLAB): finite empirical peak delay, 0 < delay <= 200 ms,
N_responses = N_stimulations x probability >= 100, finite simulated TTP;
unconnected = F <= 0.25.  ``fig3`` masks add the Fig. 3 criteria (same
hemisphere, ED >= 20 mm, N_implantations >= 5).

Outputs (``--output-dir``): ``ttp_variants.csv`` (Spearman / partial Spearman
with 95 % bootstrap CIs for every variant), ``provenance.json`` and, for the
chosen ``--plot-variant`` / ``--plot-mask``, ``fig5_ttp.png`` / ``.pdf``.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402
from scipy.stats import pearsonr, rankdata, spearmanr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rc_ftract_core as rcc  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_DATA_ROOT = rcc.DEFAULT_DATA_ROOT
DEFAULT_RC_NPZ = PROJECT_ROOT / "output" / "fig3a_bootstrap_ci" / "rc_peak_xstruct_w.npz"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "fig5_ttp"
RC_RED = (0.8392, 0.1529, 0.1569)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def simulate_ttp(SC: np.ndarray, D: np.ndarray, order: np.ndarray, num_points: int = 2000) -> np.ndarray:
    """Peak time (ms) of the total outgoing |current|, [target, stimulus], cortical block in ``order``."""
    SC = np.where(np.isfinite(SC), SC, 0.0)
    D = np.where(np.isfinite(D), D, 0.0)
    G = np.zeros_like(SC)
    edges = (SC > 0) & (D > 0)
    G[edges] = SC[edges] / D[edges]
    G = 0.5 * (G + G.T)
    _, ttp, _ = rcc.simulate_rc(G, 200.0, num_points, G.shape[0])
    return ttp[np.ix_(order, order)][:360, :360]


def clean_ttp(ttp: np.ndarray) -> np.ndarray:
    out = ttp.astype(float).copy()
    out[out >= 200.0] = np.nan  # MATLAB: TTP(TTP >= 200) = NaN
    return out


def stats_with_ci(sim: np.ndarray, emp: np.ndarray, dist: np.ndarray, n_boot: int, rng: np.random.Generator):
    full = spearmanr(sim, emp).statistic
    part = rcc.partial_spearman(sim, emp, dist)
    n = sim.size
    boot_full = np.empty(n_boot)
    boot_part = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        s, e, d = sim[idx], emp[idx], dist[idx]
        boot_full[b] = pearsonr(rankdata(s), rankdata(e)).statistic
        boot_part[b] = rcc.partial_spearman(s, e, d)
    return {
        "n": n, "rho": full, "rho_lo": np.percentile(boot_full, 2.5), "rho_hi": np.percentile(boot_full, 97.5),
        "partial": part, "partial_lo": np.nanpercentile(boot_part, 2.5), "partial_hi": np.nanpercentile(boot_part, 97.5),
    }


def binned_median(x: np.ndarray, y: np.ndarray, n_bins: int = 10):
    edges = np.quantile(x, np.linspace(0, 1, n_bins + 1))
    centres, med, q25, q75 = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (x >= lo) & (x <= hi)
        if sel.sum() < 5:
            continue
        centres.append(np.median(x[sel]))
        med.append(np.median(y[sel]))
        q25.append(np.percentile(y[sel], 25))
        q75.append(np.percentile(y[sel], 75))
    return np.array(centres), np.array(med), np.array(q25), np.array(q75)


def binned_median_ci(x: np.ndarray, y: np.ndarray, n_bins: int = 10, n_boot: int = 1000, seed: int = 17):
    """Median of y within quantile bins of x with a 95 % bootstrap CI of that median (resampling pairs within the bin)."""
    rng = np.random.default_rng(seed)
    edges = np.quantile(x, np.linspace(0, 1, n_bins + 1))
    centres, med, lo_ci, hi_ci = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (x >= lo) & (x <= hi)
        if sel.sum() < 5:
            continue
        ys = y[sel]
        boots = np.median(ys[rng.integers(0, ys.size, size=(n_boot, ys.size))], axis=1)
        centres.append(np.median(x[sel]))
        med.append(np.median(ys))
        lo_ci.append(np.percentile(boots, 2.5))
        hi_ci.append(np.percentile(boots, 97.5))
    return np.array(centres), np.array(med), np.array(lo_ci), np.array(hi_ci)


DELAY_XLABEL = {"peak_delay": "Empirical delay-to-peak (ms)", "onset_delay": "Empirical onset latency (ms)",
                "latency_start": "Empirical latency start (ms)"}


def draw_panel(ax, emp, sim, result: dict, title: str, font: float, axis_max: float,
               xlabel: str = DELAY_XLABEL["peak_delay"], x_max: float | None = None, y_max: float | None = None,
               point_size: float = 9.0, point_alpha: float = 0.25, point_color: str = "0.35",
               point_style: str = "plain", band: str = "iqr", ylabel: str = "Simulated time-to-peak (ms)") -> None:
    # x_max / y_max override the common axis_max (e.g. onset latency never exceeds ~145 ms while TTP reaches 200 ms);
    # with unequal limits the panel keeps its square cell (aspect auto) and the identity line is still y = x.
    x_max = axis_max if x_max is None else x_max
    y_max = axis_max if y_max is None else y_max
    if point_style == "density":
        # colour each pair by the number of pairs in its 2-D bin (log scale); dense points drawn last
        nb = 40
        H, xe, ye = np.histogram2d(emp, sim, bins=[np.linspace(0, x_max, nb + 1), np.linspace(0, y_max, nb + 1)])
        ix = np.clip(np.digitize(emp, xe) - 1, 0, nb - 1)
        iy = np.clip(np.digitize(sim, ye) - 1, 0, nb - 1)
        dens = np.log1p(H[ix, iy])
        order = np.argsort(dens)
        cmap = LinearSegmentedColormap.from_list("blues_trunc", plt.get_cmap("Blues")(np.linspace(0.3, 1.0, 256)))
        ax.scatter(emp[order], sim[order], c=dens[order], cmap=cmap, vmin=0, s=point_size, alpha=0.9, linewidths=0,
                   zorder=2, rasterized=True)
    elif point_style == "hexbin":
        cmap = LinearSegmentedColormap.from_list("blues_trunc", plt.get_cmap("Blues")(np.linspace(0.2, 1.0, 256)))
        ax.hexbin(emp, sim, gridsize=(28, int(round(28 * y_max / x_max / 1.15))), extent=(0, x_max, 0, y_max), cmap=cmap,
                  bins="log", mincnt=1, linewidths=0.2, edgecolors="face", zorder=2)
    else:
        ax.scatter(emp, sim, s=point_size, color=point_color, alpha=point_alpha, linewidths=0, zorder=2, rasterized=True)
    n_bins = 10 if emp.size >= 1000 else 6
    if band == "ci95":
        # band = 95 % bootstrap CI of the binned median (uncertainty of the trend, not the spread of the pairs)
        cx, med, lo_b, hi_b = binned_median_ci(emp, sim, n_bins=n_bins)
        ax.fill_between(cx, lo_b, hi_b, color=RC_RED, alpha=0.35, linewidth=0, zorder=3)
        label = "Binned median (95% CI)"
    else:
        cx, med, lo_b, hi_b = binned_median(emp, sim, n_bins=n_bins)
        ax.fill_between(cx, lo_b, hi_b, color=RC_RED, alpha=0.18 if point_style == "plain" else 0.28, linewidth=0, zorder=3)
        label = "Binned median (IQR band)"
    ax.plot(cx, med, "-o", color=RC_RED, linewidth=2.2, markersize=5, zorder=4, label=label)
    lim = max(x_max, y_max)
    ax.plot([0, lim], [0, lim], "--", color="0.5", linewidth=1.2, zorder=1, label="y = x")
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_aspect("equal" if x_max == y_max else "auto")
    ax.set_xlabel(xlabel, fontsize=font + 1)
    ax.set_ylabel(ylabel, fontsize=font + 1)
    ax.set_title(title, fontsize=font + 1, fontweight="bold")
    ax.tick_params(labelsize=font - 1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.2)
    text = (f"r = {result['rho']:.2f} [{result['rho_lo']:.2f}, {result['rho_hi']:.2f}]\n"
            f"r | ED = {result['partial']:.2f} [{result['partial_lo']:.2f}, {result['partial_hi']:.2f}]")
    ax.text(0.03, 0.97, text, transform=ax.transAxes, ha="left", va="top", fontsize=font,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.8"))
    ax.legend(loc="lower right", frameon=False, fontsize=font - 3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--rc-npz", type=Path, default=DEFAULT_RC_NPZ)
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--skip-neuron", action="store_true", help="do not re-simulate the SC_data_neuron variant")
    parser.add_argument("--neuron-npz", type=Path, default=None,
                        help="saved raw-order TTP of the SC_data_neuron (F < 0.25 removed) run, skips its simulation")
    parser.add_argument("--plot-variant", default="xstruct", choices=("xstruct", "neuron_F025"))
    parser.add_argument("--plot-mask", default="fig3", choices=("matlab", "fig3", "matlab_diag", "no_ed"))
    parser.add_argument("--unconnected-def", default="F025", choices=("F025", "W0"), help="unconnected pairs: F025 = structural connection probability F <= 0.25 (connected F >= 0.75); W0 = no direct edge in the group connectome, W = 0 (connected W > 0)")
    parser.add_argument("--xstruct-file", type=Path, default=None,
                        help="X_struct .mat for the W = 0 definition (default <data-root>/yihang_ftract/X_struct_v1-19_FCO.mat)")
    parser.add_argument("--suffix", default="", help="extra suffix appended to the output file names")
    parser.add_argument("--count-variable", default="responses", choices=("responses", "stimulations"),
                        help="sampling criterion for a pair: N_responses >= 100 (default) or N_stimulations >= 100 (the Fig. 3 probability criterion)")
    parser.add_argument("--delay-variable", default="peak_delay", choices=tuple(DELAY_XLABEL),
                        help="empirical latency variable (F-TRACT <variable>__median.txt.gz)")
    parser.add_argument("--axis-max", type=float, default=150.0)
    parser.add_argument("--font-size", type=float, default=17.0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    data = rcc.load_comparison_data(args.data_root)
    order_sc = rcc.load_permutation(args.data_root / "ReOrder_Caio_SC.mat", "ReOrder_Caio_SC", 374)
    order_ft = rcc.load_permutation(args.data_root / "ReOreder_F_tract.mat", "ReOrder_F_tract", 360, column=1)
    delay_file = rcc.ftract_dir(args.data_root) / f"{args.delay_variable}__median.txt.gz"
    delay = rcc.read_ftract_gzip(delay_file)[np.ix_(order_ft, order_ft)]
    n_resp = data.n_stimulations * data.probability
    hemi = np.zeros((360, 360), bool)
    hemi[:180, :180] = True
    hemi[180:, 180:] = True
    count = n_resp if args.count_variable == "responses" else data.n_stimulations
    base_diag = np.isfinite(delay) & (delay > 0) & (delay <= 200) & (count >= 100)
    base = base_diag & ~np.eye(360, dtype=bool)
    fig3 = base & hemi & (data.distance >= 20) & (data.n_implantations >= 5)
    if args.unconnected_def == "W0":
        xs_file = args.xstruct_file or (rcc.xstruct_path(args.data_root))
        xs = rcc.load_xstruct_v8(xs_file)
        W = np.asarray(xs.W, dtype=float)[np.ix_(order_sc, order_sc)][:360, :360]
        unconnected = W == 0  # no direct edge in the group connectome (F < ~0.4 after consistency thresholding)
    else:
        unconnected = data.frequency <= 0.25
    no_ed = base & hemi & (data.n_implantations >= 5)  # Fig. 3 criteria without the ED >= 20 mm cut
    masks = {
        # as the MATLAB plot scripts: the diagonal (stimulation-site self-pairs, simulated TTP = 0) is included
        ("matlab_diag", "all"): base_diag,
        ("matlab_diag", "unconnected"): base_diag & unconnected,
        ("matlab", "all"): base,
        ("matlab", "unconnected"): base & unconnected,
        ("fig3", "all"): fig3,
        ("fig3", "unconnected"): fig3 & unconnected,
        ("no_ed", "all"): no_ed,
        ("no_ed", "unconnected"): no_ed & unconnected,
    }

    variants: dict[str, np.ndarray] = {}
    with np.load(args.rc_npz) as npz:
        variants["xstruct"] = clean_ttp(npz["ttp_ms"][np.ix_(order_sc, order_sc)][:360, :360])
    if args.neuron_npz is not None:
        with np.load(args.neuron_npz) as npz:
            variants["neuron_F025"] = clean_ttp(npz["ttp_ms"][np.ix_(order_sc, order_sc)][:360, :360])
    elif not args.skip_neuron:
        mat = loadmat(args.data_root / "SC_data_neuron_2023_yihang.mat", squeeze_me=True)
        SC, D, F = (np.asarray(mat[k], dtype=float) for k in ("SC", "D", "F"))
        SC_thr = SC.copy()
        SC_thr[F < 0.25] = 0.0  # RC_delay_scatter_refine.m line 169
        print("Simulating SC_data_neuron (F < 0.25 removed) ...", flush=True)
        variants["neuron_F025"] = clean_ttp(simulate_ttp(SC_thr, D, order_sc))

    rows = []
    results: dict[tuple[str, str, str], dict] = {}
    for variant, ttp in variants.items():
        for (mask_name, pair_set), mask in masks.items():
            valid = mask & np.isfinite(ttp)
            res = stats_with_ci(ttp[valid], delay[valid], data.distance[valid], args.n_boot, rng)
            results[(variant, mask_name, pair_set)] = res
            rows.append({"variant": variant, "mask": mask_name, "pair_set": pair_set, **res})
            print(f"{variant:12s} {mask_name:6s} {pair_set:11s} n = {res['n']:5d}  "
                  f"rho = {res['rho']:.3f} [{res['rho_lo']:.3f}, {res['rho_hi']:.3f}]  "
                  f"partial = {res['partial']:.3f} [{res['partial_lo']:.3f}, {res['partial_hi']:.3f}]", flush=True)
    suffix = ("" if args.count_variable == "responses" else "_nstim") + \
        ("" if args.delay_variable == "peak_delay" else f"_{args.delay_variable}") + \
        ("" if args.unconnected_def == "F025" else "_W0") + args.suffix
    with (args.output_dir / f"ttp_variants{suffix}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # figure for the chosen variant
    ttp = variants[args.plot_variant]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.2))
    for ax, pair_set, letter in zip(axes, ("all", "unconnected"), ("A", "B")):
        valid = masks[(args.plot_mask, pair_set)] & np.isfinite(ttp)
        res = results[(args.plot_variant, args.plot_mask, pair_set)]
        label = "All pairs" if pair_set == "all" else "Unconnected pairs"
        draw_panel(ax, delay[valid], ttp[valid], res, f"{label} (n = {res['n']:,})", args.font_size, args.axis_max,
                   xlabel=DELAY_XLABEL[args.delay_variable])
        ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontsize=args.font_size + 9, fontweight="bold",
                ha="left", va="bottom")
    fig.subplots_adjust(left=0.08, right=0.99, top=0.92, bottom=0.11, wspace=0.28)
    png = args.output_dir / f"fig5_ttp{suffix}.png"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(png.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    provenance = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": "code/python_figures/fig5_ttp_delay.py",
        "empirical_delay": str(delay_file),
        "count_variable": args.count_variable,
        "unconnected_definition": args.unconnected_def,
        "empirical_delay_sha256": sha256(delay_file),
        "rc_npz": str(args.rc_npz),
        "rc_npz_sha256": sha256(args.rc_npz),
        "ttp_definition": "time of the maximum of the total outgoing |current| (2,000 points over 200 ms); TTP >= 200 ms set to NaN",
        "masks": {
            "matlab_diag": "finite empirical peak delay, 0 < delay <= 200 ms, N_stimulations x probability >= 100, "
                           "finite TTP, diagonal (self-pairs) INCLUDED as in RC_delay_scatter_plot_Amp.m / _Prob.m "
                           "(n = 2,320 matches the old figure); unconnected = F <= 0.25",
            "matlab": "matlab_diag without the diagonal",
            "fig3": "matlab mask + same hemisphere + ED >= 20 mm + N_implantations >= 5",
        },
        "variants": {
            "xstruct": "X_struct.v8.W / ED (Fig. 3 run)",
            "neuron_F025": "SC_data_neuron_2023_yihang SC / D with SC[F < 0.25] = 0 (RC_delay_scatter_refine.m)",
        },
        "note": "Self-pairs have simulated TTP = 0 and empirical delays of ~26 ms; including them inflates the "
                "'unconnected' correlation (F = 0 on the diagonal) from ~0.2 to ~0.7.",
        "bootstrap": {"n_boot": args.n_boot, "seed": args.seed, "ci": "percentile 2.5 / 97.5"},
        "figure": {"variant": args.plot_variant, "mask": args.plot_mask, "axis_max_ms": args.axis_max},
    }
    (args.output_dir / f"provenance{suffix}.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Wrote {png} (+ .pdf)")


if __name__ == "__main__":
    main()
