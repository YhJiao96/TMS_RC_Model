#!/usr/bin/env python3
"""Fig. 5 alternative: bar comparison of delay predictors in the style of Fig. 3A.

Predictors of the empirical F-TRACT delay-to-peak: Euclidean distance (ED),
simulated RC time-to-peak (RC), weighted shortest-path length (SPL = 1/SPE)
and search information cost (SI).  Full Spearman |r| (translucent) and
distance-controlled partial |r| (opaque) with 95 % bootstrap CIs
(B = 1000, seed 17), for all / connected / unconnected pairs under the Fig. 3
inclusion criteria plus 0 < delay <= 200 ms, N_responses >= 100, no self-pairs.

Outputs: ``fig5_delay_bars.csv`` (same columns as fig3a_bootstrap_ci.csv),
``fig5_delay_bars.png`` / ``.pdf``.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402
from scipy.stats import pearsonr, rankdata  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_fig3ab_bars as bars  # noqa: E402
import rc_ftract_core as rcc  # noqa: E402
from export_fig3a_bootstrap_ci import load_comparators  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_DATA_ROOT = rcc.DEFAULT_DATA_ROOT
FIG3_DIR = PROJECT_ROOT / "output" / "fig3a_bootstrap_ci"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "fig5_ttp"

bars.COLOURS["SPL"] = bars.COLOURS["SPE"]
bars.CSV_LABEL["SPL"] = "SPL"
bars.TICK_LABEL["SPL"] = "SPL"

# F-TRACT latency variables (MNI-HCP-MMP1/<variable>__median.txt.gz)
DELAY_LABEL = {"peak_delay": "Delay", "onset_delay": "Onset delay", "latency_start": "Latency start"}


def boot_stats(x: np.ndarray, y: np.ndarray, d: np.ndarray, n_boot: int, rng, with_partial: bool) -> dict:
    full = pearsonr(rankdata(x), rankdata(y)).statistic
    part = rcc.partial_spearman(x, y, d) if with_partial else np.nan
    bf, bp = np.empty(n_boot), np.full(n_boot, np.nan)
    n = x.size
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        bf[b] = pearsonr(rankdata(x[idx]), rankdata(y[idx])).statistic
        if with_partial:
            bp[b] = rcc.partial_spearman(x[idx], y[idx], d[idx])
    return {"r_full": full, "r_full_ci_low": np.percentile(bf, 2.5), "r_full_ci_high": np.percentile(bf, 97.5),
            "r_partial_ed": part, "r_partial_ci_low": np.nanpercentile(bp, 2.5) if with_partial else np.nan,
            "r_partial_ci_high": np.nanpercentile(bp, 97.5) if with_partial else np.nan}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--rc-npz", type=Path, default=FIG3_DIR / "rc_peak_xstruct_w.npz")
    parser.add_argument("--bct-file", type=Path, default=rcc.comparator_dir(DEFAULT_DATA_ROOT) / "bct_measures.mat")
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--order", type=str, default="ED,RC,SPL,SI")
    parser.add_argument("--mask", default="fig3", choices=("fig3", "offdiag", "offdiag_hemi", "no_ed"),
                        help="fig3: Fig. 3 criteria; no_ed: fig3 without ED >= 20 mm; offdiag_hemi: off-diagonal + same hemisphere; "
                             "offdiag: off-diagonal only (plus delay range and N_responses >= 100 in every case)")
    parser.add_argument("--count-variable", default="responses", choices=("responses", "stimulations"),
                        help="sampling criterion for a pair: N_responses >= 100 (default) or N_stimulations >= 100 (the Fig. 3 probability criterion)")
    parser.add_argument("--delay-variable", default="peak_delay", choices=tuple(DELAY_LABEL),
                        help="empirical latency: peak_delay (time to CCEP peak), onset_delay (threshold crossing) "
                             "or latency_start (last local minimum before the peak)")
    parser.add_argument("--unconnected-def", default="F025", choices=("F025", "W0"), help="unconnected pairs: F025 = structural connection probability F <= 0.25 (connected F >= 0.75); W0 = no direct edge in the group connectome, W = 0 (connected W > 0)")
    parser.add_argument("--xstruct-file", type=Path, default=None)
    parser.add_argument("--suffix", default="", help="extra suffix appended to the output file names")
    parser.add_argument("--font-size", type=float, default=17.0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    order = [m.strip() for m in args.order.split(",")]
    rng = np.random.default_rng(args.seed)

    data = rcc.load_comparison_data(args.data_root)
    order_sc = rcc.load_permutation(args.data_root / "ReOrder_Caio_SC.mat", "ReOrder_Caio_SC", 374)
    order_ft = rcc.load_permutation(args.data_root / "ReOreder_F_tract.mat", "ReOrder_F_tract", 360, column=1)
    delay = rcc.read_ftract_gzip(rcc.ftract_dir(args.data_root) / f"{args.delay_variable}__median.txt.gz")
    delay = delay[np.ix_(order_ft, order_ft)]
    comp, _ = load_comparators(rcc.comparator_dir(args.data_root) / "Precomputed_Models.mat", args.bct_file, order_sc)
    with np.load(args.rc_npz) as npz:
        ttp = npz["ttp_ms"][np.ix_(order_sc, order_sc)][:360, :360].astype(float)
    ttp[ttp >= 200] = np.nan
    spe = comp["SPE"]
    spl = np.where(spe > 0, 1.0 / np.where(spe > 0, spe, 1.0), np.nan)
    predictors = {"ED": data.distance, "RC": ttp, "SPL": spl, "SI": np.abs(comp["SI"])}

    nresp = data.n_stimulations * data.probability
    eye = np.eye(360, dtype=bool)
    hemi = np.zeros((360, 360), bool)
    hemi[:180, :180] = True
    hemi[180:, 180:] = True
    count = nresp if args.count_variable == "responses" else data.n_stimulations
    core = np.isfinite(delay) & (delay > 0) & (delay <= 200) & (count >= 100) & ~eye
    base = {
        "fig3": core & hemi & (data.distance >= 20) & (data.n_implantations >= 5),
        "no_ed": core & hemi & (data.n_implantations >= 5),
        "offdiag_hemi": core & hemi,
        "offdiag": core,
    }[args.mask]
    if args.unconnected_def == "W0":
        xs_file = args.xstruct_file or (rcc.xstruct_path(args.data_root))
        xs = rcc.load_xstruct_v8(xs_file)
        W = np.asarray(xs.W, dtype=float)[np.ix_(order_sc, order_sc)][:360, :360]
        connected, unconnected = W > 0, W == 0
    else:
        connected, unconnected = data.frequency >= 0.75, data.frequency <= 0.25
    subsets = {"All pairs": base, "Connected pairs": base & connected, "Unconnected pairs": base & unconnected}

    rows = []
    table: dict[tuple[str, str], dict[str, dict[str, float]]] = {}
    for pair_set, mask in subsets.items():
        # one common pair set per panel: every predictor finite (TTP >= 200 ms is NaN for 2 pairs)
        common = mask & np.isfinite(delay)
        for metric in order:
            common &= np.isfinite(predictors[metric])
        for metric in order:
            x = predictors[metric]
            v = common
            res = boot_stats(x[v], delay[v], data.distance[v], args.n_boot, rng, with_partial=(metric != "ED"))
            rows.append({"outcome": DELAY_LABEL[args.delay_variable], "pair_set": pair_set, "n_pairs": int(v.sum()), "metric": metric,
                         "label": bars.CSV_LABEL[metric], **res})
            full, lo, hi = bars.abs_with_ci(res["r_full"], res["r_full_ci_low"], res["r_full_ci_high"])
            part, plo, phi = bars.abs_with_ci(res["r_partial_ed"], res["r_partial_ci_low"], res["r_partial_ci_high"])
            table.setdefault(("Delay", pair_set), {})[bars.CSV_LABEL[metric]] = {
                "n": int(v.sum()), "full": full, "full_lo": lo, "full_hi": hi, "part": part, "part_lo": plo, "part_hi": phi}
            print(f"{pair_set:18s} {metric:4s} n={v.sum():5d} |r|={full:.3f} [{lo:.3f},{hi:.3f}] "
                  f"partial={part:.3f} [{plo:.3f},{phi:.3f}]")
    suffix = ("" if args.count_variable == "responses" else "_nstim") + \
        ("" if args.mask == "fig3" else f"_{args.mask}") + \
        ("" if args.delay_variable == "peak_delay" else f"_{args.delay_variable}") + \
        ("" if args.unconnected_def == "F025" else "_W0") + args.suffix
    with (args.output_dir / f"fig5_delay_bars{suffix}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    top = max(max(s["full_hi"], 0 if np.isnan(s["part_hi"]) else s["part_hi"])
              for panel in table.values() for s in panel.values())
    ylim = float(np.ceil((top + 0.015) * 10) / 10)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6))
    fig.subplots_adjust(left=0.065, right=0.99, top=0.80, bottom=0.11, wspace=0.24)
    for ax, pair_set in zip(axes, subsets):
        panel = table[("Delay", pair_set)]
        n = next(iter(panel.values()))["n"]
        bars.draw_panel(ax, panel, order, ylim, f"{pair_set} (n = {n:,})", show_legend=(pair_set == "Unconnected pairs"),
                        font=args.font_size)
    axes[0].set_ylabel("|Spearman r|", fontsize=args.font_size + 1)
    pos0, pos1 = axes[0].get_position(), axes[1].get_position()
    fig.text((pos0.x0 + pos1.x1) / 2, pos0.y1 + 0.09, f"{DELAY_LABEL[args.delay_variable]}: Full vs Partial correlation (95% bootstrap CI)",
             fontsize=args.font_size + 3, fontweight="bold", ha="center", va="bottom")
    png = args.output_dir / f"fig5_delay_bars{suffix}.png"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(png.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(f"Wrote {png} (+ .pdf); ylim {ylim}")


if __name__ == "__main__":
    main()
