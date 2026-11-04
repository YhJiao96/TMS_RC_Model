#!/usr/bin/env python3
"""Compose the manuscript Fig. 5: (A, B) simulated time-to-peak vs empirical onset
delay for all / structurally unconnected pairs, (C) delay predictors compared.

Pair set (Fig. 3 criteria + delay criteria, self-pairs excluded): same hemisphere,
ED >= 20 mm, N_implantations >= 5, N_responses >= 100, 0 < delay <= 200 ms,
finite simulated TTP; unconnected = F <= 0.25.  The bar row and the scatter
annotations both read the bootstrap CSV written by ``fig5_delay_bars.py`` so the
figure quotes a single bootstrap.

Output: ``<results-dir>/fig5_full.png`` (300 dpi) and ``fig5_full.pdf``.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig5_delay_bars  # noqa: E402,F401  (registers the SPL colour and labels)
import plot_fig3ab_bars as bars  # noqa: E402
import rc_ftract_core as rcc  # noqa: E402
from fig5_ttp_delay import DELAY_XLABEL, RC_RED, clean_ttp, draw_panel  # noqa: E402

# manuscript wording for the bar-row title, keyed by the "outcome" column of the bars CSV
OUTCOME_LABEL = {"Onset delay": "Onset latency", "Delay": "Peak latency", "Latency start": "Latency start"}

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_RESULTS = PROJECT_ROOT / "output" / "fig5_ttp" / "final_onset_neuron"
DEFAULT_DATA_ROOT = rcc.DEFAULT_DATA_ROOT


SUBSETS = ("All pairs", "Connected pairs", "Unconnected pairs")


def load_spin(spin_dir: Path) -> dict:
    rows = {}
    for row in csv.DictReader((spin_dir / "spin_results.csv").open()):
        rows[(row["metric"], row["pair_set"])] = {k: (float(v) if k not in ("metric", "pair_set") else v) for k, v in row.items()}
    nulls = dict(np.load(spin_dir / "spin_nulls.npz"))
    return {"rows": rows, "nulls": nulls}


def draw_spin_panel(ax, spin: dict, font: float) -> None:
    """RC observed rho / partial rho against the spin null for the three pair sets."""
    pos, ticks, labels, obs_all, null_lo, groups = 0.0, [], [], [], [], []
    gap = 4.2   # spacing between the three pair-set groups (keeps the bold group labels apart at large fonts)
    for subset in SUBSETS:
        key = subset.replace(" ", "_")
        start = pos
        for stat, label in (("full", "full"), ("partial", "\npartial")):
            null = spin["nulls"][f"RC__{key}__{stat}"]
            row = spin["rows"][("RC", subset)]
            obs = row["rho"] if stat == "full" else row["partial"]
            pval = row["p_spin_rho"] if stat == "full" else row["p_spin_partial"]
            body = ax.violinplot([null], positions=[pos], widths=0.85, showextrema=False)["bodies"][0]
            body.set_facecolor("0.72")
            body.set_edgecolor("none")
            body.set_alpha(1.0)
            ax.plot(pos, obs, "o", color=RC_RED, markersize=9, zorder=5)
            # full: label centred above the dot; partial: label to the right of the dot, clear of the neighbouring violin
            if stat == "full":
                ax.text(pos, obs + 0.014, f"p =\n{pval:.3f}" if pval >= 0.001 else "p <\n0.001", ha="center", va="bottom",
                        fontsize=font - 2, color="black", linespacing=0.95)
            else:
                ax.text(pos + 0.22, obs + 0.004, f"p =\n{pval:.3f}" if pval >= 0.001 else "p <\n0.001", ha="left", va="bottom",
                        fontsize=font - 2, color="black", linespacing=0.95)
            ticks.append(pos)
            labels.append(label)
            obs_all.append(obs)
            null_lo.append(null.min())
            pos += 1.15
        groups.append(((start + pos - 1.15) / 2, subset.replace(" pairs", "")))
        pos += gap
    ax.axhline(0, color="0.5", linewidth=1.0, zorder=1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=font - 4, linespacing=0.9)
    ax.set_xticks([g[0] for g in groups], minor=True)
    ax.set_xticklabels([g[1] for g in groups], minor=True, fontsize=font - 8, fontweight="bold")
    ax.tick_params(axis="x", which="minor", length=0, pad=38)
    ax.set_xlim(-1.6, pos - gap + 2.2)   # room for the first p-value label on the left and the last partial label on the right
    ax.set_ylim(min(null_lo) - 0.10, max(obs_all) + 0.15)   # room below the violins for the legend and above for the two-line p labels
    ax.set_ylabel("Spearman r (RC model)", fontsize=font + 1)
    ax.tick_params(axis="y", labelsize=font - 2)
    ax.set_title("Spin null\n(1,000 rotations)", fontsize=font + 1, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="0.72", edgecolor="none", label="spin null"),
                       Line2D([], [], marker="o", color=RC_RED, linestyle="none", markersize=8, label="observed")],
              loc="lower left", ncol=2, columnspacing=1.0, frameon=False, fontsize=font - 4, handlelength=1.2)


def mark_spin(ax, panel: dict, order: list[str], rows: dict, pair_set: str, font: float) -> None:
    """asterisk above a bar whose correlation exceeds the spin null (p_spin < 0.05); ED is not spun."""
    for x, metric in enumerate(order):
        if (metric, pair_set) not in rows:
            continue
        stats, row = panel[bars.CSV_LABEL[metric]], rows[(metric, pair_set)]
        if row["p_spin_rho"] < 0.05:
            ax.text(x - bars.ERR_DX, stats["full_hi"] + 0.008, "*", ha="center", va="bottom", fontsize=font + 2, zorder=6)
        if row["p_spin_partial"] < 0.05 and not np.isnan(stats["part_hi"]):
            ax.text(x + bars.ERR_DX, stats["part_hi"] + 0.008, "*", ha="center", va="bottom", fontsize=font + 2, zorder=6)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--ttp-npz", type=Path, default=PROJECT_ROOT / "output" / "fig5_ttp" / "ttp_neuron_F025_raw.npz")
    parser.add_argument("--bars-csv", default="fig5_delay_bars_onset_delay_neuron.csv")
    parser.add_argument("--delay-variable", default="onset_delay", choices=tuple(DELAY_XLABEL))
    parser.add_argument("--count-variable", default="responses", choices=("responses", "stimulations"),
                        help="sampling criterion for a pair: N_responses >= 100 (default) or N_stimulations >= 100 (the Fig. 3 probability criterion)")
    parser.add_argument("--order", default="ED,RC,SPL,SI")
    parser.add_argument("--axis-max", type=float, default=200.0)
    parser.add_argument("--x-max", type=float, default=None, help="x-axis limit of the scatter panels (default: --axis-max)")
    parser.add_argument("--y-max", type=float, default=None, help="y-axis limit of the scatter panels (default: --axis-max)")
    parser.add_argument("--point-size", type=float, default=9.0)
    parser.add_argument("--point-alpha", type=float, default=0.25)
    parser.add_argument("--point-color", default="0.35")
    parser.add_argument("--point-style", default="plain", choices=("plain", "density", "hexbin"),
                        help="plain: uniform points; density: points coloured by local pair density; hexbin: hexagonal density map")
    parser.add_argument("--band", default="iqr", choices=("iqr", "ci95"),
                        help="band around the binned median: interquartile range of the pairs (iqr) or 95%% bootstrap CI of the median (ci95)")
    parser.add_argument("--output-name", default=None, help="output file stem (default fig5_full_spin / fig5_full / fig5_scatter_spin)")
    parser.add_argument("--font-size", type=float, default=17.0)
    parser.add_argument("--spin-dir", type=Path, default=None,
                        help="folder written by fig5_spin_test.py; adds a spin-null panel (C) and marks bars with p_spin < 0.05")
    parser.add_argument("--no-bars", action="store_true",
                        help="with --spin-dir: one row only (A/B scatter + C spin null), no bar row; output fig5_scatter_spin.png")
    args = parser.parse_args()
    if args.no_bars and args.spin_dir is None:
        parser.error("--no-bars requires --spin-dir")
    font = args.font_size
    order = [m.strip() for m in args.order.split(",")]

    # --- data for the scatter panels (same masks as fig5_ttp_delay.py / fig5_delay_bars.py, mask "fig3")
    data = rcc.load_comparison_data(args.data_root)
    order_sc = rcc.load_permutation(args.data_root / "ReOrder_Caio_SC.mat", "ReOrder_Caio_SC", 374)
    order_ft = rcc.load_permutation(args.data_root / "ReOreder_F_tract.mat", "ReOrder_F_tract", 360, column=1)
    delay = rcc.read_ftract_gzip(rcc.ftract_dir(args.data_root) / f"{args.delay_variable}__median.txt.gz")
    delay = delay[np.ix_(order_ft, order_ft)]
    with np.load(args.ttp_npz) as npz:
        ttp = clean_ttp(npz["ttp_ms"][np.ix_(order_sc, order_sc)][:360, :360])
    n_resp = data.n_stimulations * data.probability
    hemi = np.zeros((360, 360), bool)
    hemi[:180, :180] = True
    hemi[180:, 180:] = True
    count = n_resp if args.count_variable == "responses" else data.n_stimulations
    fig3 = (np.isfinite(delay) & (delay > 0) & (delay <= 200) & (count >= 100) & ~np.eye(360, dtype=bool)
            & hemi & (data.distance >= 20) & (data.n_implantations >= 5) & np.isfinite(ttp))
    masks = {"All pairs": fig3, "Unconnected pairs": fig3 & (data.frequency <= 0.25)}

    # --- statistics (one bootstrap for the whole figure)
    table = bars.load(args.results_dir / args.bars_csv)
    outcome = next(iter(table))[0]
    outcome_label = OUTCOME_LABEL.get(outcome, outcome)  # row title wording (CSV keeps the original outcome name)
    for pair_set, mask in masks.items():
        n_csv = table[(outcome, pair_set)]["RC"]["n"]
        assert mask.sum() == n_csv, f"{pair_set}: mask n = {mask.sum()} but CSV n = {n_csv}"

    spin = load_spin(args.spin_dir) if args.spin_dir else None
    if spin is None:
        fig = plt.figure(figsize=(16, 14.2))
        outer = fig.add_gridspec(2, 1, height_ratios=[6.6, 4.3], left=0.065, right=0.99, top=0.935, bottom=0.055, hspace=0.34)
        top = outer[0, 0].subgridspec(1, 2, wspace=0.28)
    elif args.no_bars:
        fig = plt.figure(figsize=(16, 5.9))
        top = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.15], left=0.06, right=0.99, top=0.88, bottom=0.165, wspace=0.34)
    else:
        # uniform 2 x 3 grid: equal columns in both rows, equal row heights (the square scatter panels fill their cells)
        fig = plt.figure(figsize=(16, 12.6))
        outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], left=0.07, right=0.99, top=0.93, bottom=0.06, hspace=0.62)
        top = outer[0, 0].subgridspec(1, 3, wspace=0.36)
    narrow = spin is not None and not args.no_bars   # three panels per row: two-line titles, short axis labels
    for c, (pair_set, letter) in enumerate(zip(masks, "AB")):
        ax = fig.add_subplot(top[0, c])
        rc = table[(outcome, pair_set)]["RC"]
        result = {"rho": rc["full"], "rho_lo": rc["full_lo"], "rho_hi": rc["full_hi"],
                  "partial": rc["part"], "partial_lo": rc["part_lo"], "partial_hi": rc["part_hi"]}
        mask = masks[pair_set]
        title = f"{pair_set}\n(n = {rc['n']:,})" if narrow else f"{pair_set} (n = {rc['n']:,})"
        draw_panel(ax, delay[mask], ttp[mask], result, title, font, args.axis_max,
                   xlabel=("Onset latency (ms)" if narrow else DELAY_XLABEL[args.delay_variable]), x_max=args.x_max, y_max=args.y_max,
                   ylabel=("Simulated TTP (ms)" if narrow else "Simulated time-to-peak (ms)"),
                   point_size=args.point_size, point_alpha=args.point_alpha, point_color=args.point_color,
                   point_style=args.point_style, band=args.band)
        if spin is not None:
            ax.get_legend().remove()  # narrower panels: the binned-median / identity legend goes in the caption
        ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontsize=font + 9, fontweight="bold", ha="left", va="bottom")

    if spin is not None:
        ax_spin = fig.add_subplot(top[0, 2])
        draw_spin_panel(ax_spin, spin, font)
        ax_spin.text(-0.16, 1.04, "C", transform=ax_spin.transAxes, fontsize=font + 9, fontweight="bold", ha="left", va="bottom")
    if args.no_bars:
        out = args.results_dir / f"{args.output_name or 'fig5_scatter_spin'}.png"
        fig.savefig(out, dpi=300, facecolor="white")
        fig.savefig(out.with_suffix(".pdf"), facecolor="white")
        plt.close(fig)
        print(f"Wrote {out} (+ .pdf)")
        return
    bottom = outer[1, 0].subgridspec(1, 3, wspace=0.24 if spin is None else 0.34)
    axes = [fig.add_subplot(bottom[0, c]) for c in range(3)]
    ylim = bars.row_ylim(table, outcome, order, None)
    for c, pair_set in enumerate(bars.PAIR_SETS):
        panel = table[(outcome, pair_set)]
        n = next(iter(panel.values()))["n"]
        bars.draw_panel(axes[c], panel, order, ylim, f"{pair_set}\n(n = {n:,})" if narrow else f"{pair_set} (n = {n:,})", show_legend=(c == 2), font=font)
        if spin is not None:
            mark_spin(axes[c], panel, order, spin["rows"], pair_set, font)
        if narrow and c == 2 and axes[c].get_legend() is not None:
            # re-create the bar legend as one row on the row-title line (above the two-line panel title)
            leg = axes[c].get_legend()
            handles = getattr(leg, "legend_handles", None) or getattr(leg, "legendHandles")
            texts = [t.get_text() for t in leg.get_texts()]
            leg.remove()
            axes[c].legend(handles=handles, labels=texts, loc="lower right", bbox_to_anchor=(1.0, 1.25), ncol=2,
                           frameon=False, fontsize=font - 1, handlelength=1.3, columnspacing=1.0, borderaxespad=0.0)
    axes[0].set_ylabel("|Spearman r|", fontsize=font + 1)
    pos0, pos1, pos2 = axes[0].get_position(), axes[1].get_position(), axes[2].get_position()
    if spin is None:
        y_title = pos0.y1 + 0.030
        fig.text((pos0.x0 + pos1.x1) / 2, y_title, f"{outcome_label}: Full vs Partial correlation (95% bootstrap CI)",
                 fontsize=font + 3, fontweight="bold", ha="center", va="bottom")
        fig.text(pos0.x0 - 0.055, y_title, "C", fontsize=font + 9, fontweight="bold", ha="left", va="bottom")
    else:
        # letter in the same place as A-C; row title centred over the whole row, above the panel titles
        axes[0].text(-0.16, 1.04, "D", transform=axes[0].transAxes, fontsize=font + 9, fontweight="bold", ha="left", va="bottom")
        fig.text((pos0.x0 + pos1.x1) / 2, pos0.y1 + 0.085, f"{outcome_label}: full vs partial correlation",
                 fontsize=font + 3, fontweight="bold", ha="center", va="bottom")

    out = args.results_dir / f"{args.output_name or ('fig5_full' if spin is None else 'fig5_full_spin')}.png"
    fig.savefig(out, dpi=300, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(f"Wrote {out} (+ .pdf); bar ylim {ylim}")


if __name__ == "__main__":
    main()
