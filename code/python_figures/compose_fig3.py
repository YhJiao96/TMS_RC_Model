#!/usr/bin/env python3
"""Compose the full manuscript Fig. 3: A/B bar charts, C per-site map, D distribution.

A, B : drawn natively with ``plot_fig3ab_bars`` from the bootstrap CI export.
C    : the four cropped Cerebro views saved by ``render_fig3c_site_correlation.py``
       (``<map-name>_crops.npy``), plus colour bar and category legend.
D    : site-wise Spearman r vs. number of eligible targets, with a marginal
       histogram of r, from ``site_correlations.csv``.

Output: ``<output-dir>/fig3_full.png`` (300 dpi) and ``fig3_full.pdf``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

import plot_fig3ab_bars as bars  # noqa: E402
from render_fig3c_site_correlation import CMAP_NAME, NS_COLOR, legend_handles  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_BARS_CSV = PROJECT_ROOT / "output" / "fig3a_bootstrap_ci" / "fig3a_bootstrap_ci.csv"
DEFAULT_SITE_DIR = PROJECT_ROOT / "output" / "fig3c_site_corr"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "fig3_composite"


def load_sites(site_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rho, count, category = [], [], []
    for row in csv.DictReader((site_dir / "site_correlations.csv").open()):
        rho.append(float(row["spearman_rho"]) if row["spearman_rho"] else math.nan)
        count.append(int(row["n_eligible_targets"]))
        category.append(row["category"])
    return np.asarray(rho), np.asarray(count), np.asarray(category, dtype=object)


def draw_bar_rows(fig, gs_rows, table, order, font: float) -> None:
    ylims = {outcome: bars.row_ylim(table, outcome, order, None) for outcome, _, _ in bars.ROWS}
    for r, (outcome, letter, row_title) in enumerate(bars.ROWS):
        axes = [fig.add_subplot(gs_rows[r][0, c]) for c in range(3)]
        for c, pair_set in enumerate(bars.PAIR_SETS):
            panel = table[(outcome, pair_set)]
            n = next(iter(panel.values()))["n"]
            bars.draw_panel(axes[c], panel, order, ylims[outcome], f"{pair_set} (n = {n:,})",
                            show_legend=(c == 2), font=font)
        axes[0].set_ylabel("|Spearman r|", fontsize=font + 1)
        pos0, pos1 = axes[0].get_position(), axes[1].get_position()
        y_title = pos0.y1 + 0.030
        fig.text((pos0.x0 + pos1.x1) / 2, y_title, row_title, fontsize=font + 3, fontweight="bold",
                 ha="center", va="bottom")
        fig.text(pos0.x0 - 0.055, y_title, letter, fontsize=font + 9, fontweight="bold", ha="left", va="bottom")


def draw_map(fig, cell, crops: list[np.ndarray], rho_lims, counts, min_targets, alpha, font: float):
    heights = [max(crops[0].shape[0], crops[1].shape[0]), max(crops[2].shape[0], crops[3].shape[0])]
    widths = [max(crops[0].shape[1], crops[2].shape[1]), max(crops[1].shape[1], crops[3].shape[1])]
    inner = cell.subgridspec(4, 1, height_ratios=[heights[0], heights[1], 0.24 * sum(heights), 0.10 * sum(heights)],
                             hspace=0.0)
    views = inner[0:2, 0].subgridspec(2, 2, height_ratios=heights, width_ratios=widths, hspace=0.01, wspace=0.01)
    for index, crop in enumerate(crops):
        axis = fig.add_subplot(views[index // 2, index % 2])
        axis.imshow(crop)
        axis.axis("off")
    # colour bar centred under the views
    holder = fig.add_subplot(inner[2, 0])
    holder.axis("off")
    pos = holder.get_position()
    cax = fig.add_axes([pos.x0 + 0.12 * pos.width, pos.y0 + 0.62 * pos.height, 0.76 * pos.width, 0.20 * pos.height])
    scalar = plt.cm.ScalarMappable(norm=Normalize(*rho_lims), cmap=plt.get_cmap(CMAP_NAME))
    scalar.set_array([])
    colorbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
    colorbar.set_label("Spearman r (significant sites)", fontsize=font, labelpad=4)
    colorbar.ax.tick_params(labelsize=font - 1, width=1.4, length=6)
    colorbar.outline.set_linewidth(1.4)
    legend_ax = fig.add_subplot(inner[3, 0])
    legend_ax.axis("off")
    handles = legend_handles(counts, min_targets, alpha)
    for handle, label in zip(handles, (f"not tested (n = {counts['not tested']})",
                                       f"n.s. (n = {counts['n.s.']})",
                                       f"FDR q < {alpha:g} (n = {counts['significant']})")):
        handle.set_label(label)
    legend_ax.legend(handles=handles, loc="center", ncol=3, frameon=False, fontsize=font - 3,
                     handlelength=1.5, handleheight=1.0, columnspacing=1.6)


def draw_panel_d(fig, cell, rho: np.ndarray, count: np.ndarray, category: np.ndarray, alpha: float,
                 font: float):
    """Panel D: site-wise r against the number of eligible targets, with a
    marginal histogram of r on the right. Returns (scatter axis, histogram axis)."""
    tested = np.isfinite(rho)
    sig = category == "significant"
    ns = tested & ~sig
    n_tested, n_sig = int(tested.sum()), int(sig.sum())
    median = float(np.nanmedian(rho[tested]))
    sig_color = plt.get_cmap(CMAP_NAME)(0.6)
    inner = cell.subgridspec(1, 2, width_ratios=[2.9, 1.0], wspace=0.06)
    ax = fig.add_subplot(inner[0, 0])
    ax_h = fig.add_subplot(inner[0, 1], sharey=ax)

    ax.scatter(count[ns], rho[ns], s=34, color=NS_COLOR, label=f"n.s. (n = {n_tested - n_sig})", zorder=3)
    ax.scatter(count[sig], rho[sig], s=34, color=sig_color, label=f"FDR q < {alpha:g} (n = {n_sig})", zorder=2)
    ax.axhline(0.0, color="grey", linewidth=0.9, zorder=1)
    ax.axhline(median, color="black", linestyle="--", linewidth=1.3, label=f"median r = {median:.2f}", zorder=1)
    ax.set_xlabel("Eligible targets per stimulation site", fontsize=font + 1)
    ax.set_ylabel("Spearman r per stimulation site", fontsize=font + 1)
    ax.set_xlim(0, np.ceil(count[tested].max() / 10) * 10 + 5)
    ax.set_ylim(-0.4, 1.02)
    handles, labels = ax.get_legend_handles_labels()
    order = [1, 0, 2]  # significant, n.s., median
    ax.legend([handles[i] for i in order], [labels[i] for i in order], frameon=False, fontsize=font - 2,
              loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.2)
    ax.tick_params(labelsize=font - 1)

    lo = np.floor(np.nanmin(rho[tested]) * 10) / 10 - 0.05
    bins = np.arange(min(lo, 0.0), 1.0001, 0.05)
    ax_h.hist([rho[sig], rho[ns]], bins=bins, stacked=True, orientation="horizontal",
              color=[sig_color, NS_COLOR], edgecolor="white", linewidth=0.6)
    ax_h.axhline(median, color="black", linestyle="--", linewidth=1.3)
    ax_h.set_xlabel("Sites", fontsize=font + 1)
    ax_h.set_xlim(0, np.ceil(ax_h.get_xlim()[1] / 10) * 10 + 4)
    ax_h.set_xticks(np.arange(0, ax_h.get_xlim()[1], 20))
    ax_h.tick_params(labelsize=font - 1, labelleft=False, left=False)
    ax_h.spines[["top", "right", "left"]].set_visible(False)
    ax_h.spines["bottom"].set_linewidth(1.2)
    return ax, ax_h


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--bars-csv", type=Path, default=DEFAULT_BARS_CSV)
    parser.add_argument("--site-dir", type=Path, default=DEFAULT_SITE_DIR)
    parser.add_argument("--map-name", type=str, default="fig3c_site_rho_map_cortex")
    parser.add_argument("--order", type=str, default="ED,SI,RC,SPE,SC")
    parser.add_argument("--font-size", type=float, default=17.0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    order = [m.strip() for m in args.order.split(",") if m.strip()]
    font = args.font_size

    table = bars.load(args.bars_csv)
    provenance = json.loads((args.site_dir / "provenance.json").read_text())
    rho_lims = tuple(provenance["rho_colour_limits"])
    summary = provenance["summary"]
    counts = {"not tested": summary["sites_not_tested"], "n.s.": summary["sites_ns"],
              "significant": summary["sites_significant_positive"]}
    min_targets, alpha = provenance["min_targets_per_site"], provenance["alpha"]
    crops = list(np.load(args.site_dir / f"{args.map_name}_crops.npy", allow_pickle=True))
    rho, count, category = load_sites(args.site_dir)

    fig = plt.figure(figsize=(16, 17.5))
    outer = fig.add_gridspec(3, 1, height_ratios=[4.3, 4.3, 6.4], left=0.065, right=0.99, top=0.925, bottom=0.04,
                             hspace=0.42)
    gs_rows = [outer[r, 0].subgridspec(1, 3, wspace=0.24) for r in range(2)]
    draw_bar_rows(fig, gs_rows, table, order, font)

    bottom = outer[2, 0].subgridspec(1, 2, width_ratios=[1.08, 0.92], wspace=0.10)
    draw_map(fig, bottom[0, 0], crops, rho_lims, counts, min_targets, alpha, font)
    ax_sc, ax_h = draw_panel_d(fig, bottom[0, 1], rho, count, category, alpha, font)
    # shrink panel D vertically so it sits level with the brain views
    # (scatter shifts right to clear the y-label; marginal histogram keeps its right edge)
    pos = ax_sc.get_position()
    ax_sc.set_position([pos.x0 + 0.03, pos.y0 + 0.16 * pos.height, pos.width - 0.03, 0.74 * pos.height])
    pos = ax_h.get_position()
    ax_h.set_position([pos.x0, pos.y0 + 0.16 * pos.height, pos.width, 0.74 * pos.height])

    # panel letters C and D at the top-left of each bottom cell
    y_cd = outer[2, 0].get_position(fig).y1 + 0.004
    x_c = outer[2, 0].get_position(fig).x0 - 0.055
    x_d = bottom[0, 1].get_position(fig).x0 - 0.01
    fig.text(x_c, y_cd, "C", fontsize=font + 9, fontweight="bold", ha="left", va="bottom")
    fig.text(x_d, y_cd, "D", fontsize=font + 9, fontweight="bold", ha="left", va="bottom")

    png = args.output_dir / "fig3_full.png"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(png.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(f"Wrote {png} (+ .pdf)")


if __name__ == "__main__":
    main()
