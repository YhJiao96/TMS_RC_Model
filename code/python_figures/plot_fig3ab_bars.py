#!/usr/bin/env python3
"""Redraw Fig. 3A/3B (F-TRACT probability / amplitude bar charts) in Python.

Reads the bootstrap export of ``Python/export_fig3a_bootstrap_ci.py``
(``fig3a_bootstrap_ci.csv``) and reproduces the MATLAB figure
(``yihang_ftract/RC_ftract_boot_subplot.m``) with three changes requested for
the revision:

* the same bar order in all six panels (``--order``);
* predictors that cannot be evaluated in a panel (SC in the unconnected
  panels, where every SC weight is zero) are shown as an empty ``n/a`` slot;
* one y-axis per row, wide enough for every confidence interval.

Bars: translucent = full Spearman |r|, opaque = distance-controlled partial
Spearman |r|; whiskers = 95 % bootstrap CI (B = 1000, paired resampling).
Colours and transparency follow the MATLAB script.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_CSV = PROJECT_ROOT / "output" / "fig3a_bootstrap_ci" / "fig3a_bootstrap_ci.csv"

COLOURS = {  # from RC_ftract_boot_subplot.m
    "RC": (0.8392, 0.1529, 0.1569),
    "SI": (0.1216, 0.4667, 0.7059),
    "SPE": (0.0902, 0.7451, 0.8118),
    "ED": (0.9294, 0.6941, 0.1255),
    "SC": (0.5490, 0.3373, 0.2941),
    "CMY": (0.5804, 0.4039, 0.7412),
    "NE": (1.0000, 0.4980, 0.0549),
    "DE": (0.1725, 0.6275, 0.1725),
}
CSV_LABEL = {"RC": "RC", "SI": "SI", "SPE": "SPE", "ED": "ED", "SC": "W", "CMY": "CMY", "NE": "NE", "DE": "DE"}
TICK_LABEL = {"RC": "RC", "SI": "SI", "SPE": "SPE", "ED": "ED", "SC": "SC", "CMY": "CMY", "NE": "NE", "DE": "DE"}
PAIR_SETS = ("All pairs", "Connected pairs", "Unconnected pairs")
ROWS = (
    ("Probability", "A", "Prob: Full vs Partial correlation (95% bootstrap CI)"),
    ("Amplitude", "B", "Amp: Full vs Partial correlation (95% bootstrap CI)"),
)
FULL_ALPHA, PART_ALPHA, BAR_WIDTH, ERR_DX = 0.30, 0.95, 0.72, 0.08


def fnum(text: str) -> float:
    return float(text) if text not in ("", "nan") else math.nan


def abs_with_ci(r: float, lo: float, hi: float) -> tuple[float, float, float]:
    """|r| and its CI; a negative correlation flips the interval."""
    if math.isnan(r):
        return math.nan, math.nan, math.nan
    if r < 0:
        return -r, -hi, -lo
    return r, lo, hi


def load(csv_path: Path) -> dict[tuple[str, str], dict[str, dict[str, float]]]:
    table: dict[tuple[str, str], dict[str, dict[str, float]]] = {}
    for row in csv.DictReader(csv_path.open()):
        key = (row["outcome"], row["pair_set"])
        full, full_lo, full_hi = abs_with_ci(fnum(row["r_full"]), fnum(row["r_full_ci_low"]), fnum(row["r_full_ci_high"]))
        part, part_lo, part_hi = abs_with_ci(fnum(row["r_partial_ed"]), fnum(row["r_partial_ci_low"]), fnum(row["r_partial_ci_high"]))
        table.setdefault(key, {})[row["label"]] = {
            "n": int(row["n_pairs"]), "full": full, "full_lo": full_lo, "full_hi": full_hi,
            "part": part, "part_lo": part_lo, "part_hi": part_hi,
        }
    return table


def draw_panel(ax: plt.Axes, panel: dict[str, dict[str, float]], order: list[str], ylim: float,
               title: str, show_legend: bool, font: float) -> None:
    for x, metric in enumerate(order):
        stats = panel.get(CSV_LABEL[metric])
        colour = COLOURS[metric]
        if stats is None or math.isnan(stats["full"]):
            ax.text(x, 0.02 * ylim, "n/a", ha="center", va="bottom", fontsize=font, color="0.45")
            continue
        ax.bar(x, stats["full"], BAR_WIDTH, color=colour, alpha=FULL_ALPHA, edgecolor="none", zorder=2)
        ax.errorbar(x - ERR_DX, stats["full"],
                    yerr=[[stats["full"] - stats["full_lo"]], [stats["full_hi"] - stats["full"]]],
                    fmt="none", ecolor="black", elinewidth=1.0, capsize=0, zorder=4)
        if not math.isnan(stats["part"]):
            ax.bar(x, stats["part"], BAR_WIDTH, color=colour, alpha=PART_ALPHA, edgecolor="none", zorder=3)
            ax.errorbar(x + ERR_DX, stats["part"],
                        yerr=[[stats["part"] - stats["part_lo"]], [stats["part_hi"] - stats["part"]]],
                        fmt="none", ecolor="black", elinewidth=1.0, capsize=0, zorder=4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([TICK_LABEL[m] for m in order], fontsize=font)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(0, ylim)
    ax.set_yticks([round(v, 1) for v in frange(0.0, ylim + 1e-9, 0.1)])
    ax.tick_params(axis="y", labelsize=font)
    ax.set_title(title, fontsize=font + 1, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.2)
    ax.grid(axis="y", alpha=0.2, zorder=0)
    if show_legend:
        # legend above the panel, so it never overlaps a whisker
        handles = [Patch(facecolor="0.3", alpha=FULL_ALPHA, edgecolor="none", label="Full (translucent)"),
                   Patch(facecolor="0.3", alpha=PART_ALPHA, edgecolor="none", label="Partial (opaque)")]
        ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.09), ncol=1, frameon=False,
                  fontsize=font - 1, handlelength=1.3, labelspacing=0.25, borderaxespad=0.0)


def frange(start: float, stop: float, step: float) -> list[float]:
    values, v = [], start
    while v <= stop:
        values.append(v)
        v += step
    return values


def row_ylim(table, outcome: str, order: list[str], forced: float | None) -> float:
    if forced is not None:
        return forced
    top = 0.0
    for pair_set in PAIR_SETS:
        for metric in order:
            stats = table[(outcome, pair_set)].get(CSV_LABEL[metric])
            if stats is None:
                continue
            for key in ("full_hi", "part_hi"):
                if not math.isnan(stats[key]):
                    top = max(top, stats[key])
    return math.ceil((top + 0.015) * 10) / 10


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--order", type=str, default="ED,SI,RC,SPE,SC",
                        help="comma-separated bar order used in every panel")
    parser.add_argument("--ylim-prob", type=float, default=None)
    parser.add_argument("--ylim-amp", type=float, default=None)
    parser.add_argument("--font-size", type=float, default=17.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_CSV.parent / "fig3ab_bars.png",
                        help="PNG path; a PDF with the same stem is written as well")
    args = parser.parse_args()
    order = [m.strip() for m in args.order.split(",") if m.strip()]
    unknown = [m for m in order if m not in COLOURS]
    if unknown:
        raise SystemExit(f"Unknown metrics in --order: {unknown}; known: {sorted(COLOURS)}")

    table = load(args.csv)
    ylims = {"Probability": row_ylim(table, "Probability", order, args.ylim_prob),
             "Amplitude": row_ylim(table, "Amplitude", order, args.ylim_amp)}

    fig, axes = plt.subplots(2, 3, figsize=(16, 10.5))
    fig.subplots_adjust(left=0.065, right=0.99, top=0.905, bottom=0.055, hspace=0.58, wspace=0.24)
    for r, (outcome, letter, row_title) in enumerate(ROWS):
        for c, pair_set in enumerate(PAIR_SETS):
            panel = table[(outcome, pair_set)]
            n = next(iter(panel.values()))["n"]
            draw_panel(axes[r, c], panel, order, ylims[outcome], f"{pair_set} (n = {n:,})",
                       show_legend=(c == 2), font=args.font_size)
            if c == 0:
                axes[r, c].set_ylabel("|Spearman r|", fontsize=args.font_size + 1)
        # row title centred over the first two panels; the legend sits above the third
        pos0, pos1 = axes[r, 0].get_position(), axes[r, 1].get_position()
        y_title = pos0.y1 + 0.052
        fig.text((pos0.x0 + pos1.x1) / 2, y_title, row_title, fontsize=args.font_size + 3,
                 fontweight="bold", ha="center", va="bottom")
        fig.text(pos0.x0 - 0.055, y_title, letter, fontsize=args.font_size + 9, fontweight="bold",
                 ha="left", va="bottom")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, facecolor="white")
    fig.savefig(args.output.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.output} (+ .pdf); y-limits: {ylims}; order: {order}")


if __name__ == "__main__":
    main()
