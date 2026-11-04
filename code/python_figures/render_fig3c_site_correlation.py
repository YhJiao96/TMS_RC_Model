#!/usr/bin/env python3
"""Per-stimulation-site correlation map (Fig. 3C) and its distribution (Fig. 3D).

For every cortical stimulation site *i* the simulated RC peak current across
target parcels (row *i* of the peak-current matrix, used un-transposed against
F-TRACT exactly as in Fig. 3A) is correlated (Spearman) with the empirical
F-TRACT response-probability row *i*, using the Fig. 3A pair-inclusion mask
(same hemisphere, off-diagonal, ED >= 20 mm, N_implantations >= 5,
N_stimulations >= 100, finite probability).  Sites with fewer than
``--min-targets`` eligible targets are *not tested*; the remaining p-values are
Benjamini-Hochberg corrected across sites.

Outputs (``--output-dir``):
  fig3c_site_rho_map.png            four-view Cerebro map (Fig. 2 pipeline):
                                    light grey = not tested, dark grey = n.s.,
                                    colour = Spearman rho of significant sites
  fig3d_site_rho_distribution.png   histogram of site-wise rho + rho vs n_targets
  fig3d_site_rho_histogram.png      histogram only
  site_correlations.csv             per-site table
  provenance.json, README.md

The RC matrix defaults to the X_struct.v8.W run behind Fig. 3A
(``results/fig3a_bootstrap_ci_2026-10-28/rc_peak_xstruct_w.npz``).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap, Normalize  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rc_ftract_core as rcc  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_DATA_ROOT = rcc.DEFAULT_DATA_ROOT
DEFAULT_RC_NPZ = (
    PROJECT_ROOT / "output" / "fig3a_bootstrap_ci" / "rc_peak_xstruct_w.npz"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "fig3c_site_corr"

UNTESTED_COLOR = "#B8B8B8"  # anatomy layer of the Fig. 2 renderer
NS_COLOR = "#3A3A3A"
NS_BINS = 8  # leading colormap entries reserved for the n.s. colour
CMAP_NAME = "plasma"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_fdr(p: np.ndarray) -> np.ndarray:
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * m / (np.arange(m) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(m)
    tmp[order] = np.minimum(adjusted, 1.0)
    q[ok] = tmp
    return q


def per_site_spearman(
    rc: np.ndarray, emp: np.ndarray, mask: np.ndarray, min_targets: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = rc.shape[0]
    rho = np.full(n, np.nan)
    p = np.full(n, np.nan)
    count = np.zeros(n, dtype=int)
    for i in range(n):
        valid = mask[i] & np.isfinite(emp[i]) & np.isfinite(rc[i])
        count[i] = int(valid.sum())
        if count[i] < min_targets:
            continue
        x, y = rc[i, valid], emp[i, valid]
        if np.ptp(x) == 0 or np.ptp(y) == 0:
            continue
        result = spearmanr(x, y)
        rho[i], p[i] = float(result.statistic), float(result.pvalue)
    return rho, p, count


def categorise(rho: np.ndarray, q: np.ndarray, alpha: float) -> np.ndarray:
    category = np.full(rho.size, "not tested", dtype=object)
    tested = np.isfinite(rho)
    category[tested] = "n.s."
    category[tested & (q < alpha) & (rho > 0)] = "significant"
    category[tested & (q < alpha) & (rho < 0)] = "significant negative"
    return category


def build_colormap(lo: float, hi: float) -> tuple[ListedColormap, tuple[float, float], float]:
    """Colormap whose first NS_BINS entries are the n.s. colour, the rest CMAP_NAME.

    Returns (cmap, display clims, sentinel value for n.s. sites).  A site with
    rho = lo lands on the first colour entry of CMAP_NAME; the sentinel is
    below the display range and is clipped to entry 0 (the n.s. colour).
    """
    n_total = 256
    base = plt.get_cmap(CMAP_NAME)(np.linspace(0.0, 1.0, n_total - NS_BINS))
    ns = np.array(matplotlib.colors.to_rgba(NS_COLOR))[None, :].repeat(NS_BINS, axis=0)
    cmap = ListedColormap(np.vstack([ns, base]))
    step = (hi - lo) / (n_total - NS_BINS - 1)
    display_lo = lo - NS_BINS * step
    sentinel = display_lo - 10.0 * step
    return cmap, (display_lo, hi), sentinel


def render_map(
    renderer,
    values_374: np.ndarray,
    cmap: ListedColormap,
    display_clims: tuple[float, float],
    rho_lims: tuple[float, float],
    path: Path,
    counts: dict[str, int],
    min_targets: int,
    alpha: float,
    title: str | None,
    subcortex: bool = True,
) -> None:
    dense = renderer.to_dense(values_374)
    left = dense.copy()
    left[~renderer.left_mask] = np.nan
    right = dense.copy()
    right[~renderer.right_mask] = np.nan

    # Render every view on its own transparent canvas, crop it to its content
    # and compose the four crops in a tight grid: the raw Cerebro viewport
    # leaves a wide empty band between the lateral and medial rows.
    lateral_left = ((-420, 0, 0), None, None, None)
    lateral_right = ((420, 0, 0), None, None, None)
    specs = [
        (left, renderer.left_mask, lateral_left, 100, 0.78),    # left lateral
        (right, renderer.right_mask, lateral_right, 100, 0.78),  # right lateral
        (left, renderer.left_mask, lateral_right, -100, 0.76),   # left medial
        (right, renderer.right_mask, lateral_left, -100, 0.76),  # right medial
    ]
    crops = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        for index, (hemi, mask, view, separation, zoom) in enumerate(specs):
            fig_view = plt.figure(figsize=(6, 6))
            axis = fig_view.add_axes([0, 0, 1, 1])
            renderer._render_view(hemi, mask, axis, view, display_clims, cmap, separation=separation, zoom=zoom,
                                  volumetric_structures="subcortex" if subcortex else None)
            view_path = Path(tmp_dir) / f"view_{index}.png"
            fig_view.savefig(view_path, dpi=200, transparent=True)
            plt.close(fig_view)
            crops.append(crop_to_content(view_path))

    heights = [max(crops[0].shape[0], crops[1].shape[0]), max(crops[2].shape[0], crops[3].shape[0])]
    widths = [max(crops[0].shape[1], crops[2].shape[1]), max(crops[1].shape[1], crops[3].shape[1])]
    total_h, total_w = sum(heights), sum(widths)
    fig_w = 10.0
    path.parent.mkdir(parents=True, exist_ok=True)

    np.save(path.with_name(path.stem + "_crops.npy"), np.array(crops, dtype=object), allow_pickle=True)
    for with_legend, out_path in ((True, path), (False, path.with_name(path.stem + "_nolegend.png"))):
        extra = 1.9 if with_legend else 0.0  # room for colour bar and legend
        fig_h = fig_w * total_h / total_w + extra
        fig = plt.figure(figsize=(fig_w, fig_h))
        grid_bottom = (1.75 / fig_h) if with_legend else 0.01
        grid = fig.add_gridspec(2, 2, left=0.01, right=0.99, top=0.99, bottom=grid_bottom,
                                height_ratios=heights, width_ratios=widths, hspace=0.01, wspace=0.01)
        for index, crop in enumerate(crops):
            axis = fig.add_subplot(grid[index // 2, index % 2])
            axis.imshow(crop)
            axis.axis("off")
        if with_legend:
            scalar = plt.cm.ScalarMappable(norm=Normalize(*rho_lims), cmap=plt.get_cmap(CMAP_NAME))
            scalar.set_array([])
            cax = fig.add_axes([0.30, 1.12 / fig_h, 0.40, 0.2 / fig_h])
            colorbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
            colorbar.set_label("Spearman r (significant sites)", fontsize=13, labelpad=4)
            colorbar.ax.tick_params(labelsize=12, width=1.2, length=5)
            colorbar.outline.set_linewidth(1.2)
            fig.legend(handles=legend_handles(counts, min_targets, alpha), loc="lower center",
                       bbox_to_anchor=(0.5, 0.02 / fig_h), ncol=3, frameon=False, fontsize=11,
                       handlelength=1.4, columnspacing=1.4)
        if title:
            fig.suptitle(title, fontsize=14, y=0.995)
        fig.savefig(out_path, dpi=200, facecolor="white")
        plt.close(fig)


def legend_handles(counts: dict[str, int], min_targets: int, alpha: float) -> list[Patch]:
    return [
        Patch(facecolor=UNTESTED_COLOR, edgecolor="none",
              label=f"not tested, < {min_targets} targets (n = {counts['not tested']})"),
        Patch(facecolor=NS_COLOR, edgecolor="none",
              label=f"n.s. (n = {counts['n.s.']})"),
        Patch(facecolor=plt.get_cmap(CMAP_NAME)(0.6), edgecolor="none",
              label=f"FDR q < {alpha:g} (n = {counts['significant']})"),
    ]


def save_colorbar_legend(
    rho_lims: tuple[float, float],
    counts: dict[str, int],
    min_targets: int,
    alpha: float,
    path: Path,
    font_size: float,
) -> None:
    """Standalone, centre-aligned colour bar plus the three category swatches."""
    fig = plt.figure(figsize=(12.0, 3.4))
    scalar = plt.cm.ScalarMappable(norm=Normalize(*rho_lims), cmap=plt.get_cmap(CMAP_NAME))
    scalar.set_array([])
    cax = fig.add_axes([0.20, 0.60, 0.60, 0.13])
    colorbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
    colorbar.set_label("Spearman r (significant sites)", fontsize=font_size + 2, labelpad=6)
    colorbar.ax.tick_params(labelsize=font_size, width=1.6, length=7)
    colorbar.outline.set_linewidth(1.6)
    fig.legend(handles=legend_handles(counts, min_targets, alpha), loc="lower center",
               bbox_to_anchor=(0.5, 0.02), ncol=3, frameon=False, fontsize=font_size,
               handlelength=1.6, handleheight=1.1, columnspacing=1.8)
    fig.savefig(path, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)


def crop_to_content(path: Path, pad: int = 6) -> np.ndarray:
    """Load a view PNG and crop it to the bounding box of its non-white content.

    Cerebro draws an opaque white-background bitmap into the axes, so the
    alpha channel carries no information; content is any pixel that is not
    (near-)white after compositing onto white.
    """
    image = plt.imread(path)
    if image.shape[-1] == 4:
        rgb = image[..., :3] * image[..., 3:4] + (1.0 - image[..., 3:4])
    else:
        rgb = image[..., :3]
    content = rgb.min(axis=2) < 0.97
    rows = np.flatnonzero(content.any(axis=1))
    cols = np.flatnonzero(content.any(axis=0))
    if rows.size == 0 or cols.size == 0:
        return np.clip(rgb, 0.0, 1.0)
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, rgb.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, rgb.shape[1])
    return np.clip(rgb[r0:r1, c0:c1], 0.0, 1.0)


def plot_distribution(
    rho: np.ndarray,
    count: np.ndarray,
    category: np.ndarray,
    min_targets: int,
    alpha: float,
    path_two_panel: Path,
    path_hist_only: Path,
) -> None:
    tested = np.isfinite(rho)
    sig = category == "significant"
    ns = tested & ~sig
    n_tested, n_sig = int(tested.sum()), int(sig.sum())
    median = float(np.nanmedian(rho[tested]))
    sig_color = plt.get_cmap(CMAP_NAME)(0.6)
    lo = np.floor(np.nanmin(rho[tested]) * 10) / 10 - 0.05
    bins = np.arange(min(lo, 0.0), 1.0001, 0.05)

    def draw_hist(ax: plt.Axes) -> None:
        ax.hist([rho[sig], rho[ns]], bins=bins, stacked=True, color=[sig_color, NS_COLOR],
                edgecolor="white", linewidth=0.6,
                label=[f"FDR q < {alpha:g} (n = {n_sig})", f"n.s. (n = {n_tested - n_sig})"])
        ax.axvline(median, color="black", linestyle="--", linewidth=1.2)
        ax.text(median, ax.get_ylim()[1] * 0.97, f" median r = {median:.2f}", va="top", ha="left", fontsize=11)
        ax.set_xlabel("Spearman r per stimulation site", fontsize=13)
        ax.set_ylabel("Number of stimulation sites", fontsize=13)
        ax.legend(frameon=False, fontsize=11, loc="upper left")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=11)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    draw_hist(axes[0])
    ax = axes[1]
    ax.scatter(count[ns], rho[ns], s=22, color=NS_COLOR, label="n.s.", zorder=2)
    ax.scatter(count[sig], rho[sig], s=22, color=sig_color, label=f"FDR q < {alpha:g}", zorder=3)
    ax.axhline(0.0, color="grey", linewidth=0.8)
    ax.set_xlabel(f"Eligible targets per site (≥ {min_targets})", fontsize=13)
    ax.set_ylabel("Spearman r", fontsize=13)
    ax.set_title("Site-wise r vs. number of eligible targets", fontsize=13)
    ax.legend(frameon=False, fontsize=11, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(path_two_panel, dpi=200, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    draw_hist(ax)
    fig.tight_layout()
    fig.savefig(path_hist_only, dpi=200, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--rc-npz", type=Path, default=DEFAULT_RC_NPZ)
    parser.add_argument("--outcome", choices=("Probability", "Amplitude"), default="Probability")
    parser.add_argument("--min-targets", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--rho-lims", type=float, nargs=2, default=None,
                        help="colour range for significant rho (default: data range rounded to 0.1)")
    parser.add_argument("--title", type=str, default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-render", action="store_true", help="statistics and distribution only")
    parser.add_argument("--legend-font-size", type=float, default=20.0,
                        help="font size of the standalone colour bar / legend figure")
    parser.add_argument("--no-subcortex", action="store_true",
                        help="render the cortical surfaces only (no grey subcortical structures)")
    parser.add_argument("--map-name", type=str, default="fig3c_site_rho_map",
                        help="file stem of the map outputs")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = rcc.load_comparison_data(args.data_root)
    order_sc = rcc.load_permutation(args.data_root / "ReOrder_Caio_SC.mat", "ReOrder_Caio_SC", 374)
    with np.load(args.rc_npz) as npz:
        peak = npz["peak"]
    rc = np.round(peak[np.ix_(order_sc, order_sc)][:360, :360], 11)
    emp = data.probability if args.outcome == "Probability" else data.amplitude
    mask = data.masks[(args.outcome, "All pairs")]

    rho, p, count = per_site_spearman(rc, emp, mask, args.min_targets)
    q = bh_fdr(p)
    category = categorise(rho, q, args.alpha)
    counts = {key: int((category == key).sum())
              for key in ("not tested", "n.s.", "significant", "significant negative")}
    tested = np.isfinite(rho)
    summary = {
        "sites_total": 360,
        "sites_tested": int(tested.sum()),
        "sites_significant_positive": counts["significant"],
        "sites_significant_negative": counts["significant negative"],
        "sites_ns": counts["n.s."],
        "sites_not_tested": counts["not tested"],
        "sites_p_lt_alpha_uncorrected_positive": int((tested & (p < args.alpha) & (rho > 0)).sum()),
        "rho_median_tested": float(np.nanmedian(rho[tested])),
        "rho_iqr_tested": [float(np.nanpercentile(rho[tested], 25)), float(np.nanpercentile(rho[tested], 75))],
        "rho_min_tested": float(np.nanmin(rho[tested])),
        "rho_max_tested": float(np.nanmax(rho[tested])),
        "targets_per_tested_site_median": float(np.median(count[tested])),
        "targets_per_tested_site_min": int(count[tested].min()),
    }
    print(json.dumps(summary, indent=2))

    names = data.node_names[:360]
    with (args.output_dir / "site_correlations.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["site_index_1based", "node_name", "hemisphere", "n_eligible_targets",
                         "spearman_rho", "p_asymptotic", "q_bh_fdr", "category"])
        for i in range(360):
            writer.writerow([i + 1, names[i], "L" if i < 180 else "R", count[i],
                             "" if not np.isfinite(rho[i]) else f"{rho[i]:.6f}",
                             "" if not np.isfinite(p[i]) else f"{p[i]:.3e}",
                             "" if not np.isfinite(q[i]) else f"{q[i]:.3e}",
                             category[i]])

    plot_distribution(rho, count, category, args.min_targets, args.alpha,
                      args.output_dir / "fig3d_site_rho_distribution.png",
                      args.output_dir / "fig3d_site_rho_histogram.png")

    sig = category == "significant"
    if args.rho_lims is not None:
        rho_lims = (float(args.rho_lims[0]), float(args.rho_lims[1]))
    else:
        rho_lims = (float(np.floor(rho[sig].min() * 10) / 10), float(np.ceil(rho[sig].max() * 10) / 10))
    map_path = args.output_dir / f"{args.map_name}.png"
    save_colorbar_legend(rho_lims, counts, args.min_targets, args.alpha,
                         args.output_dir / "fig3c_colorbar_legend.png", args.legend_font_size)
    if not args.skip_render:
        from rc_ftract_cerebro import CerebroRenderer, Config, load_rc_data

        config = Config(model_variant="main_caio", dense_mapping_mode="matlab_sc_order")
        rc_data = load_rc_data(config)
        if rc_data.node_names[:360] != names:
            raise ValueError("Renderer node order differs from the comparison-data node order.")
        renderer = CerebroRenderer(config.atlas_file, config.cortical_atlas_file,
                                   rc_data.node_names, rc_data.raw_order, config.dense_mapping_mode)
        cmap, display_clims, sentinel = build_colormap(*rho_lims)
        values = np.full(374, np.nan)
        values[:360][category == "n.s."] = sentinel
        values[:360][category == "significant negative"] = sentinel
        values[:360][sig] = np.clip(rho[sig], rho_lims[0], rho_lims[1])
        render_map(renderer, values, cmap, display_clims, rho_lims, map_path, counts,
                   args.min_targets, args.alpha, args.title, subcortex=not args.no_subcortex)
        print(f"Wrote {map_path}")

    provenance = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": "code/python_figures/render_fig3c_site_correlation.py",
        "rc_npz": str(args.rc_npz),
        "rc_npz_sha256": sha256(args.rc_npz),
        "rc_orientation": "rows = stimulation site (un-transposed, as in Fig. 3A), cortical 360 x 360 in ReOrder_Caio_SC order",
        "empirical": f"F-tract-2 / MNI-HCP-MMP1 {args.outcome.lower()} (median amplitude if Amplitude)",
        "pair_mask": "Fig. 3A all-pairs mask: same hemisphere, off-diagonal, ED >= 20 mm, N_implantations >= 5, "
                     + ("N_stimulations >= 100, finite probability" if args.outcome == "Probability"
                        else "N_stimulations x probability >= 100, finite amplitude"),
        "min_targets_per_site": args.min_targets,
        "test": "Spearman per site, asymptotic p (scipy.stats.spearmanr), Benjamini-Hochberg across tested sites",
        "alpha": args.alpha,
        "rho_colour_limits": list(rho_lims),
        "summary": summary,
        "renderer": None if args.skip_render else
                    "rc_ftract_cerebro.CerebroRenderer (main_caio node order, matlab_sc_order dense mapping, "
                    "Cerebro inflated surfaces + subcortex in grey), colormap plasma with n.s. sentinel",
    }
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    readme = f"""# Per-site correlation map (Fig. 3C) and distribution (Fig. 3D)

Generated {provenance['generated']} by `{provenance['script']}`.

* RC matrix: `{args.rc_npz.relative_to(PROJECT_ROOT) if args.rc_npz.is_relative_to(PROJECT_ROOT) else args.rc_npz}` (X_struct.v8.W run behind Fig. 3A), rows = stimulation site.
* Empirical: F-TRACT {args.outcome.lower()}, Fig. 3A all-pairs inclusion mask applied per site.
* Site tested when >= {args.min_targets} eligible targets; Spearman rho, BH-FDR across tested sites, alpha = {args.alpha:g}.

| quantity | value |
|---|---|
| sites tested | {summary['sites_tested']} / 360 |
| significant (FDR q < {args.alpha:g}, rho > 0) | {summary['sites_significant_positive']} |
| significant negative | {summary['sites_significant_negative']} |
| n.s. | {summary['sites_ns']} |
| not tested (< {args.min_targets} targets) | {summary['sites_not_tested']} |
| p < {args.alpha:g} uncorrected, rho > 0 | {summary['sites_p_lt_alpha_uncorrected_positive']} |
| rho median (IQR), tested sites | {summary['rho_median_tested']:.3f} ({summary['rho_iqr_tested'][0]:.3f}-{summary['rho_iqr_tested'][1]:.3f}) |
| rho range, tested sites | {summary['rho_min_tested']:.3f} to {summary['rho_max_tested']:.3f} |
| eligible targets per tested site, median (min) | {summary['targets_per_tested_site_median']:.0f} ({summary['targets_per_tested_site_min']}) |

Files: `fig3c_site_rho_map.png` (light grey = not tested, dark grey = n.s., colour = rho in
[{rho_lims[0]:.1f}, {rho_lims[1]:.1f}]), `fig3c_site_rho_map_nolegend.png` (same map without colour bar / legend),
`fig3c_colorbar_legend.png` (standalone colour bar + category swatches, font size {args.legend_font_size:g}),
`fig3d_site_rho_distribution.png`, `fig3d_site_rho_histogram.png`, `site_correlations.csv`, `provenance.json`.
"""
    (args.output_dir / "README.md").write_text(readme)


if __name__ == "__main__":
    main()
