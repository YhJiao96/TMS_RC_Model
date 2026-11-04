#!/usr/bin/env python3
"""Render static DLPFC peak-current and empirical-probability brain maps."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from rc_ftract_cerebro import CerebroRenderer, Config, load_rc_data

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter

from render_three_currents_four_timepoints import (
    combine_cerebro_row,
    horizontally_crop_cerebro_image,
    pad_rows_to_common_canvas,
)


def render_compact_static_map(
    renderer: CerebroRenderer,
    values: np.ndarray,
    output_path: Path,
    *,
    title: str,
    clims: tuple[float, float],
    colorbar_label: str,
    scientific_colorbar: bool,
) -> None:
    """Render a static map with the same compact geometry as single frames."""
    font_scale = 1.5
    font_weight = "bold"
    cmap = plt.get_cmap("magma")
    dense = renderer.to_dense(values)
    left = dense.copy()
    left[~renderer.left_mask] = np.nan
    right = dense.copy()
    right[~renderer.right_mask] = np.nan

    figure = plt.figure(figsize=(7.2, 9.5), facecolor="white")
    grid = figure.add_gridspec(
        3,
        2,
        left=0.015,
        right=0.985,
        bottom=0.14,
        top=0.82,
        height_ratios=(1.0, 1.0, 0.24),
        hspace=0.12,
        wspace=0.0,
    )
    axes = np.asarray(
        [
            [figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])],
            [figure.add_subplot(grid[1, 0]), figure.add_subplot(grid[1, 1])],
        ]
    )
    renderer._render_view(
        left,
        renderer.left_mask,
        axes[0, 0],
        ((-420, 0, 0), None, None, None),
        clims,
        cmap,
        separation=100,
        zoom=0.90,
    )
    renderer._render_view(
        right,
        renderer.right_mask,
        axes[0, 1],
        ((420, 0, 0), None, None, None),
        clims,
        cmap,
        separation=100,
        zoom=0.90,
    )
    renderer._render_view(
        left,
        renderer.left_mask,
        axes[1, 0],
        ((420, 0, 0), None, None, None),
        clims,
        cmap,
        separation=-100,
        zoom=0.76,
    )
    renderer._render_view(
        right,
        renderer.right_mask,
        axes[1, 1],
        ((-420, 0, 0), None, None, None),
        clims,
        cmap,
        separation=-100,
        zoom=0.76,
    )

    row_images = []
    for row_index in range(2):
        row_images.append(
            combine_cerebro_row(
                horizontally_crop_cerebro_image(axes[row_index, 0]),
                horizontally_crop_cerebro_image(axes[row_index, 1]),
            )
        )
    row_images = pad_rows_to_common_canvas(row_images)
    for brain_axis in axes.flat:
        brain_axis.remove()
    for row_index, row_image in enumerate(row_images):
        combined_axis = figure.add_subplot(grid[row_index, :])
        combined_axis.imshow(row_image, origin="lower", interpolation="none")
        combined_axis.set_aspect("equal", adjustable="box")
        combined_axis.axis("off")

    colorbar_axis = figure.add_subplot(grid[2, :])
    colorbar_position = colorbar_axis.get_position()
    colorbar_axis.set_position(
        [0.06, colorbar_position.y0, 0.88, colorbar_position.height]
    )
    scalar = plt.cm.ScalarMappable(norm=Normalize(*clims), cmap=cmap)
    scalar.set_array([])
    colorbar = figure.colorbar(
        scalar,
        cax=colorbar_axis,
        orientation="horizontal",
    )
    formatter = ScalarFormatter(useMathText=scientific_colorbar)
    formatter.set_scientific(scientific_colorbar)
    if scientific_colorbar:
        formatter.set_powerlimits((0, 0))
    else:
        formatter.set_useOffset(False)
    colorbar.formatter = formatter
    colorbar.locator.set_params(nbins=4)
    colorbar.update_ticks()
    colorbar.ax.tick_params(labelsize=15 * font_scale, width=1.6, length=7)
    for tick_label in colorbar.ax.get_xticklabels():
        tick_label.set_fontweight(font_weight)
    offset_text = colorbar.ax.xaxis.get_offset_text()
    offset_text.set_fontsize(14 * font_scale)
    offset_text.set_fontweight(font_weight)
    colorbar.set_label(
        colorbar_label,
        fontsize=17 * font_scale,
        fontweight=font_weight,
        labelpad=4,
    )
    colorbar.outline.set_linewidth(1.6)
    figure.suptitle(
        title,
        fontsize=16 * font_scale,
        fontweight=font_weight,
        y=0.98,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, facecolor="white")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--atlas-file", type=Path, required=True)
    parser.add_argument("--cortical-atlas-file", type=Path, required=True)
    parser.add_argument("--dynamic-file", type=Path, required=True)
    parser.add_argument("--empirical-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--output-prefix",
        default="DLPFC",
        help="Prefix of the output PNG/NPZ file names (default: DLPFC).",
    )
    args = parser.parse_args()

    config = Config(
        data_root=args.data_root,
        atlas_file=args.atlas_file,
        cortical_atlas_file=args.cortical_atlas_file,
        model_variant="legacy_xstruct_v8",
        output_dir=args.output_dir,
        dense_mapping_mode="matlab_sc_order",
    )
    data = load_rc_data(config)
    renderer = CerebroRenderer(
        config.atlas_file,
        config.cortical_atlas_file,
        data.node_names,
        data.raw_order,
        config.dense_mapping_mode,
    )

    with np.load(args.dynamic_file) as dynamic:
        dynamic_names = dynamic["node_names"].astype(str).tolist()
        if dynamic_names != data.node_names:
            raise ValueError("Dynamic-result node order does not match X_struct.v8.")
        stim_indices = np.asarray(
            dynamic["stim_indices_zero_based"], dtype=int
        ).reshape(-1)
        if stim_indices.size != 1:
            raise ValueError("This static comparison expects one stimulation node.")
        stim_index = int(stim_indices[0])
        stim_label = str(dynamic["stim_labels"].reshape(-1)[0])
        if data.node_names[stim_index] != stim_label:
            raise ValueError("Stimulus index and label disagree.")
        peak_current_unmasked = np.nanmax(dynamic["current_abs_total"], axis=1)
        source_time_ms = np.asarray(dynamic["time_ms"], dtype=float)

    with np.load(args.empirical_file) as empirical:
        empirical_names = empirical["parcel_list"].astype(str).tolist()
        if empirical_names != data.node_names[: data.n_cortex]:
            raise ValueError(
                "Empirical parcel order does not match the cortical X_struct order."
            )
        probability_matrix = np.asarray(empirical["probability"], dtype=float)
        if probability_matrix.shape != (data.n_cortex, data.n_cortex):
            raise ValueError("Empirical probability must be 360 x 360.")

    # I_node_max is [receiver, stimulation], so use P[:, stim] for the
    # receiver-wise empirical map paired with this DLPFC stimulation.
    excluded_by_distance = (
        np.isfinite(data.D[:, stim_index]) & (data.D[:, stim_index] < 20.0)
    )
    excluded_by_distance[stim_index] = True
    peak_current = peak_current_unmasked.copy()
    peak_current[excluded_by_distance] = np.nan
    peak_current_log10 = np.full_like(peak_current, np.nan)
    positive_peak = np.isfinite(peak_current) & (peak_current > 0)
    peak_current_log10[positive_peak] = np.log10(peak_current[positive_peak])
    empirical_probability = np.full(len(data.node_names), np.nan, dtype=float)
    empirical_probability[: data.n_cortex] = probability_matrix[:, stim_index]
    empirical_probability[excluded_by_distance] = np.nan

    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix
    peak_path = args.output_dir / f"{prefix}_peak_total_absolute_current.png"
    peak_log10_path = (
        args.output_dir / f"{prefix}_peak_log10_total_absolute_current.png"
    )
    probability_path = args.output_dir / f"{prefix}_empirical_probability.png"
    render_compact_static_map(
        renderer,
        peak_current,
        peak_path,
        title="Peak total absolute current",
        clims=(0.0, max(float(np.nanmax(peak_current)), np.finfo(float).eps)),
        colorbar_label="Peak model current (a.u.)",
        scientific_colorbar=True,
    )
    render_compact_static_map(
        renderer,
        peak_current_log10,
        peak_log10_path,
        title="Peak total absolute current",
        clims=(
            float(np.nanmin(peak_current_log10)),
            float(np.nanmax(peak_current_log10)),
        ),
        colorbar_label="log10 peak model current (a.u.)",
        scientific_colorbar=False,
    )
    render_compact_static_map(
        renderer,
        empirical_probability,
        probability_path,
        title="Empirical response probability",
        clims=(0.0, 1.0),
        colorbar_label="Empirical probability",
        scientific_colorbar=False,
    )

    np.savez_compressed(
        args.output_dir / f"{prefix}_peak_current_empirical_probability.npz",
        peak_total_absolute_current=peak_current,
        peak_total_absolute_current_unmasked=peak_current_unmasked,
        peak_total_absolute_current_log10=peak_current_log10,
        peak_log_base=10.0,
        peak_log_nonpositive_is_nan=True,
        empirical_probability=empirical_probability,
        stim_index_zero_based=stim_index,
        stim_index_matlab=stim_index + 1,
        stim_label=stim_label,
        node_names=np.asarray(data.node_names),
        n_cortex=data.n_cortex,
        subcortex_empirical_is_nan=True,
        empirical_orientation="probability[receiver, stimulation]",
        stimulated_node_masked=True,
        distance_exclusion_threshold_mm=20.0,
        excluded_by_distance_mask=excluded_by_distance,
        excluded_node_indices_zero_based=np.flatnonzero(excluded_by_distance),
        excluded_node_indices_matlab=np.flatnonzero(excluded_by_distance) + 1,
        excluded_node_names=np.asarray(data.node_names)[excluded_by_distance],
        distance_to_stim_mm=data.D[:, stim_index],
        source_time_ms=source_time_ms,
        capacitance_scale_factor=10.0,
        current_scale_factor=1.0,
    )
    print(f"Saved {peak_path}", flush=True)
    print(f"Saved {peak_log10_path}", flush=True)
    print(f"Saved {probability_path}", flush=True)


if __name__ == "__main__":
    main()
