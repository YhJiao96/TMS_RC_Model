#!/usr/bin/env python3
"""Render Figure-C-style multi-timepoint Cerebro panels for three currents."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rc_model_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import ScalarFormatter
from PIL import Image

from rc_ftract_cerebro import (
    CerebroRenderer,
    Config,
    load_rc_data,
    nearest_time_indices,
)


METRICS = (
    (
        "total_abs",
        "Total absolute model current",
        "current_abs_total",
    ),
    (
        "positive_inflow",
        "Positive inflow model current",
        "current_positive_inflow",
    ),
    (
        "net_positive",
        "Positive-only net model current",
        "current_net",
    ),
)

# Crop only transparent margins from the shared row canvas so that every brain
# is enlarged isotropically without changing the Cerebro camera or anatomy.
BRAIN_SIZE_SCALE = 1.20
EXTERNAL_VIEW_ZOOM = 0.90
INTERNAL_VIEW_ZOOM = 0.76


def horizontally_crop_cerebro_image(
    axis: plt.Axes,
    *,
    white_threshold: int = 250,
    padding_fraction: float = 0.02,
) -> np.ndarray:
    """Return RGBA pixels with only the empty left/right columns removed."""
    if not axis.images:
        raise RuntimeError("Cerebro did not add an image to the target axes.")
    image = axis.images[-1]
    pixels = np.asarray(image.get_array())
    if pixels.ndim != 3 or pixels.shape[2] < 3:
        raise ValueError("Expected a Cerebro RGB/RGBA offscreen image.")

    foreground = np.any(pixels[:, :, :3] < white_threshold, axis=2)
    _, columns = np.where(foreground)
    if columns.size == 0:
        raise ValueError("Cerebro image contains no visible anatomy.")

    column_pad = max(2, int(round(padding_fraction * pixels.shape[1])))
    column_start = max(0, int(columns.min()) - column_pad)
    column_stop = min(pixels.shape[1], int(columns.max()) + column_pad + 1)
    cropped_rgb = np.asarray(
        pixels[:, column_start:column_stop, :3],
        dtype=np.uint8,
    )
    alpha = np.where(
        np.all(cropped_rgb >= white_threshold, axis=2),
        0,
        255,
    ).astype(np.uint8)
    return np.dstack((cropped_rgb, alpha))


def tightly_crop_cerebro_image(
    axis: plt.Axes,
    *,
    white_threshold: int = 250,
    padding_pixels: int = 2,
) -> np.ndarray:
    """Crop all four empty borders from one independently rendered view."""
    if not axis.images:
        raise RuntimeError("Cerebro did not add an image to the target axes.")
    pixels = np.asarray(axis.images[-1].get_array())
    if pixels.ndim != 3 or pixels.shape[2] < 3:
        raise ValueError("Expected a Cerebro RGB/RGBA offscreen image.")

    rgb = np.asarray(pixels[:, :, :3], dtype=np.uint8)
    foreground = np.any(rgb < white_threshold, axis=2)
    if pixels.shape[2] >= 4:
        foreground &= np.asarray(pixels[:, :, 3]) > 0
    rows, columns = np.where(foreground)
    if rows.size == 0:
        raise ValueError("Cerebro image contains no visible anatomy.")

    row_start = max(0, int(rows.min()) - padding_pixels)
    row_stop = min(rgb.shape[0], int(rows.max()) + padding_pixels + 1)
    column_start = max(0, int(columns.min()) - padding_pixels)
    column_stop = min(rgb.shape[1], int(columns.max()) + padding_pixels + 1)
    cropped_rgb = rgb[row_start:row_stop, column_start:column_stop]
    alpha = np.where(
        np.all(cropped_rgb >= white_threshold, axis=2),
        0,
        255,
    ).astype(np.uint8)
    return np.dstack((cropped_rgb, alpha))


def resize_rgba_image(image: np.ndarray, scale: float) -> np.ndarray:
    """Resize one tightly cropped view isotropically."""
    if scale <= 0:
        raise ValueError("Image scale must be positive.")
    target_size = (
        max(1, int(round(image.shape[1] * scale))),
        max(1, int(round(image.shape[0] * scale))),
    )
    resized = Image.fromarray(image, mode="RGBA").resize(
        target_size,
        resample=Image.Resampling.LANCZOS,
    )
    return np.asarray(resized)


def combine_cerebro_row(
    left_image: np.ndarray,
    right_image: np.ndarray,
    *,
    gap_pixels: int = 6,
) -> np.ndarray:
    """Join two views without rescaling either image."""
    height = max(left_image.shape[0], right_image.shape[0])
    width = left_image.shape[1] + gap_pixels + right_image.shape[1]
    combined = np.zeros((height, width, 4), dtype=np.uint8)
    left_y = (height - left_image.shape[0]) // 2
    right_y = (height - right_image.shape[0]) // 2
    combined[
        left_y : left_y + left_image.shape[0],
        : left_image.shape[1],
    ] = left_image
    right_x = left_image.shape[1] + gap_pixels
    combined[
        right_y : right_y + right_image.shape[0],
        right_x : right_x + right_image.shape[1],
    ] = right_image
    return combined


def pad_rows_to_common_canvas(row_images: list[np.ndarray]) -> list[np.ndarray]:
    """Center rows on identical canvases so their anatomical scale matches."""
    common_height = max(image.shape[0] for image in row_images)
    common_width = max(image.shape[1] for image in row_images)
    padded_rows: list[np.ndarray] = []
    for image in row_images:
        canvas = np.zeros((common_height, common_width, 4), dtype=np.uint8)
        y_start = (common_height - image.shape[0]) // 2
        x_start = (common_width - image.shape[1]) // 2
        canvas[
            y_start : y_start + image.shape[0],
            x_start : x_start + image.shape[1],
        ] = image
        padded_rows.append(canvas)
    return padded_rows


def render_time_evolution(
    renderer: CerebroRenderer,
    display_time_ms: np.ndarray,
    values: np.ndarray,
    metric_slug: str,
    display_name: str,
    output_path: Path,
    model_label: str = "",
    font_scale: float = 1.0,
    font_weight: str = "normal",
    log10_transform: bool = False,
    site_label: str = "DLPFC",
) -> None:
    """Render multiple 2x2 Cerebro view groups in one compact landscape figure."""
    if font_scale <= 0:
        raise ValueError("font_scale must be positive.")
    cmap = plt.get_cmap("magma")
    panel_count = len(display_time_ms)
    enlarged_text = font_scale > 1.0
    figure_height = 9.5 if enlarged_text else 8.5
    figure_width = 7.2 * panel_count
    figure = plt.figure(
        figsize=(figure_width, figure_height),
        facecolor="white",
    )
    outer = figure.add_gridspec(
        1,
        panel_count,
        left=0.015,
        right=0.985,
        bottom=0.14 if enlarged_text else 0.105,
        top=0.82 if enlarged_text else 0.885,
        wspace=0.035,
    )

    for panel_index, shown_time in enumerate(display_time_ms):
        frame_values = values[:, panel_index]
        if metric_slug == "net_positive":
            frame_values = np.where(frame_values > 0.0, frame_values, np.nan)

        if log10_transform:
            positive = np.isfinite(frame_values) & (frame_values > 0.0)
            logged_values = np.full_like(frame_values, np.nan, dtype=float)
            logged_values[positive] = np.log10(frame_values[positive])
            frame_values = logged_values

        finite = frame_values[np.isfinite(frame_values)]
        if log10_transform:
            if finite.size:
                lower = float(np.min(finite))
                upper = float(np.max(finite))
                if np.isclose(lower, upper):
                    lower -= 0.5
                    upper += 0.5
                clims = (lower, upper)
            else:
                clims = (-1.0, 0.0)
        else:
            upper = float(np.max(finite)) if finite.size else 1.0
            upper = max(upper, np.finfo(float).eps)
            clims = (0.0, upper)

        dense = renderer.to_dense(frame_values)
        left = dense.copy()
        left[~renderer.left_mask] = np.nan
        right = dense.copy()
        right[~renderer.right_mask] = np.nan

        inner = outer[panel_index].subgridspec(
            3,
            2,
            height_ratios=(1.0, 1.0, 0.24 if enlarged_text else 0.18),
            hspace=0.12,
            wspace=0.0,
        )
        axes = np.asarray(
            [
                [figure.add_subplot(inner[0, 0]), figure.add_subplot(inner[0, 1])],
                [figure.add_subplot(inner[1, 0]), figure.add_subplot(inner[1, 1])],
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
            zoom=EXTERNAL_VIEW_ZOOM,
        )
        renderer._render_view(
            right,
            renderer.right_mask,
            axes[0, 1],
            ((420, 0, 0), None, None, None),
            clims,
            cmap,
            separation=100,
            zoom=EXTERNAL_VIEW_ZOOM,
        )
        renderer._render_view(
            left,
            renderer.left_mask,
            axes[1, 0],
            ((420, 0, 0), None, None, None),
            clims,
            cmap,
            separation=-100,
            zoom=INTERNAL_VIEW_ZOOM,
        )
        renderer._render_view(
            right,
            renderer.right_mask,
            axes[1, 1],
            ((-420, 0, 0), None, None, None),
            clims,
            cmap,
            separation=-100,
            zoom=INTERNAL_VIEW_ZOOM,
        )

        row_images = []
        for row_index in range(2):
            left_image = resize_rgba_image(
                tightly_crop_cerebro_image(axes[row_index, 0]),
                BRAIN_SIZE_SCALE,
            )
            right_image = resize_rgba_image(
                tightly_crop_cerebro_image(axes[row_index, 1]),
                BRAIN_SIZE_SCALE,
            )
            row_images.append(
                combine_cerebro_row(
                    left_image,
                    right_image,
                )
            )
        row_images = pad_rows_to_common_canvas(row_images)
        for brain_axis in axes.flat:
            brain_axis.remove()
        for row_index, row_image in enumerate(row_images):
            combined_axis = figure.add_subplot(inner[row_index, :])
            combined_axis.imshow(row_image, origin="lower", interpolation="none")
            combined_axis.set_aspect("equal", adjustable="box")
            combined_axis.axis("off")

        colorbar_axis = figure.add_subplot(inner[2, :])
        scalar = plt.cm.ScalarMappable(norm=Normalize(*clims), cmap=cmap)
        scalar.set_array([])
        colorbar = figure.colorbar(
            scalar,
            cax=colorbar_axis,
            orientation="horizontal",
        )
        formatter = ScalarFormatter(useMathText=not log10_transform)
        formatter.set_scientific(not log10_transform)
        if log10_transform:
            formatter.set_useOffset(False)
        else:
            formatter.set_powerlimits((0, 0))
        colorbar.formatter = formatter
        colorbar.locator.set_params(nbins=4)
        colorbar.update_ticks()
        colorbar.ax.tick_params(
            labelsize=15 * font_scale,
            width=1.6,
            length=7,
        )
        for tick_label in colorbar.ax.get_xticklabels():
            tick_label.set_fontweight(font_weight)
        offset_text = colorbar.ax.xaxis.get_offset_text()
        offset_text.set_fontsize(14 * font_scale)
        offset_text.set_fontweight(font_weight)
        colorbar.set_label(
            (
                "log10 model current (a.u.)"
                if log10_transform
                else "Model current (a.u.)"
            ),
            fontsize=17 * font_scale,
            fontweight=font_weight,
            labelpad=4 if enlarged_text else 3,
        )
        colorbar.outline.set_linewidth(1.6)

        panel_x = 0.015 + (panel_index + 0.5) * (0.97 / panel_count)
        figure.text(
            panel_x,
            0.855 if enlarged_text else 0.905,
            f"T = {shown_time:.2f} ms",
            ha="center",
            va="bottom",
            fontsize=19 * font_scale,
            fontweight=font_weight,
        )

    # Mark the left-to-right temporal progression between adjacent panels.
    arrow_y = 0.52
    for boundary_index in range(1, panel_count):
        boundary_x = 0.015 + boundary_index * (0.97 / panel_count)
        figure.add_artist(
            FancyArrowPatch(
                (boundary_x - 0.010, arrow_y),
                (boundary_x + 0.010, arrow_y),
                transform=figure.transFigure,
                arrowstyle="-|>",
                mutation_scale=16 * font_scale,
                linewidth=2.2 * font_scale,
                color="#0072B2",
                zorder=30,
                clip_on=False,
            )
        )

    # Keep the transform in the colorbar label; the main title names the
    # physical current quantity only.
    title_display_name = display_name
    if panel_count == 1:
        single_display_name = title_display_name.replace(
            " model current", " current"
        )
        title = single_display_name.capitalize()
    else:
        title = f"{site_label} stimulation — {title_display_name} evolution"
    if model_label:
        title += f" ({model_label})"
    title_base_size = 16 if panel_count == 1 else 24
    figure.suptitle(
        title,
        fontsize=title_base_size * font_scale,
        fontweight=font_weight,
        y=0.98 if enlarged_text else 0.985,
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
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--times-ms",
        type=float,
        nargs="+",
        default=(1.0, 20.0, 50.0, 100.0, 200.0),
    )
    parser.add_argument("--display-time-scale", type=float, default=1.0)
    parser.add_argument("--model-label", default="")
    parser.add_argument(
        "--site-label",
        default="DLPFC",
        help="Stimulation-site name shown in the figure title (default: DLPFC).",
    )
    parser.add_argument(
        "--output-prefix",
        default="DLPFC",
        help="Prefix of the output PNG file names (default: DLPFC).",
    )
    parser.add_argument("--font-scale", type=float, default=1.0)
    parser.add_argument("--font-weight", default="normal")
    parser.add_argument(
        "--log10",
        action="store_true",
        help="Plot log10 of positive current; non-positive values become NaN.",
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        choices=[metric[0] for metric in METRICS],
        default=[metric[0] for metric in METRICS],
        help="Current maps to render (default: all three).",
    )
    parser.add_argument(
        "--mask-stimulus",
        action="store_true",
        help="Set stimulated-node values to NaN before plotting.",
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

    with np.load(args.dynamic_file) as archive:
        requested_times = np.asarray(args.times_ms, dtype=float)
        indices = nearest_time_indices(archive["time_ms"], requested_times)
        if len(indices) != len(requested_times):
            raise ValueError(
                "Two requested times map to the same simulation frame; "
                "choose distinct timepoints."
            )
        actual_times = archive["time_ms"][indices]
        display_times = requested_times * args.display_time_scale
        stim_indices = np.asarray(
            archive["stim_indices_zero_based"], dtype=int
        ).reshape(-1)
        for metric_slug, display_name, archive_key in METRICS:
            if metric_slug not in args.metrics:
                continue
            metric_values = np.asarray(
                archive[archive_key][:, indices], dtype=float
            ).copy()
            if args.mask_stimulus:
                metric_values[stim_indices, :] = np.nan
            transform_slug = "_log10" if args.log10 else ""
            output_path = args.output_dir / (
                f"{args.output_prefix}_{metric_slug}{transform_slug}_"
                f"{len(actual_times)}_timepoints.png"
            )
            render_time_evolution(
                renderer,
                display_times,
                metric_values,
                metric_slug,
                display_name,
                output_path,
                model_label=args.model_label,
                font_scale=args.font_scale,
                font_weight=args.font_weight,
                log10_transform=args.log10,
                site_label=args.site_label,
            )
            print(f"Saved {output_path}", flush=True)


if __name__ == "__main__":
    main()
