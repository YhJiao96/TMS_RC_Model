#!/usr/bin/env python3
"""RC/F-TRACT model in Python with Tian S1 Cerebro visualisation.

This is a Python counterpart to ``Matlab/RC_ftract_Caio_clean.m``.  It keeps
the current scientific conventions used by the MATLAB pipeline:

* SC, distance (D), and connection frequency (F) are read from
  ``SC_data_neuron_2023_yihang.mat``.
* ``ReOrder_Caio_SC`` is applied as a parcel permutation.
* Edges with F < 0.25 are removed.
* R = D / SC, C = 1, and stimulation is an initial voltage pulse.
* The dynamic node current used for peak-current analysis is
  ``sum_j(abs((V_i - V_j) / R_ij))``.
* ``N_imp = N_stimulations * probability`` estimates N_resp for amplitude
  thresholding.

The default command simulates a left-DLPFC pulse at ``L_p9-46v_ROI``, saves
the complete 374-node voltage/current time course, creates a dynamic summary
plot, and renders selected frames with Cerebro using the Glasser 360 + Tian
S1 atlas.  The model has 14 subcortical nodes.  Before rendering, every model
frame is restored from F-TRACT order to the raw SC/``parcel_list`` order.  The
default MATLAB-compatible mapping then expands raw SC rows 1:360 through the
Glasser-360 cortical dlabel and rows 361:374 through the corresponding 14
standard CIFTI subcortical structures into 91,282 grayordinates.

Inputs are resolved through ``rc_ftract_core`` (environment variables
RC_DATA_ROOT, RC_FTRACT_DIR, RC_ATLAS_DIR, RC_OUTPUT_ROOT).

Examples
--------
Default DLPFC dynamic analysis::

  python rc_ftract_cerebro.py

Use another parcel and skip Cerebro rendering::

  python rc_ftract_cerebro.py --stim-label L_46_ROI --skip-cerebro

Run all stimulation nodes and F-TRACT threshold analyses::

  python rc_ftract_cerebro.py --run-all-stimulations --run-ftract
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

# Keep Matplotlib's cache outside the project and make the script usable when
# the user's home directory is read-only (for example, in an offscreen job).
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rc_model_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.ticker import ScalarFormatter
from scipy.integrate import solve_ivp
from scipy.io import loadmat
from scipy.stats import pearsonr, rankdata, spearmanr

import rc_ftract_core as rcc


DEFAULT_DATA_ROOT = rcc.DEFAULT_DATA_ROOT
DEFAULT_CEREBRO_ROOT = rcc.atlas_dir(DEFAULT_DATA_ROOT)  # folder with the .dlabel.nii atlases (RC_ATLAS_DIR)
DEFAULT_ATLAS = DEFAULT_CEREBRO_ROOT / (
    "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_"
    "Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii"
)
DEFAULT_CORTICAL_ATLAS = DEFAULT_CEREBRO_ROOT / (
    "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_"
    "Areas_Group_Colors.32k_fs_LR.dlabel.nii"
)
DEFAULT_OUTPUT = rcc.DEFAULT_OUTPUT_ROOT / "rc_ftract_cerebro"


S1_CIFTI_STRUCTURE_TRANSLATION: dict[str, str] = {
    "Left-Thalamus-Proper": "CIFTI_STRUCTURE_THALAMUS_LEFT",
    "Left-Caudate": "CIFTI_STRUCTURE_CAUDATE_LEFT",
    "Left-Putamen": "CIFTI_STRUCTURE_PUTAMEN_LEFT",
    "Left-Pallidum": "CIFTI_STRUCTURE_PALLIDUM_LEFT",
    "Left-Hippocampus": "CIFTI_STRUCTURE_HIPPOCAMPUS_LEFT",
    "Left-Amygdala": "CIFTI_STRUCTURE_AMYGDALA_LEFT",
    "Left-Accumbens-area": "CIFTI_STRUCTURE_ACCUMBENS_LEFT",
    "Right-Thalamus-Proper": "CIFTI_STRUCTURE_THALAMUS_RIGHT",
    "Right-Caudate": "CIFTI_STRUCTURE_CAUDATE_RIGHT",
    "Right-Putamen": "CIFTI_STRUCTURE_PUTAMEN_RIGHT",
    "Right-Pallidum": "CIFTI_STRUCTURE_PALLIDUM_RIGHT",
    "Right-Hippocampus": "CIFTI_STRUCTURE_HIPPOCAMPUS_RIGHT",
    "Right-Amygdala": "CIFTI_STRUCTURE_AMYGDALA_RIGHT",
    "Right-Accumbens-area": "CIFTI_STRUCTURE_ACCUMBENS_RIGHT",
}


S1_NATIVE_LABEL_TRANSLATION: dict[str, tuple[str, ...]] = {
    "Left-Thalamus-Proper": ("aTHA-lh", "pTHA-lh"),
    "Left-Caudate": ("CAU-lh",),
    "Left-Putamen": ("PUT-lh",),
    "Left-Pallidum": ("GP-lh",),
    "Left-Hippocampus": ("HIP-lh",),
    "Left-Amygdala": ("AMY-lh",),
    "Left-Accumbens-area": ("NAc-lh",),
    "Right-Thalamus-Proper": ("aTHA-rh", "pTHA-rh"),
    "Right-Caudate": ("CAU-rh",),
    "Right-Putamen": ("PUT-rh",),
    "Right-Pallidum": ("GP-rh",),
    "Right-Hippocampus": ("HIP-rh",),
    "Right-Amygdala": ("AMY-rh",),
    "Right-Accumbens-area": ("NAc-rh",),
}

VALID_CURRENT_METRICS = (
    "total_abs",
    "net",
    "net_positive",
    "positive_inflow",
)


@dataclass(frozen=True)
class Config:
    data_root: Path = DEFAULT_DATA_ROOT
    ftract_root: Path = rcc.ftract_dir(DEFAULT_DATA_ROOT).parent
    atlas_file: Path = DEFAULT_ATLAS
    cortical_atlas_file: Path = DEFAULT_CORTICAL_ATLAS
    model_variant: str = "main_caio"
    dense_mapping_mode: str = "matlab_sc_order"
    output_dir: Path = DEFAULT_OUTPUT
    parcellation: str = "MNI-HCP-MMP1"
    time_scale_ms: float = 200.0
    num_points: int = 1001
    fco_threshold: float = 0.25
    n_cortex: int = 360
    stim_labels: tuple[str, ...] = ("L_p9-46v_ROI",)
    initial_voltage: float = 1.0
    split_multisite_voltage: bool = True
    simulation_solver: str = "rk45"
    rtol: float = 1e-3
    atol: float = 1e-6
    snapshot_times_ms: tuple[float, ...] = (0.0, 10.0, 50.0, 100.0, 150.0, 200.0)
    animation_step_ms: float = 10.0
    cerebro_scale_mode: str = "both"
    cerebro_current_metrics: tuple[str, ...] = ("total_abs",)
    visualization_relative_floor: float = 0.0
    make_gif: bool = False
    skip_cerebro: bool = False
    run_all_stimulations: bool = False
    run_ftract: bool = False
    peak_file: Path | None = None
    peak_variable: str = "I_node_max"


@dataclass
class RCData:
    SC: np.ndarray
    D: np.ndarray
    F: np.ndarray
    R: np.ndarray
    G: np.ndarray
    laplacian: np.ndarray
    system_matrix: np.ndarray
    node_names: list[str]
    raw_order: np.ndarray
    n_cortex: int


@dataclass
class DynamicResult:
    time_ms: np.ndarray
    voltage: np.ndarray
    current_abs_total: np.ndarray
    current_net: np.ndarray
    current_positive_inflow: np.ndarray
    stim_indices: np.ndarray
    stim_labels: list[str]


@dataclass
class FtractData:
    amplitude: np.ndarray
    probability: np.ndarray
    n_implantations: np.ndarray
    n_stimulations: np.ndarray
    n_imp: np.ndarray


def load_aligned_ystruct_v5(config: Config, data: RCData) -> FtractData:
    """Load the no-reorder empirical matrices paired with X_struct.v8."""
    path = require_file(rcc.ystruct_path(config.data_root))
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "Y_struct" not in mat or not hasattr(mat["Y_struct"], "v5"):
        raise KeyError(f"Y_struct.v5 not found in {path}")
    values = mat["Y_struct"].v5
    required = ("A", "P", "Ni", "Ns", "Nr", "parcel_list")
    missing = [name for name in required if not hasattr(values, name)]
    if missing:
        raise KeyError(f"Missing Y_struct.v5 fields in {path}: {missing}")

    matrices = {
        name: np.asarray(getattr(values, name), dtype=float)
        for name in ("A", "P", "Ni", "Ns", "Nr")
    }
    expected = (data.n_cortex, data.n_cortex)
    for name, matrix in matrices.items():
        if matrix.shape != expected:
            raise ValueError(f"Y_struct.v5.{name} is {matrix.shape}; expected {expected}.")

    empirical_names = matlab_cellstr_to_list(values.parcel_list)
    if empirical_names != data.node_names[: data.n_cortex]:
        raise ValueError(
            "Y_struct.v5 parcel_list does not match X_struct.v8 cortical order; "
            "a verified permutation is required."
        )

    estimated_responses = matrices["Ns"] * matrices["P"]
    valid = np.isfinite(matrices["Nr"]) & np.isfinite(estimated_responses)
    if not np.allclose(matrices["Nr"][valid], estimated_responses[valid]):
        raise ValueError("Y_struct.v5.Nr differs from Ns*P at finite entries.")

    # Preserve the old variable convention: N_imp means implantation count.
    return FtractData(
        amplitude=matrices["A"],
        probability=matrices["P"],
        n_implantations=matrices["Ni"],
        n_stimulations=matrices["Ns"],
        n_imp=matrices["Ni"],
    )


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return path


def matlab_cellstr_to_list(values: np.ndarray) -> list[str]:
    """Convert scipy's representation of a MATLAB cellstr to Python strings."""
    return [str(value).strip() for value in np.asarray(values).reshape(-1).tolist()]


def validate_permutation(values: np.ndarray, length: int, name: str) -> np.ndarray:
    values = np.asarray(values).reshape(-1)
    if values.size != length:
        raise ValueError(f"{name} has {values.size} entries; expected {length}.")
    if not np.all(np.isfinite(values)) or not np.all(values == np.round(values)):
        raise ValueError(f"{name} must contain finite integer indices.")
    zero_based = values.astype(int) - 1
    if not np.array_equal(np.sort(zero_based), np.arange(length)):
        raise ValueError(f"{name} is not a permutation of MATLAB indices 1:{length}.")
    return zero_based


def load_rc_data(config: Config) -> RCData:
    sc_path = require_file(config.data_root / "SC_data_neuron_2023_yihang.mat")

    mat = loadmat(sc_path, squeeze_me=True, struct_as_record=False)
    missing = {"SC", "D", "F"}.difference(mat)
    if missing:
        raise KeyError(f"Missing variables in {sc_path}: {sorted(missing)}")

    SC_raw = np.asarray(mat["SC"], dtype=float)
    D_raw = np.asarray(mat["D"], dtype=float)
    F_raw = np.asarray(mat["F"], dtype=float)
    if SC_raw.ndim != 2 or SC_raw.shape[0] != SC_raw.shape[1]:
        raise ValueError("SC must be square.")
    if D_raw.shape != SC_raw.shape or F_raw.shape != SC_raw.shape:
        raise ValueError("SC, D, and F must have identical dimensions.")

    if "parcel_list" not in mat:
        raise KeyError(f"parcel_list not found in {sc_path}")
    raw_names = np.asarray(matlab_cellstr_to_list(mat["parcel_list"]), dtype=object)
    if raw_names.size != SC_raw.shape[0]:
        raise ValueError("parcel_list length does not match SC dimensions.")

    if config.model_variant == "main_caio":
        order_path = require_file(config.data_root / "ReOrder_Caio_SC.mat")
        order_mat = loadmat(order_path, squeeze_me=True)
        if "ReOrder_Caio_SC" not in order_mat:
            raise KeyError(f"ReOrder_Caio_SC not found in {order_path}")
        order = validate_permutation(
            order_mat["ReOrder_Caio_SC"], SC_raw.shape[0], "ReOrder_Caio_SC"
        )
        SC = SC_raw[np.ix_(order, order)].copy()
        D = D_raw[np.ix_(order, order)].copy()
        F = F_raw[np.ix_(order, order)].copy()
        node_names = raw_names[order].astype(str).tolist()
        capacitance_value = 1.0
        apply_fco_threshold = True
    elif config.model_variant == "legacy_xstruct_v8":
        xstruct_v8 = rcc.load_xstruct_v8(require_file(rcc.xstruct_path(config.data_root)))
        SC = np.asarray(xstruct_v8.W, dtype=float).copy()
        D = np.asarray(xstruct_v8.ED, dtype=float).copy()
        if SC.shape != SC_raw.shape or D.shape != SC_raw.shape:
            raise ValueError(
                "X_struct.v8.W and ED must match the 374-node SC_data dimensions."
            )
        # This reproduces RC_ftract_Caio_data.m: raw 1:374 order, C=0.1,
        # and no FCO edge threshold. Keep F only as aligned metadata.
        F = F_raw.copy()
        order = np.arange(SC.shape[0], dtype=int)
        node_names = raw_names.astype(str).tolist()
        capacitance_value = 0.1
        apply_fco_threshold = False
    else:
        raise ValueError(
            "model_variant must be 'main_caio' or 'legacy_xstruct_v8'."
        )

    SC[~np.isfinite(SC)] = 0.0
    D[~np.isfinite(D)] = 0.0
    if apply_fco_threshold:
        SC[F < config.fco_threshold] = 0.0

    R = np.zeros_like(SC)
    valid_sc = SC != 0
    R[valid_sc] = D[valid_sc] / SC[valid_sc]
    R[~np.isfinite(R)] = 0.0

    G = np.zeros_like(R)
    valid_r = R > 0
    G[valid_r] = 1.0 / R[valid_r]

    laplacian = np.diag(G.sum(axis=1)) - G
    capacitance = np.full(SC.shape[0], capacitance_value, dtype=float)
    system_matrix = -(laplacian / capacitance[:, None])

    if len(node_names) != 374 or config.n_cortex != 360:
        raise ValueError(
            "This Tian S1 implementation expects 360 cortical + 14 subcortical nodes."
        )
    expected_subcortex = list(S1_NATIVE_LABEL_TRANSLATION)
    if node_names[config.n_cortex :] != expected_subcortex:
        raise ValueError(
            "Subcortical order does not match the expected Tian S1 14-node order.\n"
            f"Found: {node_names[config.n_cortex :]}"
        )

    return RCData(
        SC=SC,
        D=D,
        F=F,
        R=R,
        G=G,
        laplacian=laplacian,
        system_matrix=system_matrix,
        node_names=node_names,
        raw_order=order,
        n_cortex=config.n_cortex,
    )


def resolve_stim_indices(node_names: Sequence[str], labels: Sequence[str]) -> np.ndarray:
    indices: list[int] = []
    for label in labels:
        matches = [i for i, name in enumerate(node_names) if name == label]
        if not matches:
            suggestions = [name for name in node_names if label.lower() in name.lower()]
            raise ValueError(
                f"Stimulus label not found: {label}. Similar labels: {suggestions[:12]}"
            )
        indices.append(matches[0])
    return np.asarray(indices, dtype=int)


def current_timecourses(
    G: np.ndarray,
    laplacian: np.ndarray,
    voltage: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return total-absolute, signed-net, and gross-positive-inflow currents."""
    n, n_time = voltage.shape
    row, col = np.nonzero(G > 0)
    conductance = G[row, col]
    absolute_total = np.zeros((n, n_time), dtype=float)
    positive_inflow = np.zeros((n, n_time), dtype=float)

    for time_index in range(n_time):
        outward_edge_current = conductance * (
            voltage[row, time_index] - voltage[col, time_index]
        )
        absolute_total[:, time_index] = np.bincount(
            row, weights=np.abs(outward_edge_current), minlength=n
        )
        # Positive current into node i from neighbour j is G_ij * (V_j - V_i).
        # Sum positive edge inflows directly; this is not clip(net_current, 0).
        positive_inflow[:, time_index] = np.bincount(
            row, weights=np.maximum(-outward_edge_current, 0.0), minlength=n
        )

    signed_net = -laplacian @ voltage
    return absolute_total, signed_net, positive_inflow


def simulate_dynamic(data: RCData, config: Config) -> DynamicResult:
    stim_indices = resolve_stim_indices(data.node_names, config.stim_labels)
    initial_voltage = np.zeros(len(data.node_names), dtype=float)
    voltage_per_site = config.initial_voltage
    if config.split_multisite_voltage and stim_indices.size > 1:
        voltage_per_site /= stim_indices.size
    initial_voltage[stim_indices] = voltage_per_site

    time_ms = np.linspace(0.0, config.time_scale_ms, config.num_points)
    if config.simulation_solver == "exact_eigh":
        if not np.allclose(data.system_matrix, data.system_matrix.T, rtol=1e-10, atol=1e-12):
            raise ValueError("exact_eigh requires a symmetric RC system matrix.")
        eigenvalues, eigenvectors = np.linalg.eigh(data.system_matrix)
        modal_initial = eigenvectors.T @ initial_voltage
        voltage = eigenvectors @ (
            modal_initial[:, None] * np.exp(eigenvalues[:, None] * time_ms[None, :])
        )
    elif config.simulation_solver == "rk45":
        solution = solve_ivp(
            lambda _time, voltage: data.system_matrix @ voltage,
            (time_ms[0], time_ms[-1]),
            initial_voltage,
            method="RK45",
            t_eval=time_ms,
            rtol=config.rtol,
            atol=config.atol,
            vectorized=False,
        )
        if not solution.success:
            raise RuntimeError(f"RC integration failed: {solution.message}")
        voltage = solution.y
    else:
        raise ValueError("simulation_solver must be 'rk45' or 'exact_eigh'.")

    current_abs_total, current_net, current_positive_inflow = current_timecourses(
        data.G, data.laplacian, voltage
    )
    return DynamicResult(
        time_ms=time_ms,
        voltage=voltage,
        current_abs_total=current_abs_total,
        current_net=current_net,
        current_positive_inflow=current_positive_inflow,
        stim_indices=stim_indices,
        stim_labels=[data.node_names[index] for index in stim_indices],
    )


def simulate_all_peak_currents(data: RCData, config: Config) -> np.ndarray:
    """Compute the MATLAB-compatible I_node_max matrix for all stimulations."""
    n = len(data.node_names)
    time_ms = np.linspace(0.0, config.time_scale_ms, config.num_points)
    peak = np.zeros((n, n), dtype=float)

    for stim in range(n):
        initial_voltage = np.zeros(n, dtype=float)
        initial_voltage[stim] = config.initial_voltage
        solution = solve_ivp(
            lambda _time, voltage: data.system_matrix @ voltage,
            (time_ms[0], time_ms[-1]),
            initial_voltage,
            method="RK45",
            t_eval=time_ms,
            rtol=config.rtol,
            atol=config.atol,
        )
        if not solution.success:
            raise RuntimeError(f"RC integration failed at stim {stim + 1}: {solution.message}")
        current_abs_total, _, _ = current_timecourses(
            data.G, data.laplacian, solution.y
        )
        peak[:, stim] = np.max(current_abs_total, axis=1)
        print(f"Stim node {stim + 1}/{n} done.", flush=True)
    return peak


def read_ftract_matrix(directory: Path, base_name: str, expected_size: int = 360) -> np.ndarray:
    txt = directory / f"{base_name}.txt"
    gz = directory / f"{base_name}.txt.gz"
    if txt.is_file():
        opener = lambda: txt.open("rt", encoding="utf-8")
    elif gz.is_file():
        opener = lambda: gzip.open(gz, "rt", encoding="utf-8")
    else:
        raise FileNotFoundError(f"Neither {txt} nor {gz} exists.")

    rows: list[np.ndarray] = []
    with opener() as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(np.fromstring(line, sep=" ", dtype=float))

    if len(rows) != expected_size or any(row.size != expected_size for row in rows):
        widths = sorted({row.size for row in rows})
        raise ValueError(
            f"{base_name} must be {expected_size}x{expected_size}; "
            f"found {len(rows)} rows with widths {widths}."
        )
    return np.vstack(rows)


def load_ftract_data(config: Config) -> FtractData:
    order_path = require_file(config.data_root / "ReOreder_F_tract.mat")
    order_mat = loadmat(order_path, squeeze_me=True)
    if "ReOrder_F_tract" not in order_mat:
        raise KeyError(f"ReOrder_F_tract not found in {order_path}")
    order_table = np.asarray(order_mat["ReOrder_F_tract"])
    if order_table.ndim != 2 or order_table.shape[1] < 2:
        raise ValueError("ReOrder_F_tract must contain at least two columns.")
    order = validate_permutation(
        order_table[:, 1], config.n_cortex, "ReOrder_F_tract(:,2)"
    )

    directory = config.ftract_root / config.parcellation
    amplitude = read_ftract_matrix(directory, "amplitude__mad")
    probability = read_ftract_matrix(directory, "probability")
    n_implantations = read_ftract_matrix(directory, "N_implantations")
    n_stimulations = read_ftract_matrix(directory, "N_stimulations")

    ix = np.ix_(order, order)
    amplitude = amplitude[ix]
    probability = probability[ix]
    n_implantations = n_implantations[ix]
    n_stimulations = n_stimulations[ix]
    n_imp = n_stimulations * probability
    return FtractData(
        amplitude=amplitude,
        probability=probability,
        n_implantations=n_implantations,
        n_stimulations=n_stimulations,
        n_imp=n_imp,
    )


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 2:
        return math.nan, math.nan
    result = spearmanr(x[valid], y[valid])
    return float(result.statistic), float(result.pvalue)


def partial_spearman(x: np.ndarray, y: np.ndarray, covariate: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(covariate)
    if valid.sum() < 4:
        return math.nan, math.nan

    x_rank = rankdata(x[valid])
    y_rank = rankdata(y[valid])
    z_rank = rankdata(covariate[valid])
    design = np.column_stack((np.ones(valid.sum()), z_rank))
    x_residual = x_rank - design @ np.linalg.lstsq(design, x_rank, rcond=None)[0]
    y_residual = y_rank - design @ np.linalg.lstsq(design, y_rank, rcond=None)[0]
    if np.allclose(x_residual, 0) or np.allclose(y_residual, 0):
        return math.nan, math.nan
    result = pearsonr(x_residual, y_residual)
    return float(result.statistic), float(result.pvalue)


def threshold_sweep(
    empirical: np.ndarray,
    peak_current: np.ndarray,
    data: RCData,
    ftract: FtractData,
    thresholds: Iterable[int] = range(101),
    low_fco_only: bool = False,
) -> list[dict[str, float]]:
    n = data.n_cortex
    D = data.D[:n, :n]
    SC = data.SC[:n, :n]
    R = data.R[:n, :n]
    F = data.F[:n, :n]
    current = peak_current[:n, :n]
    hemisphere = np.block(
        [
            [np.ones((n // 2, n // 2), dtype=bool), np.zeros((n // 2, n // 2), dtype=bool)],
            [np.zeros((n // 2, n // 2), dtype=bool), np.ones((n // 2, n // 2), dtype=bool)],
        ]
    )

    rows: list[dict[str, float]] = []
    for threshold in thresholds:
        mask = (ftract.n_imp > threshold) & (D >= 20) & hemisphere
        if low_fco_only:
            mask &= F < 0.25
        observed = empirical[mask]
        distance = D[mask]
        vectors = {
            "distance": -distance,
            "rc": current[mask],
            "sc": SC[mask],
            "resistance": -R[mask],
        }
        row: dict[str, float] = {"threshold": float(threshold), "count": float(mask.sum())}
        for name, values in vectors.items():
            row[f"r_{name}"], row[f"p_{name}"] = safe_spearman(values, observed)
            if name != "distance":
                row[f"partial_r_{name}"], row[f"partial_p_{name}"] = partial_spearman(
                    observed, values, distance
                )
        rows.append(row)
    return rows


def load_peak_current(path: Path, variable: str, n_nodes: int) -> np.ndarray:
    require_file(path)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        values = np.load(path)
    elif suffix == ".npz":
        archive = np.load(path, allow_pickle=False)
        if variable not in archive:
            raise KeyError(f"{variable} not found in {path}; available: {archive.files}")
        values = archive[variable]
    elif suffix == ".mat":
        mat = loadmat(path, squeeze_me=True)
        if variable not in mat:
            raise KeyError(f"{variable} not found in {path}")
        values = mat[variable]
    else:
        raise ValueError("Peak-current file must be .npy, .npz, or non-v7.3 .mat.")
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[0] < n_nodes or values.shape[1] < n_nodes:
        raise ValueError(f"Peak-current matrix must be at least {n_nodes}x{n_nodes}.")
    return values


def save_sweep_csv(rows: Sequence[dict[str, float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_threshold_sweep(rows: Sequence[dict[str, float]], path: Path, title: str) -> None:
    threshold = np.asarray([row["threshold"] for row in rows])
    count = np.asarray([row["count"] for row in rows])
    methods = ("distance", "rc", "sc", "resistance")
    labels = {"distance": "Distance", "rc": "RC current", "sc": "SC", "resistance": "Resistance"}

    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=True)
    for method in methods:
        axes[0].plot(threshold, [row[f"r_{method}"] for row in rows], label=labels[method])
        axes[1].plot(threshold, [row[f"p_{method}"] for row in rows], label=labels[method])
    axes[0].set_ylabel("Spearman r")
    axes[0].set_title(title)
    axes[0].legend(frameon=False, ncol=2)
    axes[1].set_ylabel("p-value")
    axes[1].set_yscale("log")
    axes[2].plot(threshold, count, color="black")
    axes[2].set_ylabel("Number of pairs")
    axes[2].set_xlabel("Estimated N_resp threshold")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def positive_log_transform(current: np.ndarray) -> tuple[np.ndarray, float]:
    positive = current[np.isfinite(current) & (current > 0)]
    scale = float(np.median(positive)) if positive.size else 1.0
    if scale <= 0 or not np.isfinite(scale):
        scale = 1.0
    return np.log1p(np.maximum(current, 0) / scale), scale


def signed_log_transform(current: np.ndarray) -> tuple[np.ndarray, float]:
    magnitude = np.abs(current[np.isfinite(current) & (current != 0)])
    scale = float(np.median(magnitude)) if magnitude.size else 1.0
    if scale <= 0 or not np.isfinite(scale):
        scale = 1.0
    transformed = np.sign(current) * np.log1p(np.abs(current) / scale)
    return transformed, scale


def plot_dynamic_summary(result: DynamicResult, data: RCData, path: Path) -> None:
    transformed, scale = positive_log_transform(result.current_abs_total)
    time = result.time_ms
    subcortex = result.current_abs_total[data.n_cortex :, :]

    fig, axes = plt.subplots(3, 1, figsize=(13, 13), gridspec_kw={"height_ratios": [2.2, 1, 1.8]})
    image = axes[0].imshow(
        transformed,
        aspect="auto",
        origin="lower",
        extent=(time[0], time[-1], 1, transformed.shape[0]),
        cmap="magma",
        interpolation="nearest",
    )
    axes[0].axhline(data.n_cortex + 0.5, color="cyan", linewidth=1.2)
    axes[0].set_ylabel("Node (1–360 cortex; 361–374 subcortex)")
    axes[0].set_title(
        "DLPFC RC propagation: log1p(total absolute node current / median positive current)"
    )
    fig.colorbar(image, ax=axes[0], label=f"log1p(I / {scale:.3g})")

    cortical_total = result.current_abs_total[: data.n_cortex].sum(axis=0)
    subcortical_total = subcortex.sum(axis=0)
    axes[1].plot(time, cortical_total, label="Cortex (360 nodes)", linewidth=2)
    axes[1].plot(time, subcortical_total, label="Subcortex (14 nodes)", linewidth=2)
    axes[1].set_ylabel("Total absolute current")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.25)

    colors = plt.cm.tab20(np.linspace(0, 1, subcortex.shape[0]))
    for index, (name, color) in enumerate(
        zip(data.node_names[data.n_cortex :], colors, strict=True)
    ):
        axes[2].plot(time, subcortex[index], label=name, color=color, linewidth=1.2)
    axes[2].set_xlabel("Time (ms)")
    axes[2].set_ylabel("Absolute node current")
    axes[2].set_title("Tian S1 subcortical current dynamics")
    axes[2].grid(alpha=0.25)
    axes[2].legend(frameon=False, fontsize=8, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.suptitle("Stimulation: " + ", ".join(result.stim_labels), fontsize=15)
    fig.tight_layout(rect=(0, 0, 0.83, 0.97))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


class CerebroRenderer:
    """Restore 374-node results, expand them to 91k, and render."""

    def __init__(
        self,
        atlas_file: Path,
        cortical_atlas_file: Path,
        node_names: Sequence[str],
        model_to_raw_order: np.ndarray,
        mapping_mode: str = "matlab_sc_order",
    ):
        require_file(atlas_file)
        require_file(cortical_atlas_file)
        self.atlas_file = atlas_file
        self.cortical_atlas_file = cortical_atlas_file
        self.mapping_mode = mapping_mode
        if mapping_mode not in {"matlab_sc_order", "native_tian_s1"}:
            raise ValueError(
                "mapping_mode must be 'matlab_sc_order' or 'native_tian_s1'."
            )
        self.node_names = list(node_names)
        self.model_to_raw_order = np.asarray(model_to_raw_order, dtype=int).reshape(-1)
        if self.model_to_raw_order.size != len(self.node_names):
            raise ValueError("Model-to-raw permutation length does not match node names.")
        if not np.array_equal(
            np.sort(self.model_to_raw_order), np.arange(len(self.node_names))
        ):
            raise ValueError("Model-to-raw order is not a zero-based permutation.")

        raw_names = np.empty(len(self.node_names), dtype=object)
        raw_names[self.model_to_raw_order] = self.node_names
        self.raw_node_names = raw_names.astype(str).tolist()

        from cerebro import cerebro_brain_utils as cbu

        template = nib.load(cbu.cifti_template_file)
        self.dense_size = int(np.prod(template.shape))

        brain_models = list(template.header.get_index_map(1).brain_models)
        structure_names: list[str] = []
        for model in brain_models:
            structure_names.extend([model.brain_structure] * model.index_count)
        self.structure_names = np.asarray(structure_names)
        if self.structure_names.size != self.dense_size:
            raise ValueError(
                "Cerebro brain-model labels do not match the dense CIFTI length: "
                f"{self.structure_names.size} vs {self.dense_size}."
            )

        atlas = nib.load(atlas_file)
        self.atlas_dense_labels = np.asarray(atlas.get_fdata()).squeeze().astype(int)
        if self.atlas_dense_labels.size != self.dense_size:
            raise ValueError(
                "Atlas and Cerebro template have different dense CIFTI lengths: "
                f"{self.atlas_dense_labels.size} vs {self.dense_size}."
            )

        label_axis = atlas.header.get_axis(0)
        label_table = label_axis.label[0]
        self.atlas_name_to_key = {
            str(label_entry[0]): int(key) for key, label_entry in label_table.items()
        }

        self.cortical_dense_labels: np.ndarray | None = None
        if self.mapping_mode == "matlab_sc_order":
            from Connectome_Spatial_Smoothing import CSS as css

            self.cortical_dense_labels = np.asarray(
                css._label_surface_atlas_to_data(
                    cifti=template,
                    atlas_file=str(cortical_atlas_file),
                    subcortex=False,
                ),
                dtype=object,
            )
            if self.cortical_dense_labels.size != self.dense_size:
                raise ValueError(
                    "Cortical atlas projection does not match the 91k template."
                )

        self.node_to_dense = self._build_node_mapping()
        self.left_mask = np.asarray(["LEFT" in name for name in self.structure_names])
        self.right_mask = np.asarray(["RIGHT" in name for name in self.structure_names])

    def _build_node_mapping(self) -> list[np.ndarray]:
        mapping: list[np.ndarray] = []
        missing: list[str] = []
        occupied = np.zeros(self.dense_size, dtype=bool)
        for raw_index, node_name in enumerate(self.raw_node_names):
            if self.mapping_mode == "matlab_sc_order":
                if raw_index < 360:
                    dense_indices = np.flatnonzero(
                        self.cortical_dense_labels == node_name
                    )
                else:
                    structure_name = S1_CIFTI_STRUCTURE_TRANSLATION.get(node_name)
                    dense_indices = (
                        np.flatnonzero(self.structure_names == structure_name)
                        if structure_name is not None
                        else np.empty(0, dtype=int)
                    )
            else:
                atlas_names = S1_NATIVE_LABEL_TRANSLATION.get(node_name, (node_name,))
                missing_labels = [
                    name for name in atlas_names if name not in self.atlas_name_to_key
                ]
                if missing_labels:
                    mapping.append(np.empty(0, dtype=int))
                    missing.append(f"{node_name} -> {missing_labels}")
                    continue
                keys = np.asarray(
                    [self.atlas_name_to_key[name] for name in atlas_names], dtype=int
                )
                dense_indices = np.flatnonzero(
                    np.isin(self.atlas_dense_labels, keys)
                )
            if dense_indices.size:
                if np.any(occupied[dense_indices]):
                    raise ValueError(
                        f"Dense atlas mapping overlaps at node {node_name}."
                    )
                occupied[dense_indices] = True
                mapping.append(dense_indices)
            else:
                mapping.append(np.empty(0, dtype=int))
                missing.append(node_name)
        if missing:
            raise ValueError(
                f"Nodes missing from {self.mapping_mode} mapping: {missing}"
            )
        return mapping

    def to_dense(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=float).reshape(-1)
        if values.size != len(self.node_names):
            raise ValueError(
                f"Parcel vector has {values.size} values; expected {len(self.node_names)}."
            )
        raw_values = np.empty_like(values)
        raw_values[self.model_to_raw_order] = values
        dense = np.full(self.dense_size, np.nan, dtype=float)
        for value, indices in zip(raw_values, self.node_to_dense, strict=True):
            dense[indices] = value
        return dense

    @staticmethod
    def _render_view(
        dense: np.ndarray,
        anatomy_mask: np.ndarray,
        axis: plt.Axes,
        view: object,
        clims: tuple[float, float],
        cmap: object,
        separation: float,
        zoom: float = 0.76,
        volumetric_structures: str | None = "subcortex",
    ) -> None:
        from cerebro import cerebro_brain_viewer as cbv

        viewer = None
        try:
            viewer = cbv.Cerebro_brain_viewer(
                offscreen=True,
                background_color=(1.0, 1.0, 1.0, 0.0),
                null_color=(1.0, 1.0, 1.0, 0.0),
                no_color=(0.7, 0.7, 0.7, 1.0),
            )
            viewer.load_template_GIFTI_cortical_surface_models("inflated")
            viewer.visualize_cifti_space(
                volumetric_structures=volumetric_structures,
                cifti_expansion_scale=20,
                cifti_left_right_seperation=separation,
                volumetric_structure_offset=(0, 0, -80),
                volume_rendering="surface",
                smoothing_filter="laplacian",
                smoothing=40,
            )
            anatomy = np.full(dense.shape, np.nan, dtype=float)
            anatomy[anatomy_mask] = 0.5
            viewer.add_cifti_dscalar_layer(
                dscalar_data=anatomy,
                colormap=ListedColormap(["#B8B8B8", "#B8B8B8"]),
                clims=(0.0, 1.0),
                exclusion_color=(1.0, 1.0, 1.0, 0.0),
                opacity=1.0,
            )
            viewer.add_cifti_dscalar_layer(
                dscalar_data=dense,
                colormap=cmap,
                clims=clims,
                exclusion_color=(1.0, 1.0, 1.0, 0.0),
                opacity=0.95,
            )
            axis.axis("off")
            camera = viewer._view_to_camera_config(view)
            camera = viewer._zoom_camera_to_content(camera)
            camera["camera_pos"] = tuple(value / zoom for value in camera["camera_pos"])
            viewer.viewer.change_view(**camera)
            viewer.offscreen_draw_to_matplotlib_axes(axis)
        finally:
            if viewer is not None:
                viewer.viewer.window.destroy()

    def render_four_views(
        self,
        values: np.ndarray,
        path: Path,
        title: str,
        clims: tuple[float, float],
        cmap_name: str = "magma",
        colorbar_label: str = "Model current (a.u.)",
        large_colorbar: bool = False,
        font_scale: float = 1.0,
        font_weight: str = "normal",
    ) -> None:
        if font_scale <= 0:
            raise ValueError("font_scale must be positive.")
        dense = self.to_dense(values)
        left = dense.copy()
        left[~self.left_mask] = np.nan
        right = dense.copy()
        right[~self.right_mask] = np.nan
        cmap = plt.get_cmap(cmap_name)

        # Leave generous vertical space around every Cerebro view.  A camera
        # zoom below 1 moves the camera away from the content and prevents the
        # superior cortex from being clipped by the per-panel viewport.
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        self._render_view(
            left,
            self.left_mask,
            axes[0, 0],
            ((-420, 0, 0), None, None, None),
            clims,
            cmap,
            separation=100,
            zoom=0.78,
        )
        self._render_view(
            left,
            self.left_mask,
            axes[1, 0],
            ((420, 0, 0), None, None, None),
            clims,
            cmap,
            separation=-100,
        )
        self._render_view(
            right,
            self.right_mask,
            axes[0, 1],
            ((420, 0, 0), None, None, None),
            clims,
            cmap,
            separation=100,
            zoom=0.78,
        )
        self._render_view(
            right,
            self.right_mask,
            axes[1, 1],
            ((-420, 0, 0), None, None, None),
            clims,
            cmap,
            separation=-100,
        )

        scalar = plt.cm.ScalarMappable(norm=Normalize(*clims), cmap=cmap)
        scalar.set_array([])
        colorbar = fig.colorbar(
            scalar,
            ax=axes.ravel().tolist(),
            orientation="horizontal",
            fraction=0.055 if large_colorbar else 0.035,
            pad=0.018 if large_colorbar else 0.02,
        )
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))
        colorbar.formatter = formatter
        colorbar.update_ticks()
        if large_colorbar:
            colorbar.ax.tick_params(
                labelsize=14 * font_scale,
                width=1.5,
                length=6,
            )
            for tick_label in colorbar.ax.get_xticklabels():
                tick_label.set_fontweight(font_weight)
            offset_text = colorbar.ax.xaxis.get_offset_text()
            offset_text.set_fontsize(13 * font_scale)
            offset_text.set_fontweight(font_weight)
            colorbar.set_label(
                colorbar_label,
                fontsize=15 * font_scale,
                fontweight=font_weight,
                labelpad=3,
            )
            colorbar.outline.set_linewidth(1.5)
        else:
            colorbar.set_label(
                colorbar_label,
                fontsize=plt.rcParams["axes.labelsize"] * font_scale
                if isinstance(plt.rcParams["axes.labelsize"], (int, float))
                else None,
                fontweight=font_weight,
            )
            if font_scale != 1.0 or font_weight != "normal":
                colorbar.ax.tick_params(
                    labelsize=plt.rcParams["xtick.labelsize"] * font_scale
                    if isinstance(plt.rcParams["xtick.labelsize"], (int, float))
                    else None
                )
                for tick_label in colorbar.ax.get_xticklabels():
                    tick_label.set_fontweight(font_weight)
        fig.suptitle(
            title,
            fontsize=14 * font_scale,
            fontweight=font_weight,
            y=0.965,
        )
        fig.subplots_adjust(
            left=0.04,
            right=0.96,
            bottom=0.105,
            top=0.90,
            wspace=-0.12,
            hspace=-0.12,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, transparent=False)
        plt.close(fig)


def nearest_time_indices(time_ms: np.ndarray, requested: Sequence[float]) -> np.ndarray:
    indices = [int(np.argmin(np.abs(time_ms - value))) for value in requested]
    return np.asarray(list(dict.fromkeys(indices)), dtype=int)


def save_gif(frame_paths: Sequence[Path], output: Path, duration_ms: int = 160) -> None:
    from PIL import Image

    source_images: list[Image.Image] = []
    for frame_path in frame_paths:
        with Image.open(frame_path) as frame:
            source_images.append(frame.convert("RGB").copy())
    if not source_images:
        raise ValueError("Cannot save a GIF without frames.")

    # bbox_inches="tight" can make frame widths differ with title length.
    # Centre-pad every image to one canvas so the animation does not jitter.
    width = max(image.width for image in source_images)
    height = max(image.height for image in source_images)
    images: list[Image.Image] = []
    for source in source_images:
        canvas = Image.new("RGB", (width, height), "white")
        canvas.paste(source, ((width - source.width) // 2, (height - source.height) // 2))
        images.append(canvas)
    images[0].save(
        output,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        optimize=True,
    )


def render_cerebro_dynamics(
    result: DynamicResult,
    data: RCData,
    config: Config,
) -> list[Path]:
    renderer = CerebroRenderer(
        config.atlas_file,
        config.cortical_atlas_file,
        data.node_names,
        data.raw_order,
        config.dense_mapping_mode,
    )
    snapshot_indices = nearest_time_indices(result.time_ms, config.snapshot_times_ms)
    if config.make_gif:
        requested_animation_times = np.arange(
            0.0,
            config.time_scale_ms + config.animation_step_ms * 0.5,
            config.animation_step_ms,
        )
        animation_indices = nearest_time_indices(
            result.time_ms, requested_animation_times
        )
        render_indices = np.unique(
            np.concatenate((snapshot_indices, animation_indices))
        )
    else:
        animation_indices = np.empty(0, dtype=int)
        render_indices = snapshot_indices
    animation_index_set = set(animation_indices.tolist())

    output_paths: list[Path] = []
    render_fixed = config.cerebro_scale_mode in {"fixed", "both"}
    render_framewise = config.cerebro_scale_mode in {"framewise", "both"}

    metric_specs = {
        "total_abs": (
            result.current_abs_total,
            "Total absolute model current",
            "total_abs_current",
            False,
            "magma",
        ),
        "net": (
            result.current_net,
            "Signed net model current (positive=inflow)",
            "net_current",
            True,
            "coolwarm",
        ),
        "net_positive": (
            np.where(result.current_net > 0.0, result.current_net, np.nan),
            "Positive-only net model current",
            "net_positive_current",
            False,
            "magma",
        ),
        "positive_inflow": (
            result.current_positive_inflow,
            "Positive inflow model current",
            "positive_inflow_current",
            False,
            "magma",
        ),
    }

    for metric_name in config.cerebro_current_metrics:
        raw_values, display_name, slug, signed, cmap_name = metric_specs[metric_name]
        # Brain maps use raw model-current values so the colorbar can report
        # the actual a.u. magnitude and its scientific-notation exponent.
        display_values = raw_values
        finite = display_values[np.isfinite(display_values)]
        fixed_upper = float(np.max(np.abs(finite))) if finite.size else 1.0
        fixed_upper = max(fixed_upper, np.finfo(float).eps)
        display_floor = fixed_upper * config.visualization_relative_floor
        fixed_clims = (-fixed_upper, fixed_upper) if signed else (0.0, fixed_upper)

        fixed_paths: list[Path] = []
        framewise_paths: list[Path] = []
        for index in render_indices:
            time_ms = result.time_ms[index]
            filename = f"current_t{time_ms:07.2f}ms.png"

            if render_fixed:
                fixed_output = (
                    config.output_dir
                    / f"cerebro_{slug}_fixed_scale_frames"
                    / filename
                )
                renderer.render_four_views(
                    display_values[:, index],
                    fixed_output,
                    title=(
                        f"DLPFC {display_name.lower()} at {time_ms:.2f} ms "
                        f"(a.u.; fixed global colorbar)"
                    ),
                    clims=fixed_clims,
                    cmap_name=cmap_name,
                )
                if index in animation_index_set:
                    fixed_paths.append(fixed_output)
                output_paths.append(fixed_output)
                print(f"Rendered {fixed_output}", flush=True)

            if render_framewise:
                raw_frame = raw_values[:, index]
                raw_finite = raw_frame[np.isfinite(raw_frame)]
                raw_max = (
                    float(np.max(np.abs(raw_finite)))
                    if signed and raw_finite.size
                    else float(np.max(raw_finite))
                    if raw_finite.size
                    else 0.0
                )
                below_display_floor = (
                    display_floor > 0.0 and raw_max <= display_floor
                )
                frame_values = (
                    np.zeros_like(raw_frame) if below_display_floor else raw_frame
                )
                frame_finite = frame_values[np.isfinite(frame_values)]
                frame_upper = (
                    float(np.max(np.abs(frame_finite))) if frame_finite.size else 1.0
                )
                frame_upper = max(
                    frame_upper, display_floor, np.finfo(float).eps
                )
                frame_clims = (
                    (-frame_upper, frame_upper) if signed else (0.0, frame_upper)
                )
                framewise_output = (
                    config.output_dir
                    / f"cerebro_{slug}_framewise_frames"
                    / filename
                )
                renderer.render_four_views(
                    frame_values,
                    framewise_output,
                    title=(
                        f"DLPFC {display_name.lower()} at {time_ms:.2f} ms "
                        f"(a.u.; below display floor {display_floor:.3g})"
                        if below_display_floor
                        else
                        f"DLPFC {display_name.lower()} at {time_ms:.2f} ms "
                        f"(a.u.; framewise |max|={raw_max:.3g})"
                        if signed
                        else f"DLPFC {display_name.lower()} at {time_ms:.2f} ms "
                        f"(a.u.; framewise max={raw_max:.3g})"
                    ),
                    clims=frame_clims,
                    cmap_name=cmap_name,
                )
                if index in animation_index_set:
                    framewise_paths.append(framewise_output)
                output_paths.append(framewise_output)
                print(f"Rendered {framewise_output}", flush=True)

        if config.make_gif and render_fixed:
            fixed_gif = config.output_dir / f"DLPFC_{slug}_fixed_scale.gif"
            save_gif(fixed_paths, fixed_gif)
            output_paths.append(fixed_gif)
            print(f"Saved {fixed_gif}", flush=True)
        if config.make_gif and render_framewise:
            framewise_gif = config.output_dir / f"DLPFC_{slug}_framewise.gif"
            save_gif(framewise_paths, framewise_gif)
            output_paths.append(framewise_gif)
            print(f"Saved {framewise_gif}", flush=True)
    return output_paths


def save_dynamic_result(result: DynamicResult, data: RCData, config: Config) -> Path:
    output = config.output_dir / "DLPFC_dynamic_current.npz"
    output.parent.mkdir(parents=True, exist_ok=True)
    raw_node_names = np.empty(len(data.node_names), dtype=object)
    raw_node_names[data.raw_order] = data.node_names
    np.savez_compressed(
        output,
        time_ms=result.time_ms,
        voltage=result.voltage,
        current_abs_total=result.current_abs_total,
        current_net=result.current_net,
        current_positive_inflow=result.current_positive_inflow,
        stim_indices_zero_based=result.stim_indices,
        stim_indices_matlab=result.stim_indices + 1,
        stim_labels=np.asarray(result.stim_labels),
        node_names=np.asarray(data.node_names),
        model_to_raw_order_zero_based=data.raw_order,
        model_to_raw_order_matlab=data.raw_order + 1,
        raw_node_names=raw_node_names.astype(str),
        dense_mapping=config.dense_mapping_mode,
        model_variant=config.model_variant,
        effective_capacitance=(
            0.1 if config.model_variant == "legacy_xstruct_v8" else 1.0
        ),
        fco_threshold_applied=(config.model_variant == "main_caio"),
    )
    return output


def parse_float_list(text: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in text.split(",") if item.strip())


def parse_label_list(text: str) -> tuple[str, ...]:
    values = tuple(item.strip() for item in text.split(",") if item.strip())
    if not values:
        raise argparse.ArgumentTypeError("At least one stimulus label is required.")
    return values


def parse_metric_list(text: str) -> tuple[str, ...]:
    values = tuple(item.strip() for item in text.split(",") if item.strip())
    invalid = sorted(set(values).difference(VALID_CURRENT_METRICS))
    if not values or invalid:
        raise argparse.ArgumentTypeError(
            "Current metrics must be a comma-separated subset of: "
            + ", ".join(VALID_CURRENT_METRICS)
        )
    return tuple(dict.fromkeys(values))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--ftract-root", type=Path, default=None)
    parser.add_argument("--atlas-file", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--cortical-atlas-file", type=Path, default=None)
    parser.add_argument(
        "--model-variant",
        choices=("main_caio", "legacy_xstruct_v8"),
        default="main_caio",
        help=(
            "Use the manuscript Caio-reordered SC_data model, or reproduce the "
            "older dynamic branch with X_struct.v8, C=0.1, no F threshold, "
            "and raw 1:374 order."
        ),
    )
    parser.add_argument(
        "--dense-mapping-mode",
        choices=("matlab_sc_order", "native_tian_s1"),
        default="matlab_sc_order",
        help=(
            "91k projection: raw SC 1:360 plus 14 CIFTI structures "
            "(default), or the combined atlas-native Tian labels."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stim-label", type=parse_label_list, default=("L_p9-46v_ROI",))
    parser.add_argument("--time-scale-ms", type=float, default=200.0)
    parser.add_argument("--num-points", type=int, default=1001)
    parser.add_argument(
        "--snapshot-times",
        type=parse_float_list,
        default=(0, 10, 50, 100, 150, 200),
    )
    parser.add_argument(
        "--animation-step-ms",
        type=float,
        default=10.0,
        help="Cerebro animation interval in ms (default: 10, giving 0:10:200).",
    )
    parser.add_argument(
        "--cerebro-scale-mode",
        choices=("fixed", "framewise", "both"),
        default="both",
        help="Render fixed, per-frame adaptive, or both colorbar modes.",
    )
    parser.add_argument(
        "--cerebro-current-metrics",
        type=parse_metric_list,
        default=("total_abs",),
        help=(
            "Comma-separated current maps: total_abs, net, net_positive, "
            "positive_inflow "
            "(default: total_abs)."
        ),
    )
    parser.add_argument(
        "--visualization-relative-floor",
        type=float,
        default=0.0,
        help=(
            "For framewise maps, render frames below this fraction of the "
            "metric's global maximum as zero to avoid magnifying roundoff."
        ),
    )
    parser.add_argument("--fco-threshold", type=float, default=0.25)
    parser.add_argument("--initial-voltage", type=float, default=1.0)
    parser.add_argument("--rtol", type=float, default=1e-3)
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument(
        "--simulation-solver",
        choices=("rk45", "exact_eigh"),
        default="rk45",
        help=(
            "RK45 reproduces the MATLAB-style integration. exact_eigh evaluates "
            "the symmetric linear RC system by eigendecomposition and avoids "
            "late-time adaptive-solver noise."
        ),
    )
    parser.add_argument(
        "--no-split-multisite-voltage",
        action="store_true",
        help="Give each multi-site stimulus the full initial voltage instead of splitting it equally.",
    )
    parser.add_argument("--skip-cerebro", action="store_true")
    parser.add_argument("--make-gif", action="store_true")
    parser.add_argument("--run-all-stimulations", action="store_true")
    parser.add_argument("--run-ftract", action="store_true")
    parser.add_argument("--peak-file", type=Path, default=None)
    parser.add_argument("--peak-variable", default="I_node_max")
    return parser


def config_from_args(args: argparse.Namespace) -> Config:
    ftract_root = args.ftract_root or rcc.ftract_dir(args.data_root).parent
    cortical_atlas_file = args.cortical_atlas_file or rcc.atlas_dir(args.data_root) / (
        "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_"
        "Areas_Group_Colors.32k_fs_LR.dlabel.nii"
    )
    return Config(
        data_root=args.data_root,
        ftract_root=ftract_root,
        atlas_file=args.atlas_file,
        cortical_atlas_file=cortical_atlas_file,
        model_variant=args.model_variant,
        dense_mapping_mode=args.dense_mapping_mode,
        output_dir=args.output_dir,
        time_scale_ms=args.time_scale_ms,
        num_points=args.num_points,
        fco_threshold=args.fco_threshold,
        stim_labels=tuple(args.stim_label),
        initial_voltage=args.initial_voltage,
        simulation_solver=args.simulation_solver,
        rtol=args.rtol,
        atol=args.atol,
        split_multisite_voltage=not args.no_split_multisite_voltage,
        snapshot_times_ms=tuple(args.snapshot_times),
        animation_step_ms=args.animation_step_ms,
        cerebro_scale_mode=args.cerebro_scale_mode,
        cerebro_current_metrics=tuple(args.cerebro_current_metrics),
        visualization_relative_floor=args.visualization_relative_floor,
        make_gif=args.make_gif,
        skip_cerebro=args.skip_cerebro,
        run_all_stimulations=args.run_all_stimulations,
        run_ftract=args.run_ftract,
        peak_file=args.peak_file,
        peak_variable=args.peak_variable,
    )


def serializable_config(config: Config) -> dict[str, object]:
    values = asdict(config)
    for key, value in list(values.items()):
        if isinstance(value, Path):
            values[key] = str(value)
    if config.model_variant == "legacy_xstruct_v8":
        values["effective_model_settings"] = {
            "connectivity_file": str(rcc.xstruct_path(config.data_root)),
            "connectivity_variable": "X_struct.v8.W",
            "distance_variable": "X_struct.v8.ED",
            "capacitance": 0.1,
            "fco_threshold_applied": False,
            "node_order": "raw SC/parcel_list order (1:374)",
            "empirical_file": str(rcc.ystruct_path(config.data_root)),
            "empirical_variable": "Y_struct.v5",
            "empirical_reordered": False,
        }
    else:
        values["effective_model_settings"] = {
            "connectivity_file": str(
                config.data_root / "SC_data_neuron_2023_yihang.mat"
            ),
            "connectivity_variable": "SC",
            "distance_variable": "D",
            "capacitance": 1.0,
            "fco_threshold_applied": True,
            "node_order": "ReOrder_Caio_SC model order",
        }
    return values


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = config_from_args(args)
    if config.num_points < 2:
        raise ValueError("num_points must be at least 2.")
    if config.time_scale_ms <= 0:
        raise ValueError("time_scale_ms must be positive.")
    if config.animation_step_ms <= 0:
        raise ValueError("animation_step_ms must be positive.")
    config.output_dir.mkdir(parents=True, exist_ok=True)
    with (config.output_dir / "run_config.json").open("w", encoding="utf-8") as stream:
        json.dump(serializable_config(config), stream, indent=2)

    print("Loading and aligning RC inputs...", flush=True)
    data = load_rc_data(config)
    if config.model_variant == "legacy_xstruct_v8" and not rcc.ystruct_path(config.data_root).is_file():
        print(
            f"Y_struct not found at {rcc.ystruct_path(config.data_root)}; skipping the aligned "
            "empirical export (F-TRACT-derived file, not redistributed).",
            flush=True,
        )
    elif config.model_variant == "legacy_xstruct_v8":
        empirical = load_aligned_ystruct_v5(config, data)
        empirical_path = config.output_dir / "Y_struct_v5_aligned_empirical.npz"
        np.savez_compressed(
            empirical_path,
            amplitude=empirical.amplitude,
            probability=empirical.probability,
            n_implantations=empirical.n_implantations,
            n_stimulations=empirical.n_stimulations,
            n_responses=empirical.n_stimulations * empirical.probability,
            parcel_list=np.asarray(data.node_names[: data.n_cortex]),
        )
        print(
            "Validated no-reorder X_struct.v8/Y_struct.v5 cortical alignment.",
            flush=True,
        )
        print(f"Saved {empirical_path}", flush=True)
    stim_indices = resolve_stim_indices(data.node_names, config.stim_labels)
    for index in stim_indices:
        print(
            f"Stimulus: Python index {index}, MATLAB index {index + 1}, label {data.node_names[index]}",
            flush=True,
        )

    print("Simulating DLPFC dynamic response...", flush=True)
    dynamic = simulate_dynamic(data, config)
    dynamic_path = save_dynamic_result(dynamic, data, config)
    summary_path = config.output_dir / "DLPFC_dynamic_current_summary.png"
    plot_dynamic_summary(dynamic, data, summary_path)
    print(f"Saved {dynamic_path}", flush=True)
    print(f"Saved {summary_path}", flush=True)

    if not config.skip_cerebro:
        render_cerebro_dynamics(dynamic, data, config)

    peak_current: np.ndarray | None = None
    if config.peak_file is not None:
        peak_current = load_peak_current(config.peak_file, config.peak_variable, len(data.node_names))
    elif config.run_all_stimulations:
        peak_current = simulate_all_peak_currents(data, config)
        peak_path = config.output_dir / "I_node_max_python.npy"
        np.save(peak_path, peak_current)
        print(f"Saved {peak_path}", flush=True)

    if config.run_ftract:
        if peak_current is None:
            raise ValueError(
                "--run-ftract requires --run-all-stimulations or --peak-file."
            )
        ftract = load_ftract_data(config)
        all_pairs = threshold_sweep(
            ftract.amplitude, peak_current, data, ftract, low_fco_only=False
        )
        low_fco = threshold_sweep(
            ftract.amplitude, peak_current, data, ftract, low_fco_only=True
        )
        save_sweep_csv(all_pairs, config.output_dir / "amplitude_threshold_all_pairs.csv")
        save_sweep_csv(low_fco, config.output_dir / "amplitude_threshold_low_fco.csv")
        plot_threshold_sweep(
            all_pairs,
            config.output_dir / "amplitude_threshold_all_pairs.png",
            "Amplitude: all eligible same-hemisphere pairs",
        )
        plot_threshold_sweep(
            low_fco,
            config.output_dir / "amplitude_threshold_low_fco.png",
            "Amplitude: low-FCO same-hemisphere pairs",
        )

    print("Done.", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
