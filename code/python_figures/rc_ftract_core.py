#!/usr/bin/env python3
"""Shared inputs, statistics and RC solver for the manuscript figure scripts.

This module collects the few functions that the Fig. 3 and Fig. 5 scripts in
this folder share:

* paths to the 374-node connectome inputs (``data/ftract_374/``), the F-TRACT
  atlas tables and the comparator-model matrices, overridable through the
  environment variables ``RC_DATA_ROOT``, ``RC_FTRACT_DIR``, ``RC_OUTPUT_ROOT``,
  ``RC_ATLAS_DIR`` and ``RC_SPHERE_DIR``;
* loaders for the MATLAB permutation vectors, the ``X_struct.v8`` connectome
  (small export ``X_struct_v8_W_ED_FCO.mat`` or the original
  ``X_struct_v1-19_FCO.mat``) and the F-TRACT ``*.txt.gz`` tables;
* the pair-selection masks of the manuscript (same hemisphere, Euclidean
  distance >= 20 mm, N_implantations >= 5, N_stimulations >= 100 for
  probability, N_responses >= 100 for amplitude; connected F >= 0.75,
  unconnected F <= 0.25);
* Spearman and distance-partial Spearman correlations;
* the linear RC solver ``simulate_rc``: ``dV/dt = -L V`` with the weighted
  graph Laplacian ``L`` of the conductance matrix, uniform capacitance
  ``C = 1`` and the initial condition ``V0(stim) = 1``, solved exactly by
  eigendecomposition. It returns, for every stimulation site, the peak over
  time of the total absolute node current ``sum_j |G_ij (V_i - V_j)|`` and the
  time of that peak (time-to-peak, TTP).

The model has no leak, delay or other extension; the conductance is
``G_ij = W_ij / ED_ij`` on existing edges as defined in the manuscript.
"""

from __future__ import annotations

import gzip
import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix
from scipy.stats import pearsonr, rankdata, spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_ROOT = Path(os.environ.get("RC_DATA_ROOT", REPO_ROOT / "data" / "ftract_374"))
DEFAULT_OUTPUT_ROOT = Path(os.environ.get("RC_OUTPUT_ROOT", REPO_ROOT / "output"))
XSTRUCT_SMALL = "X_struct_v8_W_ED_FCO.mat"
XSTRUCT_FULL = "X_struct_v1-19_FCO.mat"
FTRACT_PARCELLATION = "MNI-HCP-MMP1"


def _first_existing(candidates: list[Path]) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def ftract_dir(data_root: Path | None = None) -> Path:
    """Folder holding the F-TRACT ``MNI-HCP-MMP1`` tables (``*.txt.gz``).

    Set ``RC_FTRACT_DIR`` to the folder downloaded from https://f-tract.eu;
    otherwise ``<data_root>/F-tract-2/MNI-HCP-MMP1`` and two fallbacks are tried.
    """
    if os.environ.get("RC_FTRACT_DIR"):
        return Path(os.environ["RC_FTRACT_DIR"])
    root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    return _first_existing(
        [
            root / "F-tract-2" / FTRACT_PARCELLATION,
            root / FTRACT_PARCELLATION,
            root / "F-tract" / FTRACT_PARCELLATION,
        ]
    )


def comparator_dir(data_root: Path | None = None) -> Path:
    """Folder with ``bct_measures.mat`` and ``Precomputed_Models.mat``."""
    root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    return _first_existing([root / "comparators", root / "yihang_ftract", root])


def sphere_dir(data_root: Path | None = None) -> Path:
    """Folder with ``left_sphere_32k.mat`` / ``right_sphere_32k.mat`` (spin test)."""
    if os.environ.get("RC_SPHERE_DIR"):
        return Path(os.environ["RC_SPHERE_DIR"])
    root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    return _first_existing([REPO_ROOT / "data" / "tms_9656", root])


def atlas_dir(data_root: Path | None = None) -> Path:
    """Folder with the Glasser / Tian ``.dlabel.nii`` atlases (not redistributed)."""
    if os.environ.get("RC_ATLAS_DIR"):
        return Path(os.environ["RC_ATLAS_DIR"])
    return Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT


def xstruct_path(data_root: Path | None = None) -> Path:
    """Small ``X_struct_v8_W_ED_FCO.mat`` if present, else the full ``X_struct`` file."""
    root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    return _first_existing(
        [root / XSTRUCT_SMALL, root / "yihang_ftract" / XSTRUCT_FULL, root / XSTRUCT_FULL]
    )


YSTRUCT_FILE = "Y_struct_v1-v17_PAL_res.mat"


def ystruct_path(data_root: Path | None = None) -> Path:
    """``Y_struct_v1-v17_PAL_res.mat`` (F-TRACT-derived; not redistributed, see README)."""
    root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    return _first_existing([root / YSTRUCT_FILE, root / "yihang_ftract" / YSTRUCT_FILE])


@dataclass(frozen=True)
class XStructV8:
    """``X_struct.v8``: W (374 x 374 SC), ED (374 x 374 mm), FCO (360 x 360)."""

    W: np.ndarray
    ED: np.ndarray
    FCO: np.ndarray | None


def load_xstruct_v8(path: Path) -> XStructV8:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "W" in mat and "ED" in mat:
        fco = np.asarray(mat["FCO"], dtype=float) if "FCO" in mat else None
        return XStructV8(np.asarray(mat["W"], dtype=float), np.asarray(mat["ED"], dtype=float), fco)
    if "X_struct" not in mat or not hasattr(mat["X_struct"], "v8"):
        raise KeyError(f"Neither W/ED nor X_struct.v8 found in {path}")
    v8 = mat["X_struct"].v8
    fco = np.asarray(v8.FCO, dtype=float) if hasattr(v8, "FCO") else None
    return XStructV8(np.asarray(v8.W, dtype=float), np.asarray(v8.ED, dtype=float), fco)


@dataclass
class ComparisonData:
    """F-TRACT matrices and masks in the analysis (ReOrder_Caio_SC) order, cortex only."""

    conductance: np.ndarray
    distance: np.ndarray
    frequency: np.ndarray
    probability: np.ndarray
    amplitude: np.ndarray
    n_implantations: np.ndarray
    n_stimulations: np.ndarray
    node_names: list[str]
    masks: dict[tuple[str, str], np.ndarray]


def matlab_cellstr(values: np.ndarray) -> list[str]:
    return [str(value).strip() for value in np.asarray(values, dtype=object).reshape(-1).tolist()]


def load_permutation(path: Path, variable: str, length: int, column: int | None = None) -> np.ndarray:
    """MATLAB 1-based permutation vector -> 0-based NumPy index array."""
    values = np.asarray(loadmat(path, squeeze_me=True)[variable])
    if column is not None:
        values = values[:, column]
    values = values.reshape(-1)
    zero_based = values.astype(int) - 1
    if values.size != length or not np.array_equal(np.sort(zero_based), np.arange(length)):
        raise ValueError(f"{variable} in {path} is not a permutation of 1:{length}.")
    return zero_based


def read_ftract_gzip(path: Path, expected_size: int = 360) -> np.ndarray:
    """Read one F-TRACT ``*.txt.gz`` table (comment lines start with '#')."""
    rows: list[list[float]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(value) for value in line.split()])
    matrix = np.asarray(rows, dtype=float)
    if matrix.shape != (expected_size, expected_size):
        raise ValueError(f"{path} has shape {matrix.shape}; expected 360x360.")
    return matrix


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 3 or np.unique(x[valid]).size < 2 or np.unique(y[valid]).size < 2:
        return math.nan
    return float(spearmanr(x[valid], y[valid]).statistic)


def partial_spearman(x: np.ndarray, y: np.ndarray, distance: np.ndarray) -> float:
    """Spearman correlation of x and y after regressing the ranks on ranked distance."""
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(distance)
    x = x[valid]
    y = y[valid]
    distance = distance[valid]
    if x.size < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return math.nan
    design = np.column_stack([np.ones(x.size), rankdata(distance)])
    ranked_x = rankdata(x)
    ranked_y = rankdata(y)
    residual_x = ranked_x - design @ np.linalg.lstsq(design, ranked_x, rcond=None)[0]
    residual_y = ranked_y - design @ np.linalg.lstsq(design, ranked_y, rcond=None)[0]
    if np.ptp(residual_x) == 0 or np.ptp(residual_y) == 0:
        return math.nan
    return float(pearsonr(residual_x, residual_y).statistic)


def build_conductance(W: np.ndarray, ED: np.ndarray) -> np.ndarray:
    """``G_ij = W_ij / ED_ij`` on edges with W > 0 and ED > 0, symmetrised."""
    W = np.where(np.isfinite(W), W, 0.0)
    ED = np.where(np.isfinite(ED), ED, 0.0)
    G = np.zeros_like(W, dtype=float)
    edges = (W > 0) & (ED > 0)
    G[edges] = W[edges] / ED[edges]
    return 0.5 * (G + G.T)


EXPECTED_PAIR_COUNTS = {
    ("Probability", "All pairs"): 8916,
    ("Probability", "Connected pairs"): 4427,
    ("Probability", "Unconnected pairs"): 2285,
    ("Amplitude", "All pairs"): 1594,
    ("Amplitude", "Connected pairs"): 1157,
    ("Amplitude", "Unconnected pairs"): 185,
}


def load_comparison_data(data_root: Path | None = None, ftract_root: Path | None = None) -> ComparisonData:
    """Load the F-TRACT matrices and the manuscript's pair masks.

    Inputs under ``data_root``: ``ReOrder_Caio_SC.mat``, ``ReOreder_F_tract.mat``,
    ``SC_data_neuron_2023_yihang.mat`` (F, parcel_list) and the ``X_struct.v8``
    connectome (W, ED). F-TRACT tables come from ``ftract_root`` (default
    :func:`ftract_dir`). All matrices are returned in ReOrder_Caio_SC order,
    restricted to the 360 cortical parcels.
    """
    root = Path(data_root) if data_root is not None else DEFAULT_DATA_ROOT
    ftract = Path(ftract_root) if ftract_root is not None else ftract_dir(root)
    order_sc = load_permutation(root / "ReOrder_Caio_SC.mat", "ReOrder_Caio_SC", 374)
    order_ftract = load_permutation(root / "ReOreder_F_tract.mat", "ReOrder_F_tract", 360, column=1)
    x_struct = load_xstruct_v8(xstruct_path(root))
    structural = loadmat(root / "SC_data_neuron_2023_yihang.mat", squeeze_me=True, struct_as_record=False)

    index = np.ix_(order_sc, order_sc)
    distance_full = x_struct.ED[index]
    frequency_full = np.asarray(structural["F"], dtype=float)[index]
    conductance = build_conductance(x_struct.W, x_struct.ED)[index]
    raw_names = np.asarray(matlab_cellstr(structural["parcel_list"]), dtype=object)
    node_names = raw_names[order_sc].astype(str).tolist()

    empirical = {}
    for variable in ("N_implantations", "N_stimulations", "probability", "amplitude__median"):
        values = read_ftract_gzip(ftract / f"{variable}.txt.gz")
        empirical[variable] = values[np.ix_(order_ftract, order_ftract)]

    distance = distance_full[:360, :360]
    frequency = frequency_full[:360, :360]
    probability = empirical["probability"]
    amplitude = empirical["amplitude__median"]
    n_implantations = empirical["N_implantations"]
    n_stimulations = empirical["N_stimulations"]

    same_hemisphere = np.zeros((360, 360), dtype=bool)
    same_hemisphere[:180, :180] = True
    same_hemisphere[180:, 180:] = True
    base = (n_implantations >= 5) & (~np.eye(360, dtype=bool)) & (distance >= 20) & same_hemisphere
    outcome_masks = {
        "Probability": base & (n_stimulations >= 100) & np.isfinite(probability),
        "Amplitude": base & (n_stimulations * probability >= 100) & np.isfinite(amplitude),
    }
    masks: dict[tuple[str, str], np.ndarray] = {}
    for outcome, outcome_mask in outcome_masks.items():
        masks[(outcome, "All pairs")] = outcome_mask
        masks[(outcome, "Connected pairs")] = outcome_mask & (frequency >= 0.75)
        masks[(outcome, "Unconnected pairs")] = outcome_mask & (frequency <= 0.25)
    actual_counts = {key: int(mask.sum()) for key, mask in masks.items()}
    if actual_counts != EXPECTED_PAIR_COUNTS:
        raise ValueError(f"Pair masks do not reproduce the manuscript counts: {actual_counts}")

    return ComparisonData(
        conductance=conductance,
        distance=distance,
        frequency=frequency,
        probability=probability,
        amplitude=amplitude,
        n_implantations=n_implantations,
        n_stimulations=n_stimulations,
        node_names=node_names,
        masks=masks,
    )


def simulate_rc(
    conductance: np.ndarray,
    duration_ms: float,
    num_points: int,
    n_stimulations: int,
    chunk_size: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Linear RC network, uniform C = 1, V0(stim) = 1, exact eigendecomposition.

    Returns ``(peak, ttp_ms, time_ms)`` with ``peak[target, stimulus]`` the
    maximum over time of the total absolute node current
    ``sum_j |G_ij (V_i - V_j)|`` and ``ttp_ms`` the time of that maximum, for the
    first ``n_stimulations`` nodes stimulated one at a time.
    """
    G = np.asarray(conductance, dtype=float)
    if not np.allclose(G, G.T, atol=1e-12):
        raise ValueError("Exact eigendecomposition requires a symmetric conductance matrix.")
    n = G.shape[0]
    time_ms = np.linspace(0.0, duration_ms, num_points)
    laplacian = np.diag(G.sum(axis=1)) - G
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    eigenvalues = np.maximum(eigenvalues, 0.0)
    modal_decay = np.exp(-np.outer(time_ms, eigenvalues))

    edge_i, edge_j = np.where(np.triu(G, 1) > 0)
    edge_g = G[edge_i, edge_j]
    edge_count = edge_g.size
    incidence = coo_matrix(
        (
            np.ones(2 * edge_count),
            (np.concatenate([edge_i, edge_j]), np.concatenate([np.arange(edge_count), np.arange(edge_count)])),
        ),
        shape=(n, edge_count),
    ).tocsr()

    peaks = np.zeros((n, n_stimulations), dtype=np.float64)
    ttp = np.zeros_like(peaks)
    for stimulus in range(n_stimulations):
        voltage = (modal_decay * eigenvectors[stimulus, :][None, :]) @ eigenvectors.T
        for start in range(0, num_points, chunk_size):
            stop = min(start + chunk_size, num_points)
            chunk_voltage = voltage[start:stop]
            branch = np.abs(chunk_voltage[:, edge_i] - chunk_voltage[:, edge_j]) * edge_g[None, :]
            current = np.asarray(branch @ incidence.T)
            local_index = np.argmax(current, axis=0)
            local_peak = current[local_index, np.arange(n)]
            update = local_peak > peaks[:, stimulus]
            peaks[update, stimulus] = local_peak[update]
            ttp[update, stimulus] = time_ms[start + local_index[update]]
        if stimulus == 0 or (stimulus + 1) % 40 == 0:
            print(f"Computed stimulus {stimulus + 1}/{n_stimulations}", flush=True)
    return peaks, ttp, time_ms
