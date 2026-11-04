#!/usr/bin/env python3
"""Recompute the no-leak RC peak-current matrix of the manuscript's F-TRACT
analysis from ``X_struct.v8.W`` / ``X_struct.v8.ED`` and save it.

Model (identical to ``rc_visualization_demo.m`` / ``RC_FTRACT_ODE_V23.m``):
``G = SC / D`` on existing edges, uniform capacitance ``C = 1``,
``dV/dt = -L V`` with ``V0(stim) = 1``, 200 ms, ``num_points`` samples;
``peak[target, stim]`` = maximum over time of the total outgoing |current|
of ``target``. The integration uses the exact eigendecomposition of
``rc_ftract_core.simulate_rc`` (validated against the
MATLAB ode45 run to 1e-13).

Why this exists: the manuscript's Fig. 3A RC bars come from a run on
``X_struct.v8.W`` (the same connectome as the graph comparators), while the
only saved RC matrix (``yihang_ftract/rc_outputs_minimal.mat``) was run on
``SC_data_neuron_2023_yihang.mat``, whose edge weights differ by up to 0.24.
This script restores the figure's version.

Output: ``<output-dir>/rc_peak_xstruct_w.npz`` with ``peak`` and ``ttp_ms``
(374 x 374, raw ``X_struct`` parcel order, [target, stimulus]) and a JSON
sidecar with provenance.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

import rc_ftract_core as rcc

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_X_FILE = rcc.xstruct_path(rcc.DEFAULT_DATA_ROOT)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--x-file", type=Path, default=DEFAULT_X_FILE)
    parser.add_argument("--duration-ms", type=float, default=200.0)
    parser.add_argument("--num-points", type=int, default=2000)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "output" / "fig3a_bootstrap_ci",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    x = rcc.load_xstruct_v8(args.x_file)
    W = np.array(x.W, dtype=float)
    ED = np.array(x.ED, dtype=float)
    W[~np.isfinite(W)] = 0.0
    ED[~np.isfinite(ED)] = 0.0
    G = np.zeros_like(W)
    edges = (W > 0) & (ED > 0)
    G[edges] = W[edges] / ED[edges]
    G = 0.5 * (G + G.T)

    peaks, ttp, _ = rcc.simulate_rc(G, args.duration_ms, args.num_points, G.shape[0])
    out = args.output_dir / "rc_peak_xstruct_w.npz"
    np.savez_compressed(out, peak=peaks, ttp_ms=ttp)
    sidecar = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": "code/python_figures/simulate_rc_peak_xstruct.py",
        "x_file": str(args.x_file),
        "x_file_sha256": sha256(args.x_file),
        "connectivity": "X_struct.v8.W",
        "distance": "X_struct.v8.ED",
        "conductance": "G_ij = W_ij / ED_ij on edges with W > 0 and ED > 0",
        "capacitance": "uniform, C = 1",
        "initial_condition": "V0(stim) = 1, no external current",
        "duration_ms": args.duration_ms,
        "num_points": args.num_points,
        "peak_definition": "max over time of sum_j |G_ij (V_i - V_j)| for node i",
        "orientation": "[target, stimulus], raw X_struct parcel order (374 nodes)",
        "edges": int(np.count_nonzero(np.triu(G, 1))),
    }
    (args.output_dir / "rc_peak_xstruct_w.json").write_text(
        json.dumps(sidecar, indent=2) + "\n"
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
