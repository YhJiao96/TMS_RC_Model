#!/usr/bin/env python3
"""Fig. 3A of the manuscript: RC model and graph comparators vs F-TRACT, with
95% bootstrap confidence intervals.

The manuscript figure was drawn by ``yihang_ftract/RC_ftract_boot_subplot.m``
from workspace variables that were never saved as a table. This script
rebuilds the same six panels (probability / amplitude x all / connected /
unconnected pairs) from the saved MATLAB matrices and exports:

* ``fig3a_bootstrap_ci.csv``         point estimates and percentile bootstrap
                                     CIs of the full and distance-partial
                                     Spearman correlations, per model;
* ``fig3a_paired_differences.csv``   paired bootstrap differences in |r|
                                     between the RC model and each comparator
                                     (same resamples for both models);
* ``fig3a_rc_transposed.csv``        sensitivity check with the RC matrix
                                     transposed ([stimulation, target]);
* ``provenance.json`` and ``README.md``.

Inputs (all pre-existing, nothing is re-simulated):

* RC matrix, one of
  - ``--rc-npz``: ``rc_peak_xstruct_w.npz`` from
    ``simulate_rc_peak_xstruct.py`` (X_struct.v8.W / ED, uniform C = 1,
    V0(stim) = 1, 200 ms / 2,000 samples). This reproduces the figure's RC
    bars to three decimals and is the version to report.
  - default: ``yihang_ftract/rc_outputs_minimal.mat`` ``I_node_max`` from
    ``rc_visualization_demo.m`` (same model but SC from
    ``SC_data_neuron_2023_yihang.mat``; 0.001-0.012 lower correlations).
  In both cases the matrix is used *as stored* (rows = recorded node,
  columns = stimulated node), exactly as the MATLAB plot script did.
* ``yihang_ftract/Precomputed_Models.mat``   ``SI``, ``SPE``, ``CMY`` and
  ``Empirical_SC`` (= ``X_struct.v8.W``) in the raw ``X_struct`` order.
* optional ``--bct-file``: a MAT file with ``Diff_eff``, ``Nav_eff``, ``SI``,
  ``SPE``, ``CMY`` recomputed with the Brain Connectivity Toolbox from
  ``X_struct.v8.W`` / ``X_struct.v8.ED`` (see ``bct_measures.m`` in the
  results folder). When given, the comparators are taken from it.
* F-TRACT-2 ``MNI-HCP-MMP1`` matrices, masks and parcel orders through
  ``rc_ftract_core.load_comparison_data`` (reproduces the figure's
  pair counts 8,916 / 4,427 / 2,285 and 1,594 / 1,157 / 185).

Run from the repository root::

    PYTHONPATH=Python python3 Python/simulate_rc_peak_xstruct.py
    PYTHONPATH=Python python3 Python/export_fig3a_bootstrap_ci.py \
        --rc-npz results/fig3a_bootstrap_ci_2026-10-28/rc_peak_xstruct_w.npz \
        --bct-file results/fig3a_bootstrap_ci_2026-10-28/bct_measures.mat
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path

import warnings

import h5py
import numpy as np
from scipy.io import loadmat
from scipy.stats import (
    ConstantInputWarning,
    NearConstantInputWarning,
    pearsonr,
    rankdata,
    spearmanr,
)

warnings.simplefilter("ignore", ConstantInputWarning)
warnings.simplefilter("ignore", NearConstantInputWarning)

import rc_ftract_core as rcc

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_DATA_ROOT = rcc.DEFAULT_DATA_ROOT

PANELS = [
    ("Probability", "All pairs"),
    ("Probability", "Connected pairs"),
    ("Probability", "Unconnected pairs"),
    ("Amplitude", "All pairs"),
    ("Amplitude", "Connected pairs"),
    ("Amplitude", "Unconnected pairs"),
]
LABELS = {
    "I_node_max": "RC",
    "SI": "SI",
    "SPE": "SPE",
    "CMY": "CMY",
    "Nav_eff": "NE",
    "Diff_eff": "DE",
    "SC": "W",
    "ED": "ED",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--bct-file", type=Path, default=rcc.comparator_dir(DEFAULT_DATA_ROOT) / "bct_measures.mat")
    parser.add_argument(
        "--rc-npz",
        type=Path,
        default=rcc.DEFAULT_OUTPUT_ROOT / "fig3a_bootstrap_ci" / "rc_peak_xstruct_w.npz",
        help=(
            "peak matrix from simulate_rc_peak_xstruct.py (374 x 374, raw "
            "X_struct order, [target, stimulus]); when omitted, I_node_max "
            "from yihang_ftract/rc_outputs_minimal.mat is used"
        ),
    )
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=rcc.DEFAULT_OUTPUT_ROOT / "fig3a_bootstrap_ci",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_rc_matrix(path: Path, order_ft: np.ndarray) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        matrix = np.array(handle["I_node_max"]).T.astype(float)
    cortex = matrix[:360, :360]
    return cortex[np.ix_(order_ft, order_ft)]


def load_comparators(
    precomputed: Path, bct_file: Path | None, order_sc: np.ndarray
) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    pre = loadmat(precomputed, squeeze_me=True)
    sources: dict[str, str] = {}
    full: dict[str, np.ndarray] = {}
    if bct_file is not None:
        bct = loadmat(bct_file, squeeze_me=True)
        for key in ("SI", "SPE", "CMY", "Nav_eff", "Diff_eff"):
            full[key] = np.asarray(bct[key], dtype=float)
            sources[key] = str(bct_file)
        sc_full = np.asarray(bct["W"], dtype=float)
        sources["SC"] = f"{bct_file} (W = X_struct.v8.W)"
        # Precomputed_Models.mat stores SI with the sign flipped and zeros on
        # the diagonal; compare |values| off the diagonal only.
        off_diagonal = ~np.eye(full["SI"].shape[0], dtype=bool)
        for key in ("SI", "SPE", "CMY"):
            a = np.abs(full[key])[off_diagonal]
            b = np.abs(np.asarray(pre[key], float))[off_diagonal]
            finite = np.isfinite(a) & np.isfinite(b)
            diff = float(np.max(np.abs(a[finite] - b[finite])))
            if diff > 1e-8:
                raise ValueError(
                    f"{key} in {bct_file} differs from Precomputed_Models.mat "
                    f"(max abs off-diagonal diff {diff:.3g})."
                )
    else:
        for key in ("SI", "SPE", "CMY"):
            full[key] = np.asarray(pre[key], dtype=float)
            sources[key] = str(precomputed)
        sc_full = np.asarray(pre["Empirical_SC"], dtype=float)
        sources["SC"] = f"{precomputed} (Empirical_SC = X_struct.v8.W)"
    full["SC"] = sc_full
    index = np.ix_(order_sc, order_sc)
    cortex = {key: value[index][:360, :360] for key, value in full.items()}
    return cortex, sources


def partial_pearson(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    design = np.column_stack([np.ones_like(z), z])
    beta_x = np.linalg.lstsq(design, x, rcond=None)[0]
    beta_y = np.linalg.lstsq(design, y, rcond=None)[0]
    rx = x - design @ beta_x
    ry = y - design @ beta_y
    if np.ptp(rx) == 0 or np.ptp(ry) == 0:
        return float("nan")
    return float(pearsonr(rx, ry).statistic)


def spearman_full_partial(
    sim: np.ndarray, emp: np.ndarray, dist: np.ndarray, with_partial: bool
) -> tuple[float, float]:
    if np.ptp(sim) == 0 or np.ptp(emp) == 0:
        return float("nan"), float("nan")
    r_sim, r_emp = rankdata(sim), rankdata(emp)
    full = float(pearsonr(r_sim, r_emp).statistic)
    partial = (
        partial_pearson(r_sim, r_emp, rankdata(dist)) if with_partial else float("nan")
    )
    return full, partial


def bootstrap_panel(
    metrics: dict[str, np.ndarray],
    emp: np.ndarray,
    dist: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Return per-metric arrays of bootstrap full / partial Spearman r.

    The same resample indices are used for every metric so that differences
    between metrics are paired.
    """
    n = emp.size
    index = rng.integers(0, n, size=(n_boot, n))
    boot_full = {key: np.empty(n_boot) for key in metrics}
    boot_part = {key: np.empty(n_boot) for key in metrics}
    for b in range(n_boot):
        idx = index[b]
        r_emp = rankdata(emp[idx])
        r_dist = rankdata(dist[idx])
        design = np.column_stack([np.ones(n), r_dist])
        pinv = np.linalg.pinv(design)
        res_emp = r_emp - design @ (pinv @ r_emp)
        for key, sim in metrics.items():
            r_sim = rankdata(sim[idx])
            if np.ptp(r_sim) == 0:
                boot_full[key][b] = np.nan
                boot_part[key][b] = np.nan
                continue
            boot_full[key][b] = pearsonr(r_sim, r_emp).statistic
            if key == "ED":
                boot_part[key][b] = np.nan
                continue
            res_sim = r_sim - design @ (pinv @ r_sim)
            if np.ptp(res_sim) == 0 or np.ptp(res_emp) == 0:
                boot_part[key][b] = np.nan
            else:
                boot_part[key][b] = pearsonr(res_sim, res_emp).statistic
    return boot_full, boot_part


def safe_std(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    return float(np.std(finite)) if finite.size else float("nan")


def ci(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan"), float("nan")
    lo, hi = np.quantile(finite, [0.025, 0.975])
    return float(lo), float(hi)


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float, digits: int = 2) -> str:
    return "n/a" if not np.isfinite(value) else f"{value:.{digits}f}"


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = rcc.load_comparison_data(args.data_root)
    order_sc = rcc.load_permutation(
        args.data_root / "ReOrder_Caio_SC.mat", "ReOrder_Caio_SC", 374
    )
    order_ft = rcc.load_permutation(
        args.data_root / "ReOreder_F_tract.mat", "ReOrder_F_tract", 360, column=1
    )
    analysis_root = args.data_root / "yihang_ftract"
    rc_file = analysis_root / "rc_outputs_minimal.mat"
    precomputed = rcc.comparator_dir(args.data_root) / "Precomputed_Models.mat"

    if args.rc_npz is None:
        rc = load_rc_matrix(rc_file, order_ft)
        rc_provenance = {
            "file": str(rc_file),
            "sha256": sha256(rc_file),
            "variable": "I_node_max[1:360,1:360]",
            "model": "rc_visualization_demo.m: SC/D from SC_data_neuron_2023_yihang.mat, R = D./SC, G = 1./R, uniform C = 1, V0(stim) = 1, ode45, 200 ms / 2,000 samples; I_node_max = max over time of the total outgoing |current| of each node",
            "parcel_order": "cortical block is in F-TRACT atlas order; permuted with ReOrder_F_tract(:,2) to match the empirical matrices",
        }
    else:
        rc_file = args.rc_npz
        peak = np.asarray(np.load(rc_file)["peak"], dtype=float)
        rc = np.round(peak[np.ix_(order_sc, order_sc)][:360, :360], 11)
        sidecar = rc_file.with_suffix(".json")
        rc_provenance = {
            "file": str(rc_file),
            "sha256": sha256(rc_file),
            "variable": "peak[:360,:360]",
            "model": "simulate_rc_peak_xstruct.py (see sidecar JSON): X_struct.v8.W / X_struct.v8.ED, uniform C = 1, V0(stim) = 1, exact eigendecomposition, 200 ms / 2,000 samples; same peak definition as rc_visualization_demo.m",
            "sidecar": json.loads(sidecar.read_text()) if sidecar.exists() else None,
            "parcel_order": "raw X_struct order permuted with ReOrder_Caio_SC to match the empirical matrices",
        }
    comparators, sources = load_comparators(precomputed, args.bct_file, order_sc)
    metrics: dict[str, np.ndarray] = {"I_node_max": rc}
    for key in ("SI", "SPE", "CMY", "Nav_eff", "Diff_eff", "SC"):
        if key in comparators:
            metrics[key] = comparators[key]
    metrics["ED"] = data.distance
    metrics_transposed = dict(metrics)
    metrics_transposed["I_node_max"] = rc.T

    rows: list[dict[str, object]] = []
    diff_rows: list[dict[str, object]] = []
    transposed_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(args.seed)
    for outcome, pair_set in PANELS:
        mask = data.masks[(outcome, pair_set)].copy()
        for sim in metrics.values():
            mask &= np.isfinite(sim)
        emp = (data.probability if outcome == "Probability" else data.amplitude)[mask]
        dist = data.distance[mask]
        n = int(mask.sum())
        panel = {key: sim[mask] for key, sim in metrics.items()}

        point: dict[str, tuple[float, float]] = {}
        for key, sim in panel.items():
            point[key] = spearman_full_partial(sim, emp, dist, key != "ED")
        boot_full, boot_part = bootstrap_panel(panel, emp, dist, args.n_boot, rng)

        for key, sim in panel.items():
            full, part = point[key]
            full_lo, full_hi = ci(boot_full[key])
            part_lo, part_hi = ci(boot_part[key])
            full_p = float(spearmanr(sim, emp).pvalue)
            rows.append(
                {
                    "outcome": outcome,
                    "pair_set": pair_set,
                    "n_pairs": n,
                    "metric": key,
                    "label": LABELS[key],
                    "r_full": full,
                    "r_full_ci_low": full_lo,
                    "r_full_ci_high": full_hi,
                    "abs_r_full": abs(full),
                    "r_full_p_asymptotic": full_p if np.isfinite(full) else np.nan,
                    "r_partial_ed": part,
                    "r_partial_ci_low": part_lo,
                    "r_partial_ci_high": part_hi,
                    "abs_r_partial": abs(part) if np.isfinite(part) else np.nan,
                    "boot_se_full": safe_std(boot_full[key]),
                    "boot_se_partial": safe_std(boot_part[key]),
                }
            )

        for other in panel:
            if other == "I_node_max":
                continue
            d_full = np.abs(boot_full["I_node_max"]) - np.abs(boot_full[other])
            d_full_lo, d_full_hi = ci(d_full)
            if other == "ED":
                d_part = np.full(args.n_boot, np.nan)
            else:
                d_part = np.abs(boot_part["I_node_max"]) - np.abs(boot_part[other])
            d_part_lo, d_part_hi = ci(d_part)
            diff_rows.append(
                {
                    "outcome": outcome,
                    "pair_set": pair_set,
                    "n_pairs": n,
                    "comparison": f"RC minus {LABELS[other]}",
                    "metric_b": other,
                    "delta_abs_r_full": abs(point["I_node_max"][0])
                    - abs(point[other][0]),
                    "delta_full_ci_low": d_full_lo,
                    "delta_full_ci_high": d_full_hi,
                    "prob_rc_higher_full": (
                        float(np.nanmean(d_full > 0))
                        if np.isfinite(d_full).any()
                        else np.nan
                    ),
                    "delta_abs_r_partial": (
                        abs(point["I_node_max"][1]) - abs(point[other][1])
                        if other != "ED"
                        else np.nan
                    ),
                    "delta_partial_ci_low": d_part_lo,
                    "delta_partial_ci_high": d_part_hi,
                    "prob_rc_higher_partial": (
                        float(np.nanmean(d_part > 0))
                        if other != "ED" and np.isfinite(d_part).any()
                        else np.nan
                    ),
                }
            )

        rc_t = metrics_transposed["I_node_max"][mask]
        full_t, part_t = spearman_full_partial(rc_t, emp, dist, True)
        transposed_rows.append(
            {
                "outcome": outcome,
                "pair_set": pair_set,
                "n_pairs": n,
                "rc_as_stored_r_full": point["I_node_max"][0],
                "rc_as_stored_r_partial": point["I_node_max"][1],
                "rc_transposed_r_full": full_t,
                "rc_transposed_r_partial": part_t,
            }
        )
        print(
            f"{outcome:11s} {pair_set:17s} n={n:5d}  "
            + "  ".join(
                f"{LABELS[k]}={abs(point[k][0]):.3f}/{abs(point[k][1]):.3f}"
                if np.isfinite(point[k][1])
                else f"{LABELS[k]}={abs(point[k][0]):.3f}/--"
                for k in panel
            )
        )

    write_csv(rows, args.output_dir / "fig3a_bootstrap_ci.csv")
    write_csv(diff_rows, args.output_dir / "fig3a_paired_differences.csv")
    write_csv(transposed_rows, args.output_dir / "fig3a_rc_transposed.csv")

    provenance = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": str(Path(__file__).resolve().relative_to(PROJECT_ROOT)),
        "n_boot": args.n_boot,
        "seed": args.seed,
        "ci_method": "percentile, 2.5 and 97.5 quantiles of the bootstrap distribution",
        "correlation": "Spearman; partial Spearman = Pearson on rank residuals after regressing on Euclidean distance (MATLAB partialcorr 'Spearman' equivalent)",
        "resampling": "pairs resampled with replacement; identical resamples for all models within a panel (paired differences)",
        "rc_matrix": {
            **rc_provenance,
            "orientation": "as stored: rows = recorded node, columns = stimulated node, correlated with F-TRACT [stimulation, target] without transposing, as in RC_ftract_boot_subplot.m",
        },
        "comparator_sources": sources,
        "comparator_definitions": {
            "SI": "search_information(W, 1./W)",
            "SPE": "1./distance_wei_floyd(1./W)",
            "CMY": "expm(S^-1/2 W S^-1/2), S = diag(strength)",
            "Nav_eff": "1./PL_wei from navigation_wu(1./W, ED)",
            "Diff_eff": "1./mean_first_passage_time(W)",
            "SC": "X_struct.v8.W",
            "ED": "X_struct.v8.ED",
        },
        "empirical": "F-tract-2/MNI-HCP-MMP1 probability, amplitude__median, N_stimulations, N_implantations via rc_matlab_comparison_leak.load_comparison_data",
        "masks": {
            "base": "same hemisphere, off-diagonal, ED >= 20 mm, N_implantations >= 5",
            "probability": "N_stimulations >= 100",
            "amplitude": "N_stimulations * probability >= 100",
            "connected": "connection frequency F >= 0.75",
            "unconnected": "F <= 0.25",
            "counts": {f"{o} / {p}": int(data.masks[(o, p)].sum()) for o, p in PANELS},
        },
        "notes": [
            "Bootstrap CIs are regenerated with a fixed seed; they are not bit-identical to the MATLAB error bars (different random stream) but agree to ~0.005.",
            "The figure's RC bars were produced from a run on X_struct.v8.W; rc_outputs_minimal.mat was run on SC_data_neuron_2023_yihang.mat (edge weights differ by up to 0.24), which lowers the RC correlations by 0.001-0.012, most in the probability/unconnected panel. Use --rc-npz with the matrix from simulate_rc_peak_xstruct.py to reproduce the figure.",
        ],
    }
    (args.output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )

    readme = [
        "# Fig. 3A correlations with 95% bootstrap CIs",
        "",
        f"Generated {provenance['generated']} by `{provenance['script']}` "
        f"(B = {args.n_boot}, seed = {args.seed}, percentile CIs).",
        "",
        "Full = Spearman r between model prediction and F-TRACT outcome; "
        "Partial = partial Spearman r controlling Euclidean distance. "
        "Values are signed (SI is a cost and ED a distance, so both correlate "
        "negatively; the figure plots |r|). RC = `I_node_max` from "
        "`rc_outputs_minimal.mat` (uniform capacitance, R = D/SC). "
        "n/a = predictor constant within the pair set (W is zero for all "
        "unconnected pairs) or not defined (partial of ED on ED).",
        "",
    ]
    for outcome, pair_set in PANELS:
        sub = [r for r in rows if r["outcome"] == outcome and r["pair_set"] == pair_set]
        sub.sort(key=lambda r: -float(np.nan_to_num(r["abs_r_full"], nan=-1.0)))
        readme += [
            f"## {outcome}, {pair_set} (n = {sub[0]['n_pairs']:,})",
            "",
            "| Model | Full r [95% CI] | Partial r [95% CI] |",
            "|---|---|---|",
        ]
        for r in sub:
            readme.append(
                f"| {r['label']} | "
                + (
                    f"{fmt(r['r_full'])} [{fmt(r['r_full_ci_low'])}, {fmt(r['r_full_ci_high'])}] | "
                    if np.isfinite(r["r_full"])
                    else "n/a | "
                )
                + (
                    f"{fmt(r['r_partial_ed'])} [{fmt(r['r_partial_ci_low'])}, "
                    f"{fmt(r['r_partial_ci_high'])}] |"
                    if np.isfinite(r["r_partial_ed"])
                    else "n/a |"
                )
            )
        readme.append("")
        dsub = [
            d for d in diff_rows if d["outcome"] == outcome and d["pair_set"] == pair_set
        ]
        readme += [
            "Paired bootstrap differences in |r| (RC minus comparator; "
            "P = share of resamples with RC higher):",
            "",
            "| Comparison | Δ full [95% CI] | P(RC > other) | Δ partial [95% CI] | P(RC > other) |",
            "|---|---|---|---|---|",
        ]
        for d in dsub:
            readme.append(
                f"| {d['comparison']} | {fmt(d['delta_abs_r_full'], 3)} "
                f"[{fmt(d['delta_full_ci_low'], 3)}, {fmt(d['delta_full_ci_high'], 3)}] | "
                f"{fmt(d['prob_rc_higher_full'])} | "
                + (
                    f"{fmt(d['delta_abs_r_partial'], 3)} "
                    f"[{fmt(d['delta_partial_ci_low'], 3)}, {fmt(d['delta_partial_ci_high'], 3)}] | "
                    f"{fmt(d['prob_rc_higher_partial'])} |"
                    if np.isfinite(d["delta_abs_r_partial"])
                    else "n/a | n/a |"
                )
            )
        readme.append("")
    readme += [
        "## RC orientation check",
        "",
        "| Outcome | Pair set | as stored full / partial | transposed full / partial |",
        "|---|---|---|---|",
    ]
    for t in transposed_rows:
        readme.append(
            f"| {t['outcome']} | {t['pair_set']} | "
            f"{fmt(t['rc_as_stored_r_full'], 3)} / {fmt(t['rc_as_stored_r_partial'], 3)} | "
            f"{fmt(t['rc_transposed_r_full'], 3)} / {fmt(t['rc_transposed_r_partial'], 3)} |"
        )
    readme.append("")
    (args.output_dir / "README.md").write_text("\n".join(readme))
    print(f"Wrote {args.output_dir}")


if __name__ == "__main__":
    main()
