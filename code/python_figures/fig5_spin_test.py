#!/usr/bin/env python3
"""Spatial-permutation (spin) null for the Fig. 5 delay correlations.

Parcel-level spin (Vasa et al. 2018): Glasser parcel centroids on the fs_LR 32k
sphere are rotated by a random rotation (mirrored for the right hemisphere) and
matched one-to-one to the original centroids (Hungarian assignment).  Each spin
is a permutation of the 360 cortical nodes; the model matrix (RC time-to-peak,
or a comparator) is permuted on rows and columns together, the empirical latency,
Euclidean distance and the pair mask stay fixed, and Spearman / ED-partial
Spearman correlations are recomputed.  p = (1 + #null >= observed) / (1 + n_perm).

Outputs (``<output-dir>``): ``spin_perms.npy`` (n_perm x 360, in the reordered
node space), ``spin_nulls.npz`` (null distributions), ``spin_results.csv``, ``spin_null.png``.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402
from scipy.optimize import linear_sum_assignment  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rc_ftract_core as rcc  # noqa: E402
from export_fig3a_bootstrap_ci import load_comparators  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # repository root
DEFAULT_DATA_ROOT = rcc.DEFAULT_DATA_ROOT
FIG3_DIR = PROJECT_ROOT / "output" / "fig3a_bootstrap_ci"
DLABEL = "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii"
RC_RED = (0.8392, 0.1529, 0.1569)


def parcel_sphere_centroids(root: Path) -> tuple[np.ndarray, list[str]]:
    """Mean fs_LR-32k sphere coordinate of every Glasser parcel, keys 1..360 -> rows 0..359."""
    img = nib.load(rcc.atlas_dir(root) / DLABEL)
    keys = np.asarray(img.dataobj)[0].astype(int)
    names = {k: v[0] for k, v in img.header.get_axis(0).label[0].items()}
    coords = np.zeros((keys.size, 3))
    for name, sl, model in img.header.get_axis(1).iter_structures():
        verts = loadmat(rcc.sphere_dir(root) / ("left_sphere_32k.mat" if "LEFT" in name else "right_sphere_32k.mat"))
        xyz = np.asarray(verts["verticesl" if "LEFT" in name else "verticesr"], float)
        coords[sl] = xyz[model.vertex]
    cent = np.zeros((360, 3))
    for k in range(1, 361):
        c = coords[keys == k].mean(0)
        cent[k - 1] = 100.0 * c / np.linalg.norm(c)
    return cent, [names[k] for k in range(1, 361)]


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def spin_permutations(cent: np.ndarray, is_left: np.ndarray, n_perm: int, rng: np.random.Generator) -> np.ndarray:
    """n_perm x 360 permutations (perm[k] = original parcel whose value node k receives)."""
    mirror = np.diag([-1.0, 1.0, 1.0])
    left, right = np.flatnonzero(is_left), np.flatnonzero(~is_left)
    perms = np.zeros((n_perm, 360), dtype=int)
    for p in range(n_perm):
        rot = random_rotation(rng)
        for idx, R in ((left, rot), (right, mirror @ rot @ mirror)):
            rotated = cent[idx] @ R.T
            cost = np.linalg.norm(rotated[:, None, :] - cent[idx][None, :, :], axis=2)
            rows, cols = linear_sum_assignment(cost)
            perms[p, idx[rows]] = idx[cols]
    return perms


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    return rcc.partial_spearman(x, y, z)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--rc-npz", type=Path, default=FIG3_DIR / "rc_peak_xstruct_w.npz")
    parser.add_argument("--bct-file", type=Path, default=rcc.comparator_dir(DEFAULT_DATA_ROOT) / "bct_measures.mat")
    parser.add_argument("--delay-variable", default="onset_delay", choices=("peak_delay", "onset_delay", "latency_start"))
    parser.add_argument("--count-variable", default="stimulations", choices=("responses", "stimulations"))
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # --- data in the reordered (F-TRACT-aligned) node space
    data = rcc.load_comparison_data(args.data_root)
    order_sc = rcc.load_permutation(args.data_root / "ReOrder_Caio_SC.mat", "ReOrder_Caio_SC", 374)
    order_ft = rcc.load_permutation(args.data_root / "ReOreder_F_tract.mat", "ReOrder_F_tract", 360, column=1)
    ixs = np.ix_(order_sc, order_sc)
    delay = rcc.read_ftract_gzip(rcc.ftract_dir(args.data_root) / f"{args.delay_variable}__median.txt.gz")
    delay = delay[np.ix_(order_ft, order_ft)]
    with np.load(args.rc_npz) as npz:
        ttp = npz["ttp_ms"][ixs][:360, :360].astype(float)
    ttp[ttp >= 200] = np.nan
    comp, _ = load_comparators(rcc.comparator_dir(args.data_root) / "Precomputed_Models.mat", args.bct_file, order_sc)
    spl = np.where(comp["SPE"] > 0, 1.0 / np.where(comp["SPE"] > 0, comp["SPE"], 1.0), np.nan)
    predictors = {"RC": ttp, "SPL": spl, "SI": np.abs(comp["SI"])}

    # --- parcel centroids in the same node space; raw SC order == dlabel key order (checked by name)
    cent_raw, names_raw = parcel_sphere_centroids(args.data_root)
    parcel_list = [str(s) for s in loadmat(args.data_root / "SC_data_neuron_2023_yihang.mat", squeeze_me=True)["parcel_list"][:360]]
    assert parcel_list == names_raw, "raw SC parcel order differs from the dlabel key order"
    cent = cent_raw[order_sc[:360]]
    names = [names_raw[i] for i in order_sc[:360]]
    is_left = np.array([n.startswith("L_") for n in names])
    assert is_left[:180].all() == is_left[:180].any() and is_left[180:].all() == is_left[180:].any(), "hemisphere blocks"
    perms = spin_permutations(cent, is_left, args.n_perm, rng)
    np.save(args.output_dir / "spin_perms.npy", perms)
    moved = np.mean(perms != np.arange(360))
    print(f"{args.n_perm} spins generated (seed {args.seed}); fraction of nodes relabelled per spin: {moved:.3f}")

    # --- fixed pair mask (Fig. 3 criteria + delay criteria, self-pairs excluded)
    count = data.n_stimulations if args.count_variable == "stimulations" else data.n_stimulations * data.probability
    hemi = np.zeros((360, 360), bool)
    hemi[:180, :180] = True
    hemi[180:, 180:] = True
    base = (np.isfinite(delay) & (delay > 0) & (delay <= 200) & (count >= 100) & ~np.eye(360, dtype=bool)
            & hemi & (data.distance >= 20) & (data.n_implantations >= 5))
    subsets = {"All pairs": base, "Connected pairs": base & (data.frequency >= 0.75),
               "Unconnected pairs": base & (data.frequency <= 0.25)}

    rows, nulls = [], {}
    for metric, mat in predictors.items():
        for pair_set, mask in subsets.items():
            m = mask & np.isfinite(mat)
            obs_full = spearmanr(mat[m], delay[m]).statistic
            obs_part = partial_spearman(mat[m], delay[m], data.distance[m])
            null_full, null_part = np.empty(args.n_perm), np.empty(args.n_perm)
            for p in range(args.n_perm):
                mp = mat[np.ix_(perms[p], perms[p])]
                mm = mask & np.isfinite(mp)
                null_full[p] = spearmanr(mp[mm], delay[mm]).statistic
                null_part[p] = partial_spearman(mp[mm], delay[mm], data.distance[mm])
            p_full = (1 + np.sum(null_full >= obs_full)) / (1 + args.n_perm)
            p_part = (1 + np.sum(null_part >= obs_part)) / (1 + args.n_perm)
            p_part_2s = (1 + np.sum(np.abs(null_part) >= abs(obs_part))) / (1 + args.n_perm)
            nulls[(metric, pair_set)] = (null_full, null_part, obs_full, obs_part)
            rows.append({"metric": metric, "pair_set": pair_set, "n_pairs": int(m.sum()), "rho": obs_full,
                         "null_rho_mean": null_full.mean(), "null_rho_q95": np.percentile(null_full, 95), "p_spin_rho": p_full,
                         "partial": obs_part, "null_partial_mean": null_part.mean(), "null_partial_q95": np.percentile(null_part, 95),
                         "null_partial_q2.5": np.percentile(null_part, 2.5), "null_partial_q97.5": np.percentile(null_part, 97.5),
                         "p_spin_partial": p_part, "p_spin_partial_two_sided": p_part_2s})
            print(f"{metric:4s} {pair_set:18s} n={m.sum():5d}  rho={obs_full:.3f} null mean {null_full.mean():.3f} q95 {np.percentile(null_full, 95):.3f} p={p_full:.4f}  |  "
                  f"partial={obs_part:.3f} null mean {null_part.mean():.3f} [2.5–97.5%: {np.percentile(null_part, 2.5):.3f}, {np.percentile(null_part, 97.5):.3f}] p={p_part:.4f}", flush=True)
    np.savez_compressed(args.output_dir / "spin_nulls.npz",
                        **{f"{k[0]}__{k[1].replace(' ', '_')}__full": v[0] for k, v in nulls.items()},
                        **{f"{k[0]}__{k[1].replace(' ', '_')}__partial": v[1] for k, v in nulls.items()})
    with (args.output_dir / "spin_results.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for c, pair_set in enumerate(subsets):
        for r, (stat, label) in enumerate(((2, "Spearman ρ"), (3, "partial ρ | ED"))):
            ax = axes[r, c]
            null_full, null_part, obs_full, obs_part = nulls[("RC", pair_set)]
            null, obs = (null_full, obs_full) if stat == 2 else (null_part, obs_part)
            ax.hist(null, bins=40, color="0.7", edgecolor="white")
            ax.axvline(obs, color=RC_RED, linewidth=2.2)
            pv = [x for x in rows if x["metric"] == "RC" and x["pair_set"] == pair_set][0]
            ax.set_title(f"{pair_set}: {label}\nobserved {obs:.3f}, p_spin = {pv['p_spin_rho' if stat == 2 else 'p_spin_partial']:.3f}", fontsize=11)
            ax.set_xlabel(label)
            ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"RC time-to-peak vs F-TRACT {args.delay_variable}: spin null ({args.n_perm} parcel-level spins)", fontsize=13)
    fig.tight_layout()
    fig.savefig(args.output_dir / "spin_null.png", dpi=200, facecolor="white")
    print(f"Wrote {args.output_dir / 'spin_results.csv'} and spin_null.png")


if __name__ == "__main__":
    main()
