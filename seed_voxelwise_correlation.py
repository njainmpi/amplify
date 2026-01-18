#!/usr/bin/env python3
import argparse
import numpy as np
import nibabel as nib
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ============================================================
# Argument parser
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Seed-based voxelwise correlation from 4D functional data"
    )
    p.add_argument("--func", required=True, help="4D functional NIfTI (x,y,z,t)")
    p.add_argument("--seed-mask", required=True, help="Seed ROI mask NIfTI")
    p.add_argument("--target-mask", required=True, help="Target mask NIfTI")
    p.add_argument("--start", type=int, required=True, help="Start index (inclusive)")
    p.add_argument("--stop",  type=int, required=True, help="Stop index (exclusive)")
    p.add_argument("--output-prefix", default="seed_corr", help="Output prefix")
    return p.parse_args()

# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    start, stop = args.start, args.stop

    # --------------------------------------------------------
    # Load functional data
    # --------------------------------------------------------
    func_img = nib.load(args.func)
    func = func_img.get_fdata()  # x,y,z,t

    if func.ndim != 4:
        raise ValueError("Functional image must be 4D")

    T = func.shape[3]

    if not (0 <= start < stop <= T):
        raise ValueError("Invalid start/stop indices")

    # --------------------------------------------------------
    # Load masks
    # --------------------------------------------------------
    seed_mask = nib.load(args.seed_mask).get_fdata().astype(bool)
    target_mask = nib.load(args.target_mask).get_fdata().astype(bool)

    if seed_mask.shape != func.shape[:3]:
        raise ValueError("Seed mask shape mismatch")

    if target_mask.shape != func.shape[:3]:
        raise ValueError("Target mask shape mismatch")

    # --------------------------------------------------------
    # Extract voxel indices
    # --------------------------------------------------------
    seed_vox_idx = np.column_stack(np.where(seed_mask))
    target_vox_idx = np.column_stack(np.where(target_mask))

    np.savetxt(
        f"{args.output_prefix}_seed_voxel_indices.txt",
        seed_vox_idx, fmt="%d"
    )
    np.savetxt(
        f"{args.output_prefix}_target_voxel_indices.txt",
        target_vox_idx, fmt="%d"
    )

    # --------------------------------------------------------
    # Extract time series
    # --------------------------------------------------------
    ts_seed = func[seed_mask].T      # shape: T × Nseed
    ts_target = func[target_mask].T  # shape: T × Nvox

    # --------------------------------------------------------
    # Mean seed time series
    # --------------------------------------------------------
    mean_ts_seed = ts_seed.mean(axis=1)
    np.savetxt(
        f"{args.output_prefix}_mean_seed_timeseries.txt",
        mean_ts_seed, fmt="%.6f"
    )

    x_ts = mean_ts_seed[start:stop]

    # --------------------------------------------------------
    # Voxelwise correlation
    # --------------------------------------------------------
    n_vox = ts_target.shape[1]
    rvals = np.zeros(n_vox, dtype=np.float32)
    pvals = np.zeros(n_vox, dtype=np.float32)

    for i in range(n_vox):
        y_ts = ts_target[start:stop, i]
        r, p = scipy.stats.pearsonr(x_ts, y_ts)
        rvals[i] = r
        pvals[i] = p

    # --------------------------------------------------------
    # Fisher z-transform
    # --------------------------------------------------------
    eps = np.finfo(np.float32).eps
    rvals_clip = np.clip(rvals, -1 + eps, 1 - eps)
    zvals = np.arctanh(rvals_clip)

    # --------------------------------------------------------
    # Save combined output
    # --------------------------------------------------------
    out = np.column_stack((target_vox_idx, rvals, pvals, zvals))

    np.savetxt(
        f"{args.output_prefix}_voxelwise_correlation.txt",
        out,
        fmt="%d\t%d\t%d\t%.6f\t%.6e\t%.6f",
        header="x\ty\tz\tr\tp\tz_fisher"
    )

    print("Rows written:", out.shape[0])
    print("Mask voxels:", np.sum(target_mask))

   
# ============================================================
# -------------------- PLOTTING (PNG, COMBINED) ---------------
# ============================================================

    group_ids = out[:, 2].astype(int)   # z-slice indices
    y_all = rvals
    x_vox = np.arange(len(y_all))
    unique_groups = np.unique(group_ids)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 13
    })

    # ------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------
    n_slices_max = 16
    n_groups = min(len(unique_groups), n_slices_max)

    fig = plt.figure(figsize=(26, 16))
    gs = GridSpec(
        nrows=1 + int(np.ceil(n_groups / 8)),
        ncols=8,
        height_ratios=[2.8] + [1.6] * int(np.ceil(n_groups / 8)),
        hspace=0.5
    )

    # ------------------------------------------------------------
    # Helper: scatter by correlation range
    # ------------------------------------------------------------
    def scatter_by_range(ax, x, y):
        ax.scatter(x[y == 0], y[y == 0], color="black", s=14, label="0")
        ax.scatter(x[y > 0.5], y[y > 0.5], color="red", s=14, label="0.5–1.0")
        ax.scatter(x[(y > 0) & (y <= 0.5)], y[(y > 0) & (y <= 0.5)],
                color="orange", s=14, label="0–0.5")
        ax.scatter(x[(y >= -0.5) & (y < 0)], y[(y >= -0.5) & (y < 0)],
                color="greenyellow", s=14, label="-0.5–0")
        ax.scatter(x[y < -0.5], y[y < -0.5], color="blue", s=14, label="-1––0.5")

        ax.axhline(0, color="gray", lw=0.8, alpha=0.6)

    # ------------------------------------------------------------
    # TOP PANEL — All voxels
    # ------------------------------------------------------------
    ax_top = fig.add_subplot(gs[0, :])
    scatter_by_range(ax_top, x_vox, y_all)

    ax_top.set_title("Voxel-wise Correlation Coefficient — All Mask Voxels", fontsize=16)
    ax_top.set_ylabel("Pearson r")
    ax_top.set_xlabel("Voxel index")

    handles, labels = ax_top.get_legend_handles_labels()
    ax_top.legend(handles, labels, ncol=5, frameon=False, loc="upper right")

    # ------------------------------------------------------------
    # SLICE-WISE PANELS
    # ------------------------------------------------------------
    for idx, g in enumerate(unique_groups[:n_slices_max]):
        row = 1 + idx // 8
        col = idx % 8

        ax = fig.add_subplot(gs[row, col], sharex=ax_top, sharey=ax_top)

        mask = group_ids == g
        scatter_by_range(ax, x_vox[mask], y_all[mask])

        ax.set_title(f"Slice z = {g}", fontsize=11)
        ax.tick_params(axis="x", labelbottom=False)
        ax.tick_params(axis="y", labelsize=9)

    # ------------------------------------------------------------
    # Finalize + save
    # ------------------------------------------------------------
    fig.tight_layout()
    out_png = f"{args.output_prefix}_corr_coeff_all_slices.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    print(f"[OK] Combined slice-wise PNG saved → {out_png}")



if __name__ == "__main__":
    main()
