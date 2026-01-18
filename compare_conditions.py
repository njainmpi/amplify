#!/usr/bin/env python3
import argparse
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ============================================================
# Argument parsing
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Compare seed-based voxelwise connectivity across conditions"
    )
    p.add_argument("--before", required=True, help="Before condition correlation file")
    p.add_argument("--during", required=True, help="During condition correlation file")
    p.add_argument("--after",  required=True, help="After condition correlation file")
    p.add_argument("--output-prefix", default="comparison", help="Output prefix")
    return p.parse_args()

# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------
    before = np.loadtxt(args.before)
    during = np.loadtxt(args.during)
    after  = np.loadtxt(args.after)

    # Columns: x y z r p z_fisher
    z_before = before[:, 5]
    z_during = during[:, 5]
    z_after  = after[:, 5]

    assert z_before.shape == z_during.shape == z_after.shape, "Voxel mismatch"

    voxel_z = before[:, 2].astype(int)  # z-slice index

    # --------------------------------------------------------
    # Global statistics (mask-averaged)
    # --------------------------------------------------------
    mean_before = np.mean(z_before)
    mean_during = np.mean(z_during)
    mean_after  = np.mean(z_after)

    t_db, p_db = stats.ttest_rel(z_during, z_before)
    t_ab, p_ab = stats.ttest_rel(z_after,  z_before)
    t_ad, p_ad = stats.ttest_rel(z_after,  z_during)

    # Save statistics
    with open(f"{args.output_prefix}_global_stats.txt", "w") as f:
        f.write("Condition\tMean_FisherZ\n")
        f.write(f"Before\t{mean_before:.6f}\n")
        f.write(f"During\t{mean_during:.6f}\n")
        f.write(f"After\t{mean_after:.6f}\n\n")

        f.write("Comparison\tT\tP\n")
        f.write(f"During_vs_Before\t{t_db:.4f}\t{p_db:.4e}\n")
        f.write(f"After_vs_Before\t{t_ab:.4f}\t{p_ab:.4e}\n")
        f.write(f"After_vs_During\t{t_ad:.4f}\t{p_ad:.4e}\n")

    print("[OK] Global statistics saved")

    # --------------------------------------------------------
    # Difference maps
    # --------------------------------------------------------
    dz_db = z_during - z_before
    dz_ab = z_after  - z_before
    dz_ad = z_after  - z_during

    # --------------------------------------------------------
    # Publication-ready plots
    # --------------------------------------------------------
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14
    })

    # =========================
    # FIGURE 1: Global paired bar plot
    # =========================
    fig1, ax = plt.subplots(figsize=(8, 6))

    means = [mean_before, mean_during, mean_after]
    labels = ["Before", "During", "After"]
    x = np.arange(3)

    ax.bar(x, means, color=["gray", "orange", "green"], alpha=0.8)

    for zb, zd, za in zip(z_before, z_during, z_after):
        ax.plot(x, [zb, zd, za], color="black", alpha=0.15)

    ax.set_ylabel("Mean Fisher z")
    ax.set_title("Seed-based Functional Connectivity")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    fig1.tight_layout()
    fig1.savefig(f"{args.output_prefix}_global_connectivity.png", dpi=300)
    plt.close(fig1)

    # =========================
    # FIGURE 2: Δz distributions
    # =========================
    fig2, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    diffs = [dz_db, dz_ab, dz_ad]
    titles = ["During − Before", "After − Before", "After − During"]

    for ax, d, t in zip(axes, diffs, titles):
        ax.hist(d, bins=80, color="steelblue", alpha=0.85)
        ax.axvline(0, color="black", lw=1)
        ax.set_title(t)
        ax.set_xlabel("Δ Fisher z")

    axes[0].set_ylabel("Voxel count")

    fig2.tight_layout()
    fig2.savefig(f"{args.output_prefix}_deltaZ_distributions.png", dpi=300)
    plt.close(fig2)

    # =========================
    # FIGURE 3: Slice-wise Δz scatter
    # =========================
    unique_slices = np.unique(voxel_z)[:16]

    fig3 = plt.figure(figsize=(26, 14))
    gs = GridSpec(3, 8, hspace=0.4)

    for idx, sl in enumerate(unique_slices):
        row = idx // 8
        col = idx % 8
        ax = fig3.add_subplot(gs[row, col])

        mask = voxel_z == sl
        ax.scatter(
            np.arange(np.sum(mask)),
            dz_db[mask],
            s=14,
            color="crimson"
        )

        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(f"Slice z={sl}", fontsize=11)
        ax.set_ylim(-1, 1)
        ax.tick_params(labelbottom=False)

    fig3.suptitle("Voxel-wise Δ Fisher z (During − Before)", fontsize=18)
    fig3.tight_layout(rect=[0, 0, 1, 0.95])
    fig3.savefig(f"{args.output_prefix}_slicewise_deltaZ.png", dpi=300)
    plt.close(fig3)

    print("[OK] All comparison figures saved")

# ============================================================
if __name__ == "__main__":
    main()
