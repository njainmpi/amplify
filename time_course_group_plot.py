#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot group average curves from group CSVs (expects columns: x, Group_Mean, Group_SEM)."
    )
    p.add_argument("files", nargs="+", help="Group CSV files (e.g., groupA.csv groupB.csv groupC.csv)")
    p.add_argument("--labels", nargs="*", default=None,
                   help="Custom labels for curves (same length as files). Defaults to filenames.")
    p.add_argument("--out-fig", default="groups_average_plot.png", help="Output timecourse figure filename")
    p.add_argument("--out-csv", default=None,
                   help="If set, save a combined CSV with x, means, and SEMs.")
    p.add_argument("--out-scatter", default="groups_summary_scatter.png",
                   help="Output summary scatter figure filename")
    p.add_argument("--out-summary", default=None,
                   help="If set, save a summary CSV with Group, Mean, Max.")
    p.add_argument("--font-size", type=int, default=14,
                   help="Base font size for plot (default: 14)")
    return p.parse_args()

def load_group_file(path: Path, label: str):
    """Load one group CSV; return DataFrame with x, <label>, <label>_SEM."""
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if not all(k in cols for k in ["x", "group_mean", "group_sem"]):
        raise SystemExit(f"{path} must contain columns 'x', 'Group_Mean', and 'Group_SEM'.")

    df = df[[cols["x"], cols["group_mean"], cols["group_sem"]]].copy()
    df.rename(columns={
        cols["x"]: "x",
        cols["group_mean"]: label,
        cols["group_sem"]: f"{label}_SEM"
    }, inplace=True)
    df.sort_values("x", inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def main():
    args = parse_args()
    files = [Path(f) for f in args.files]
    labels = args.labels if args.labels else [f.stem for f in files]
    if args.labels and len(args.labels) != len(files):
        raise SystemExit("Error: --labels must match number of files")

    # Load each group
    dfs = [load_group_file(f, lab) for f, lab in zip(files, labels)]

    # Merge all on x
    aligned = dfs[0]
    for d in dfs[1:]:
        aligned = aligned.merge(d, on="x", how="inner")

    # Set font sizes
    plt.rcParams.update({
        "font.size": args.font_size,
        "axes.titlesize": args.font_size + 2,
        "axes.labelsize": args.font_size + 1,
        "xtick.labelsize": args.font_size - 1,
        "ytick.labelsize": args.font_size - 1,
        "legend.fontsize": args.font_size - 1
    })

    # =============== Timecourse Plot ===============
    plt.figure(figsize=(12, 6))
    x = aligned["x"].to_numpy()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, lab in enumerate(labels):
        y = aligned[lab].to_numpy()
        sem = aligned[f"{lab}_SEM"].to_numpy()
        color = colors[i % len(colors)]

        plt.plot(x, y, linewidth=2, label=lab, color=color)
        if np.any(~np.isnan(sem)):
            plt.fill_between(x, y - sem, y + sem, color=color, alpha=0.35)

    plt.xlabel("Time (in sec)")
    plt.ylabel("Percent Signal Change")
    plt.title("Group Average Timecourses")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(args.out_fig, dpi=200)
    plt.show()
    print(f"Saved timecourse figure: {args.out_fig}")

    # Optionally save combined averages + SEMs
    if args.out_csv:
        aligned.to_csv(args.out_csv, index=False)
        print(f"Saved combined averages CSV: {args.out_csv}")

    # =============== Summary Scatter (Mean vs Max) ===============
    # For each group's column, compute overall mean and max
    summary_rows = []
    for lab in labels:
        y = aligned[lab].to_numpy()
        g_mean = float(np.nanmean(y))
        g_max  = float(np.nanmax(y))
        summary_rows.append({"Group": lab, "Mean": g_mean, "Max": g_max})

    summary_df = pd.DataFrame(summary_rows)

    # Scatter plot: one category per group, plot Mean (circle) and Max (triangle)
    plt.figure(figsize=(8, 5))
    xpos = np.arange(len(labels))
    for i, lab in enumerate(labels):
        color = colors[i % len(colors)]
        row = summary_df.iloc[i]
        # Mean
        plt.scatter(xpos[i], row["Mean"], s=80, marker="o", color=color, label=(lab if i == 0 else None))
        # Max
        plt.scatter(xpos[i], row["Max"], s=80, marker="^", facecolors="none", edgecolors=color, linewidths=2)

    # Build custom legend entries
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], marker='o', color='none', label='Group Mean', markerfacecolor='k', markersize=8),
        Line2D([0], [0], marker='^', color='k', label='Group Max', markerfacecolor='none', markersize=8)
    ]
    plt.legend(handles=legend_elems, loc="best")

    plt.xticks(xpos, labels, rotation=0)
    plt.ylabel("Percent Signal Change")
    plt.title("Group Summary: Mean vs Max of Group_Mean")
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(args.out_scatter, dpi=200)
    plt.show()
    print(f"Saved summary scatter figure: {args.out_scatter}")

    # Optional: save summary stats
    if args.out_summary:
        summary_df.to_csv(args.out_summary, index=False)
        print(f"Saved summary CSV: {args.out_summary}")

if __name__ == "__main__":
    main()
