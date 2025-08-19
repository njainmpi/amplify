#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ==== CONFIG ====
WINDOW = 60  # moving-mean window (samples)
OUT_CSV = "group_avg_sem_window60.csv"
OUT_TXT = "group_avg_sem_window60.txt"
OUT_FIG = "group_avg_sem_window60.png"

def main(files):
    half = WINDOW // 2
    subjects = []

    for f in files:
        df = pd.read_csv(f, sep="\t", header=None, names=["x", "y"])
        df["y_smooth"] = df["y"].rolling(WINDOW, center=True).mean()
        # trim edges to avoid NaNs from rolling mean
        df = df.iloc[half: len(df) - half].reset_index(drop=True)
        subjects.append(df[["x", "y_smooth"]].rename(columns={"y_smooth": Path(f).stem}))

    # Align all subjects by x
    aligned = subjects[0]
    for s in subjects[1:]:
        aligned = aligned.merge(s, on="x", how="inner")

    x = aligned["x"].to_numpy()
    Y = aligned.drop(columns=["x"]).to_numpy()

    # Group stats
    mean_curve = np.mean(Y, axis=1)
    sem_curve = np.std(Y, axis=1, ddof=1) / np.sqrt(Y.shape[1])

    # Add group stats to dataframe
    aligned["Group_Mean"] = mean_curve
    aligned["Group_SEM"] = sem_curve

    # === PLOT ===
    plt.figure(figsize=(12, 6))
    for i, col in enumerate(aligned.columns[1:-2], start=1):  # individual subjects only
        plt.plot(x, aligned[col], alpha=0.35, linewidth=1, label=f"Subject {i}")
    plt.plot(x, mean_curve, linewidth=2.5, label="Group average (smoothed)")
    plt.fill_between(x, mean_curve - sem_curve, mean_curve + sem_curve,
                     alpha=0.25, label="SEM")

    plt.xlabel("X (timepoints)")
    plt.ylabel("Signal (a.u.)")
    plt.title(f"Group average ± SEM (window={WINDOW})")
    plt.legend(loc="best")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(OUT_FIG, dpi=200)
    plt.show()

    # === SAVE FILES ===
    aligned.to_csv(OUT_CSV, index=False)          # CSV for Excel
    aligned.to_csv(OUT_TXT, index=False, sep="\t") # TXT for scripts
    print(f"Saved: {OUT_CSV}, {OUT_TXT}, and {OUT_FIG}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py file1.txt file2.txt file3.txt ...")
        sys.exit(1)
    main(sys.argv[1:])
