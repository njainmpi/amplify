#!/usr/bin/env python3
import sys, os, math, argparse
import numpy as np
import matplotlib.pyplot as plt

def load_data(filepath):
    with open(filepath, "r") as f:
        return np.array([float(line.strip()) for line in f if line.strip()])

def compute_baseline(x):
    # Use samples 79..518 (inclusive of 79, exclusive of 519) if available, else whole trace
    return np.mean(x[79:519]) if x.size > 550 else np.mean(x)

def plot_voxels(file_list, output_path, tr_sec=1.0, cols=4):
    num_voxels = len(file_list)
    if num_voxels == 0:
        raise ValueError("No input files provided. Pass one or more .txt files after the output path.")

    rows = math.ceil(num_voxels / cols)
    fig, axs = plt.subplots(rows, cols, figsize=(cols * 8, rows * 5), squeeze=False)
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    # For stats writing
    values_needed = int(20 * 60 / tr_sec)

    for idx, f in enumerate(file_list):
        data = load_data(f)
        time_min = (np.arange(data.size) * tr_sec) / 60.0

        baseline = compute_baseline(data)
        normalized = ((data - baseline) / baseline) * 100.0

        r, c = divmod(idx, cols)
        ax = axs[r][c]
        label = os.path.splitext(os.path.basename(f))[0]
        ax.plot(time_min, normalized, label=label)
        ax.set_title(label, fontsize=12)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("% Change")
        # Shade 10–20 min if that range exists
        if time_min.size > 0:
            shade_start, shade_end = 10.0, 20.0
            if time_min[-1] >= shade_start:
                ax.axvspan(shade_start, min(shade_end, time_min[-1]), color='gray', alpha=0.3)
        ax.grid(True)

        # ---- Per-file stats over last 20 minutes ----
        if data.size >= values_needed and values_needed > 0:
            last = data[-values_needed:]
            mean_val = float(np.mean(last))
            p95_val = float(np.percentile(last, 95))
            base_in = os.path.splitext(f)[0]
            with open(base_in + "_mean.txt", "w") as g:
                g.write(f"{mean_val:.4f}\n")
            with open(base_in + "_p95.txt", "w") as g:
                g.write(f"{p95_val:.4f}\n")
            print(f"[{label}] Saved mean → {base_in}_mean.txt; p95 → {base_in}_p95.txt")
        else:
            print(f"[{label}] Not enough samples for last-20-min stats (need {values_needed}, have {data.size}).")

    # Turn off any unused subplots
    for idx in range(num_voxels, rows * cols):
        r, c = divmod(idx, cols)
        axs[r][c].axis('off')

    fig.suptitle("Grouped Signal Plot (Subplots per Voxel)", fontsize=24)
    fig.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()
    print(f"Saved figure: {output_path}")

def main():
    ap = argparse.ArgumentParser(description="Plot normalized % change per voxel and compute last-20-min stats.")
    ap.add_argument("output_svg", help="Output SVG figure path.")
    ap.add_argument("inputs", nargs="+", help="One or more text files with one value per line.")
    ap.add_argument("--tr", dest="tr", type=float, default=1.0, help="Repetition time in seconds (default: 1.0)")
    ap.add_argument("--cols", dest="cols", type=int, default=4, help="Number of subplot columns (default: 4)")
    args = ap.parse_args()

    # Basic sanity checks
    inputs = [f for f in args.inputs if os.path.isfile(f)]
    if len(inputs) == 0:
        raise FileNotFoundError("None of the provided input files exist.")

    plot_voxels(inputs, args.output_svg, tr_sec=args.tr, cols=args.cols)

if __name__ == "__main__":
    main()
