#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import math

def load_data(filepath):
    with open(filepath, "r") as file:
        return np.array([float(line.strip()) for line in file if line.strip()])

def plot_voxels(file_list, output_path):
    num_voxels = len(file_list)
    cols = 4  # Number of columns in subplot grid
    rows = math.ceil(num_voxels / cols)

    fig, axs = plt.subplots(rows, cols, figsize=(cols * 8, rows * 5), squeeze=False)
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    for idx, f in enumerate(file_list):
        row, col = divmod(idx, cols)
        ax = axs[row][col]

        data = load_data(f)
        time = [(i + 1) / 60 for i in range(len(data))]  # TR = 1s, in minutes

        baseline = np.mean(data[79:519]) if len(data) > 550 else np.mean(data)
        normalized = ((data - baseline) / baseline) * 100

        label = os.path.splitext(os.path.basename(f))[0]
        ax.plot(time, normalized, label=label)
        ax.set_title(label, fontsize=12)
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("% Change")
        ax.axvspan(10, 20, color='gray', alpha=0.3)
        ax.grid(True)

    # Turn off unused subplots
    for idx in range(num_voxels, rows * cols):
        row, col = divmod(idx, cols)
        axs[row][col].axis('off')

    fig.suptitle("Grouped Signal Plot (Subplots per Voxel)", fontsize=24)
    fig.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_voxel_grouped.py <output_file.svg> <input_file1> [input_file2 ...]")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    plot_voxels(input_files, output_file)
