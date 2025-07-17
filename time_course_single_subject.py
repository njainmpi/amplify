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
    

    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    plot_voxels(input_files, output_file)
    
    base_name  = os.path.splitext(output_file)[0]          # strip ".svg"
    txt_name   = f"{base_name}.txt"


    # Number of values in last 20 minutes
    values_needed = int(20 * 60 / sampling_interval)

    if total_values < values_needed:
        raise ValueError("Not enough data to cover 20 minutes.")

    last_20_min_values = values[-values_needed:]

    # Compute statistics
    mean_val = np.mean(last_20_min_values)
    perc_95_val = np.percentile(last_20_min_values, 95)

    
    # Create filenames for separate files
    mean_file = txt_name.replace('.txt', '_mean.txt')
    perc_file = txt_name.replace('.txt', '_p95.txt')

    # Save mean value
    with open(mean_file, 'w') as f_mean:
        f_mean.write(f"{mean_val:.4f}\n")

    # Save 95th percentile value
    with open(perc_file, 'w') as f_p95:
        f_p95.write(f"{perc_95_val:.4f}\n")

    print(f"Saved mean to {mean_file}")
    print(f"Saved 95th percentile to {perc_file}")