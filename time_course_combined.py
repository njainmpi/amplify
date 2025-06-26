#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse

def load_data(filepath):
    values = []
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            try:
                values.append(float(line))
            except ValueError:
                continue  # skip non-numeric lines
    return np.array(values)

def moving_average(data, window_size=5):
    return np.convolve(data, np.ones(window_size) / window_size, mode='same')

def plot_average_voxels(file_list, stimulation_range=(10, 20), title="MRI Signal Change (%)",
                        smooth=False, window_size=5, output_file=None, show_plot=True):
    if len(file_list) < 2:
        print("Please provide at least two valid text files.")
        return

    plt.figure(figsize=(14, 7))
    all_normalized = []
    time = None

    for f in file_list:
        data = load_data(f)
        time = np.arange(1, len(data) + 1) / 60.0
        baseline = np.mean(data[79:519]) if len(data) > 550 else np.mean(data)
        normalized = ((data - baseline) / baseline) * 100
        if smooth:
            normalized = moving_average(normalized, window_size)
        all_normalized.append(normalized)

        label = os.path.splitext(os.path.basename(f))[0]
        plt.plot(time, normalized, linestyle=':', alpha=0.6, label=label)

    all_normalized = np.stack(all_normalized)
    avg_normalized = np.mean(all_normalized, axis=0)
    sem_normalized = np.std(all_normalized, axis=0, ddof=1) / np.sqrt(len(file_list))

    plt.plot(time, avg_normalized, color='black', linewidth=2.5, label='Average')
    plt.fill_between(time, avg_normalized - sem_normalized, avg_normalized + sem_normalized,
                     color='red', alpha=0.2, label='±SEM')

    plt.axvspan(stimulation_range[0], stimulation_range[1], color='gray', alpha=0.3, label='Injection Window')
    plt.title(title, fontsize=18)
    plt.xlabel("Time (minutes)")
    plt.ylabel("MRI Signal Change (%)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")

    if show_plot:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot MRI signal change with optional smoothing and export.")
    parser.add_argument("files", nargs='+', help="Input .txt files")
    parser.add_argument("--smooth", action="store_true", help="Apply moving average smoothing")
    parser.add_argument("--window", type=int, default=5, help="Window size for moving average (default=5)")
    parser.add_argument("--output", type=str, help="Output filename (e.g., plot.svg, plot.png)")
    parser.add_argument("--no-show", action="store_true", help="Do not display the plot (save only)")

    args = parser.parse_args()

    input_files = [f for f in args.files if f.endswith(".txt")]
    if len(input_files) < 2:
        print("Please provide at least two valid .txt files.")
        sys.exit(1)

    plot_average_voxels(
        input_files,
        smooth=args.smooth,
        window_size=args.window,
        output_file=args.output,
        show_plot=not args.no_show
    )
