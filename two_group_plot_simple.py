import sys
import numpy as np
import matplotlib.pyplot as plt

# Check input
if len(sys.argv) < 3 or (len(sys.argv) - 1) % 2 != 0:
    print("Usage: python plot.py mean1.txt sem1.txt [mean2.txt sem2.txt ...]")
    sys.exit(1)

# Load all input pairs
num_groups = (len(sys.argv) - 1) // 2
data_groups = []
# Custom labels for each group
labels = ['FAP_with_FAP_AAV', 'PBS_with_FAP_AAV', 'PBS_with_Control_AAV']

# Load all input pairs with custom labels
for i in range(num_groups):
    mean_file = sys.argv[1 + 2 * i]
    sem_file = sys.argv[2 + 2 * i]
    mean = np.loadtxt(mean_file)
    sem = np.loadtxt(sem_file)
    label = labels[i] if i < len(labels) else f"Group_{i+1}"
    data_groups.append((mean, sem, label))
# X-axis setup
time = np.arange(1, len(data_groups[0][0]) + 1) / 60.0
stimulation_range = (10, 20)
title = "MRI Signal Change (%)"

# Plot
plt.figure(figsize=(10, 6))

for mean, sem, label in data_groups:
    plt.plot(time, mean, label=label, linewidth=2)
    plt.fill_between(time, mean - sem, mean + sem, alpha=0.3)

# Injection window
plt.axvspan(stimulation_range[0], stimulation_range[1], color='gray', alpha=0.3, label='Injection Window')

# Labels and formatting
plt.title(title, fontsize=18)
plt.xlabel("Time (minutes)", fontsize=20)
plt.ylabel("MRI Signal Change (%)", fontsize=20)
plt.grid(True)
plt.legend()
plt.xlim(3, 37)
plt.ylim(-10, 40)

# Clean up axes
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
