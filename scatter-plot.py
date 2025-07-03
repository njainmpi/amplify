import numpy as np
import matplotlib.pyplot as plt

# Data matrix: rows = groups, columns = subjects
data_matrix = np.array([
    [0.0303616067864270,  0.0271333932135730],
    [0.0146175528942116, -0.0354095409181637],
    [-0.0763563013972056, -0.0419760518962076]
])

# Transpose so each subject is a row
data_matrix_T = data_matrix.T

# Group labels and subject labels
group_labels = ['FAP_with_FAP_AAV', 'PBS_with_FAP_AAV', 'PBS_with_Control_AAV']
subject_labels = ['Subject 1', 'Subject 2']

# Colors (you can customize)
colors = {
    'Group 1': '#1f77b4',
    'Group 2': '#ff7f0e',
    'Group 3': '#2ca02c'
}
# Create figure
fig, ax = plt.subplots(figsize=(6, 5))

# Plot each subject across all groups
for i, subject_data in enumerate(data_matrix_T):
    ax.scatter([1, 2, 3], subject_data, s=80, color=colors[i], label=subject_labels[i])

# Formatting
ax.set_xlim(0.5, 3.5)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(group_labels, rotation=15)
ax.set_ylabel('Signal Change (%)')
ax.set_title("Pre Injection")
ax.grid(True)
ax.legend(loc='best')

plt.tight_layout()
plt.show()
