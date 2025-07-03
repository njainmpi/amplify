import numpy as np
import matplotlib.pyplot as plt

# Data matrix (3 groups x 2 subjects)
data_matrix = np.array([
    [25.8667584859219,	8.00011681834696],
    [0.776932118074478,	3.81829819436876],
    [-5.00025374296094,	-25.1086727847411]  # PBS_with_Control_AAV
])




# Transpose so each subject's data is in a row
data_matrix_T = data_matrix.T

# Group and subject labels
group_labels = ['FAP_with_FAP_AAV', 'PBS_with_FAP_AAV', 'PBS_with_Control_AAV']
# subject_labels = ['Subject 1', 'Subject 2']

# Group-specific colors (matched from your line plot)
group_colors = {
    'FAP_with_FAP_AAV': '#1f77b4',        # Blue
    'PBS_with_FAP_AAV': '#ff7f0e',        # Orange
    'PBS_with_Control_AAV': '#2ca02c'     # Green
}

# Create plot
fig, ax = plt.subplots(figsize=(6, 5))

# Plot each subject across all groups using group-specific colors
for i, subject_data in enumerate(data_matrix_T):
    ax.scatter(
        [1, 2, 3],                  # x = group positions
        subject_data,              # y = subject values
        s=80,
        color=[group_colors[label] for label in group_labels],
        # label=subject_labels[i]
    )

# Formatting
ax.set_xlim(0.5, 3.5)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(group_labels, rotation=15)
ax.set_ylabel('Signal Change (%)')
ax.set_title("Post Injection")
ax.grid(True)
ax.legend(loc='best')

plt.tight_layout()
plt.show()
