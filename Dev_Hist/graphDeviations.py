import pandas as pd
import matplotlib.pyplot as plt

# Load the data from CSV
df = pd.read_csv("SyncDev_20_81.csv")

# Extracting data from the DataFrame
Deviation = df[' Deviation']
Load = df['Load']

# Transform the Load values
Load_transformed = Load / 2

# Calculate the midpoint of the dataset
midpoint = len(Load_transformed) // 2

# Plotting the data
plt.figure(figsize=(16, 6), dpi=100)
# Plot the first half of the data with larger points and connected lines
plt.plot(Load_transformed[:midpoint+1], Deviation[:midpoint+1], marker='.', markersize=8, linestyle='-', linewidth=1, color='blue', label='First Half')

# Plot the second half of the data with larger points and connected lines
plt.plot(Load_transformed[midpoint:], Deviation[midpoint:], marker='.', markersize=8, linestyle='-', linewidth=1, color='red', label='Second Half')

# Generate custom labels for x-axis
# custom_labels = [15 + 20*i for i in range(len(Load))]
# custom_labels = [30*i for i in range(len(Load))]
# Set custom ticks and labels for x-axis
# plt.xticks(Load[::10], custom_labels[::10])  # Show every 10th label

# Labeling and titling the plot
plt.xlabel('Load')
plt.ylabel('Deviations')
plt.title('Deviations vs Load')
plt.grid(True)
plt.legend()


# Saving and showing the plot
plt.savefig('SyncDev_20_81.png')
plt.show()
