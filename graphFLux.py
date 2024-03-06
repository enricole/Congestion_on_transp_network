import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("SyncFlux_10_37.csv")
df2 = pd.read_csv("SyncDev_10_37.csv")

print(df2.columns)
# Separate Load and Flux data
Flux = df[' Flux']
Deviation = df2[' Deviation']
Load = df['Load']

# Transform the Load values
Load_transformed = Load / 2


# Calculate the midpoint of the dataset
midpoint = len(Load_transformed) // 2

# Create a figure and axis object with twinx for secondary y-axis
fig, ax1 = plt.subplots(figsize=(6, 11), dpi=100)
ax2 = ax1.twinx()

# Plot the flux data on the left y-axis
# Plot the first half of the data with larger points and connected lines
ax1.plot(Load_transformed[:midpoint+1], Flux[:midpoint+1], marker='.', markersize=8, linestyle='-', linewidth=1, color='blue', label='First Half')

# Plot the second half of the data with larger points and connected lines
ax1.plot(Load_transformed[midpoint:], Flux[midpoint:], marker='.', markersize=8, linestyle='-', linewidth=1, color='red', label='Second Half')


# Plot the first half of the data with larger points and connected lines
ax2.plot(Load_transformed[:midpoint+1], Deviation[:midpoint+1], marker='.', markersize=8, linestyle='--', linewidth=1, color='blue', label='First Half')

# Plot the second half of the data with larger points and connected lines
ax2.plot(Load_transformed[midpoint:], Deviation[midpoint:], marker='.', markersize=8, linestyle='--', linewidth=1, color='red', label='Second Half')
# # Plot the standard deviation on the right y-axis
# ax2.plot(Load_transformed, Deviation, color='red', linestyle='--', label='Flux Standard Deviation')

# Customize the plot
ax1.set_xlabel('Number of vehicles per node')
ax1.set_ylabel('Flux')
ax2.set_ylabel('Flux Standard Deviation')
ax1.set_title('Synchronous Flux')
ax1.grid(True)

# Show legend for both axes
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left')

# Save and show the plot
plt.savefig('SyncFlux_10_37.png')
plt.show()
