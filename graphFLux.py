import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("SyncFlux_10_91.csv")

# Separate Load and Flux data
Flux = df[' Flux']
Load = df['Load']

# Transform the Load values
Load_transformed = Load / 5

# Calculate the midpoint of the dataset
midpoint = len(Load_transformed) // 2

plt.figure(figsize=(16, 6), dpi=100)
# Plot the first half of the data with larger points and connected lines
plt.plot(Load_transformed[:midpoint+1], Flux[:midpoint+1], marker='.', markersize=8, linestyle='-', linewidth=1, color='blue', label='First Half')

# Plot the second half of the data with larger points and connected lines
plt.plot(Load_transformed[midpoint:], Flux[midpoint:], marker='.', markersize=8, linestyle='-', linewidth=1, color='red', label='Second Half')


# plt.xticks(Load_transformed)
# Customize the plot
plt.xlabel('Number of vehicles per node')
plt.ylabel('Flux')
plt.title('Synchronous Flux')
plt.grid(True)
plt.legend()

# Save and show the plot
plt.savefig('SyncFlux_10_91.png')
plt.show()
