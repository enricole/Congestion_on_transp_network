import pandas as pd
import matplotlib.pyplot as plt

# Load the data from CSV
df = pd.read_csv("AsyncDev_20_81.csv")

# Extracting data from the DataFrame
deviation = df[' Deviation']
Load = df['Load']

# Plotting the data
plt.figure(figsize=(11, 6), dpi=100)
plt.plot(Load, deviation, marker='.', linestyle='-', label='Data')

# Generate custom labels for x-axis
custom_labels = [0.5 * i for i in range(len(Load))]

# Set custom ticks and labels for x-axis
plt.xticks(Load, custom_labels, rotation=90)  # Show all labels, rotated for better readability

# Labeling and titling the plot
plt.xlabel('Number of vehicles per node')
plt.ylabel('Deviations')
plt.title('Deviations vs Load')
plt.grid(True)

# Saving and showing the plot
plt.savefig('AsyncDev_20_81.png')
plt.show()
