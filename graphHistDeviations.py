import pandas as pd
import matplotlib.pyplot as plt

# Load the data from CSV
df = pd.read_csv("SyncHist_10.csv")

# Extracting data from the DataFrame
Occurence = df[' Occurence']
Deviation = df['Deviation']

# Plotting the data
plt.figure(figsize=(16, 6), dpi=100)
# Plot the first half of the data with larger points and connected lines
plt.plot(Deviation, Occurence, marker='.', markersize=8, linestyle='-', linewidth=1, color='blue', label='First Half')

# Generate custom labels for x-axis
custom_labels = [f'{Deviation[i]:.2f} - {Deviation[i+1]:.2f}' for i in range(len(Deviation)-1)]  # Generating labels as "0.00 - 0.50", "0.50 - 1.00", etc.
# Set custom ticks and labels for x-axis
plt.xticks(Deviation[:-1], custom_labels)  # Exclude the last point for the ticks

# Labeling and titling the plot
plt.xlabel('Deviation')
plt.ylabel('Occurrence')
plt.title('Occurences vs Deviations')
plt.grid(True)

# Saving and showing the plot
plt.savefig('SyncHist_10.png')
plt.show()
