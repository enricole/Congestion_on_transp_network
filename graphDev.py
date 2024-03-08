import pandas as pd
import matplotlib.pyplot as plt

# Load the data from CSV
df = pd.read_csv("AsyncDeviations_8.csv")

# Extracting data from the DataFrame
deviation = df[' deviation']
iterations = df['iterations']

# Plotting the data
plt.figure(figsize=(16, 6), dpi=100)
plt.plot(iterations, deviation, marker='.', linestyle='-', label='Data')

# Generate custom labels for x-axis
# custom_labels = [15 + 20*i for i in range(len(iterations))]
custom_labels = [30*i for i in range(len(iterations))]
# Set custom ticks and labels for x-axis
plt.xticks(iterations[::10], custom_labels[::10])  # Show every 10th label

# Labeling and titling the plot
plt.xlabel('Iterations')
plt.ylabel('Deviations')
plt.title('Deviations vs Iterations')
plt.grid(True)

# Saving and showing the plot
plt.savefig('AsyncDeviations_8.png')
plt.show()
