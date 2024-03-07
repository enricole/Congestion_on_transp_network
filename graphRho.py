import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data from CSV
df = pd.read_csv("SyncRho_15.csv")

# Extracting data from the DataFrame
Probability = df[' Probability']
State = df['State']

# Calculate the weighted standard deviation
weighted_mean = np.average(State, weights=Probability)
weighted_variance = np.average((State - weighted_mean) ** 2, weights=Probability)
weighted_std_dev = np.sqrt(weighted_variance)

# Plotting the data
plt.figure(figsize=(12, 6), dpi=100)
plt.plot(State, Probability, marker='.', linestyle='-', label='Data')

# Setting the y-axis scale starting from 0.4
plt.ylim(0, 0.2)


# Annotate the plot with the mean value and standard deviation
plt.text(15, 0.15, f'Mean State: {weighted_mean:.2f}', ha='center', va='top')
plt.text(15.15, 0.14, f'Standard Deviation: {weighted_std_dev:.2f}', ha='center', va='top')


# Labeling and titling the plot
plt.xlabel('State')
plt.ylabel('Probability')
plt.title('Synchronous Density Distribution')
plt.grid(True)

# Saving and showing the plot
plt.savefig('SyncRho_15.png')
plt.show()
