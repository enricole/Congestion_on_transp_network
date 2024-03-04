import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data from CSV
df = pd.read_csv("SyncRho_8.csv")

# Extracting data from the DataFrame
Probability = df[' Probability']
State = df['State']

# Calculate the mean value of the distribution
mean_state = np.average(State, weights=Probability)


# Plotting the data
plt.figure(figsize=(12, 6), dpi=100)
plt.plot(State, Probability, marker='.', linestyle='-', label='Data')

# Setting the y-axis scale starting from 0.4
plt.ylim(0, 0.35)

# Annotate the plot with the mean value
plt.text(8, 0.3, f'Mean State: {mean_state:.2f}', ha='center', va='bottom')

# Labeling and titling the plot
plt.xlabel('State')
plt.ylabel('Probability')
plt.title('Synchronous Density Distribution')
plt.grid(True)

# Saving and showing the plot
plt.savefig('SyncRho_8.png')
plt.show()
