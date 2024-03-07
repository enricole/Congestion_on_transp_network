import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("AsyncFlux_20_81.csv")

# Separate Load and Flux data
Flux = df[' Flux']
Load = df['Load']

# Modifica i valori di Load
Load_modified = Load * 0.5

plt.figure(figsize=(10, 6), dpi=100)
# Plot the flux data on the left y-axis
# Plot the first half of the data with larger points and connected lines
plt.plot(Load_modified, Flux, marker='.', markersize=8, linestyle='-', linewidth=1, color='blue')

# Customize the plot
plt.xlabel('Number of vehicles per node')
plt.ylabel('Flux')
plt.title('Asynchronous Flux')
plt.grid(True)
plt.legend()

# Save and show the plot
plt.savefig('AsyncFlux_20_81.png')
plt.show()
