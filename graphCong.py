import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("AsyncCong10_33_.csv")

# Separate Load and Cong data
Cong = df[' Congested']
Load = df['Load']

# Transform the Load values
Load_transformed = Load / 2

# Calculate the midpoint of the dataset
midpoint = len(Load_transformed) // 2

plt.figure(figsize=(6, 6), dpi=100)
# Plot the first half of the data with larger points and connected lines
plt.plot(Load_transformed[:midpoint+1], Cong[:midpoint+1], marker='.', markersize=8, linestyle='-', linewidth=1, color='blue', label='First Half')

# Plot the second half of the data with larger points and connected lines
plt.plot(Load_transformed[midpoint:], Cong[midpoint:], marker='.', markersize=8, linestyle='-', linewidth=1, color='red', label='Second Half')


# plt.xticks(Load_transformed)
# Customize the plot
plt.xlabel('Number of vehicles per node')
plt.ylabel('Number of congested')
plt.title('Asynchronous Congestions')
plt.grid(True)
plt.legend()

# Save and show the plot
plt.savefig('AsyncCong10_33_.png')
plt.show()
