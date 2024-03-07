import csv
import numpy as np
import matplotlib.pyplot as plt

def read_csv(filename):
    data = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            data.append(float(row[0]))  # Assuming your data is in the first column of the CSV file
    return data

def count_occurrences(data, intervals):
    counts = [0] * (len(intervals) - 1)
    for value in data:
        for i in range(len(intervals) - 1):
            if intervals[i] <= value < intervals[i + 1]:
                counts[i] += 1
                break
    return counts

def plot_histogram(intervals, counts):
    interval_midpoints = [((intervals[i] + intervals[i+1]) / 2) for i in range(len(intervals)-1)]
    plt.bar(interval_midpoints, counts, width=np.diff(intervals), align='center', edgecolor='black')
    plt.title('Occurence vs Deviations')
    plt.xlabel('Deviations')
    plt.ylabel('Occurence')
    plt.xticks(interval_midpoints, intervals[:-1], rotation=45)  # Set x-ticks to interval midpoints and labels to interval start points
    plt.grid(True)
    plt.savefig('AsyncHist_1.png')  # Save the plot as 'histogram.png' in the current directory
    plt.show()

def main():
    # Replace 'data.csv' with the path to your CSV file
    filename = 'AsyncHist_1.csv'
    
    # Define intervals
    intervals = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]  # Adjust intervals as needed
    
    # Read data from CSV file
    data = read_csv(filename)
    
    # Count occurrences in intervals
    counts = count_occurrences(data, intervals)
    
    # Plot histogram
    plot_histogram(intervals, counts)

if __name__ == "__main__":
    main()
