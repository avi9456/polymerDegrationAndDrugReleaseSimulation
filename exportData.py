import numpy as np
import csv
import os

def exportData(data,filename):
    if not isinstance(data, np.ndarray):
        raise TypeError("Input must be a NumPy array.")
    if not isinstance(filename, str) or not filename.endswith(".csv"):
        raise ValueError("Filename must be a string ending with '.csv'.")
    try:
        # Ensure target directory exists
        directory = os.path.dirname(filename)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)

        # Save using numpy.savetxt (handles numeric arrays)
        np.savetxt(filename, data, delimiter=",", fmt="%s")
        print(f"Data saved to {filename}")

    except Exception as e:
        print(f"Error saving CSV: {e}")

