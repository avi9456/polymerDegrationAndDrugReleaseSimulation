# from numba import jit, cuda
import numpy as np
import csv
import os
import pandas as pd
import json 
# @jit(target_backend='cuda') 
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
        print(f"Error saving json: {e}")

def exportDataPandas(releaseData,ndrug,trapped,polymerSize,polymerCon,drugCon, filename):
    # Combine all the input data into a single DataFrame
    releaseData_df = pd.DataFrame({ 
        "MCstep": list(releaseData[1:,0].astype(float)),
        "releasePercent": list(releaseData[1:,1].astype(float)),
        "cumDrugRelease": list(releaseData[1:,2].astype(float)),
        "drugRelease": list(releaseData[1:,3].astype(float)),
        "remainDrug": list(releaseData[1:,4].astype(float)),
        })
    data = ({
        "dataInfo": {
            "polymerCon": polymerCon,
            'polymerSize': polymerSize,
            "ndrug": ndrug,
            "trapped": trapped,
            "drugCon": drugCon,
        },
        "releaseData": releaseData_df.to_dict(orient='list'),
        # {
        #     "MCstep": list(releaseData[1:,0].astype(float)),
        #     "releasePercent": list(releaseData[1:,1].astype(float)),
        #     "cumDrugRelease": list(releaseData[1:,2].astype(float)),
        #     "drugRelease": list(releaseData[1:,3].astype(float)),
        #     "remainDrug": list(releaseData[1:,4].astype(float)),
        # },
    })

    # if not isinstance(data, (np.ndarray, pd.DataFrame)):
    #     raise TypeError("Input must be a NumPy array or a pandas DataFrame.")
    if not isinstance(filename, str) or not filename.endswith(".json"):
        raise ValueError("Filename must be a string ending with '.json'.")
    try:
        # Ensure target directory exists
        directory = os.path.dirname(filename)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)

        # Convert to DataFrame if it's a NumPy array
        # if isinstance(data, np.ndarray):
        #     data = pd.DataFrame(data)
        # data = data.to_dict(orient='index')  # Convert DataFrame to dictionary for JSON serialization
        # Save using pandas to_json (handles both numeric and non-numeric data)
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4,)
        print(f"Data saved to {filename}")

    except Exception as e:
        print(f"Error saving JSON: {e}")


