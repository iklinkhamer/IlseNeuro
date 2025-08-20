# -*- coding: utf-8 -*-
"""
Created on Fri May 23 12:05:41 2025

@author: Ilse Klinkhamer
"""

# Now you can inspect, for example:
# print(mice['Rotterdam'])
from scipy.io import loadmat
import numpy as np
import os

def recurse_struct(data, indent=0):
    spacing = '  ' * indent
    if isinstance(data, dict):
        for key, value in data.items():
            print(f"{spacing}{key}:", file=out)
            recurse_struct(value, indent + 1)
    elif hasattr(data, '__dict__'):
        for key, value in vars(data).items():
            print(f"{spacing}{key}:", file=out)
            recurse_struct(value, indent + 1)
    elif isinstance(data, np.ndarray):
        print(f"{spacing}[array shape {data.shape}]", file=out)
    else:
        print(f"{spacing}{data}", file=out)

# Load the .mat file
directory = r"C:\Users\Ilse Klinkhamer\Downloads\metrics"
file = "n-1_results_Narrow.mat"
file_path = os.path.join(directory, file)
mat_data = loadmat(file_path, struct_as_record=False, squeeze_me=True)

# Write everything to a text file
with open(os.path.join(directory,'output.txt'), 'w') as out:
    recurse_struct(mat_data)


a = 1
