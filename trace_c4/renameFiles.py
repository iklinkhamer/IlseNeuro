#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 15:50:54 2024

@author: no1
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 17:41:19 2024

@author: no1
"""
import os
import shutil

channels = list(range(1, 33))
mouse_name = "Quimper"

# Define paths

directory = f"/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/{mouse_name}/"

# Get all folders with mouse_name in their name and without "copy"
mouse_folders = [
    folder for folder in os.listdir(directory)
    if os.path.isdir(os.path.join(directory, folder)) 
    and mouse_name in folder 
    and "copy" not in folder
]

# Loop through each processed folder
for foldername in mouse_folders:
    print(f"Processing folder: {foldername}")
    
    # Define the path for the `cell_type_classification` folder
    source_folder = os.path.join(directory, foldername, "c4", "res_10perc_old_figs_log_cwin2000")
    
    # Define the new folder where it will be moved
    destination_folder = os.path.join(directory, foldername, "c4", "res_5perc_old_figs_log_cwin2000")
    os.makedirs(destination_folder, exist_ok=True)
    # Check if the cell_type_classification folder exists
    if os.path.exists(source_folder):    
        # Move the folder
        shutil.move(source_folder, destination_folder)
        print(f"Moved {destination_folder} folder to: {destination_folder}")
    else:
        print(f"{source_folder} folder does not exist in {foldername}")
    
 