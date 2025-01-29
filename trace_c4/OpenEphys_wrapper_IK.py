#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:14:00 2024

@author: Ilse Klinkhamer
"""

import OpenEphys
import os
import shutil

channels = list(range(1, 33))
mouse_name = "Yosemite"

# Define paths

directory = f"/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/{mouse_name}/"

# Get all folders with mouse_name in their name and without "copy"
mouse_folders = [
    folder for folder in os.listdir(directory)
    if os.path.isdir(os.path.join(directory, folder)) 
    and mouse_name in folder 
    and "copy" not in folder
]

mouse_folders.sort()


# Print the list of matching folders
print(mouse_folders)

for foldername in mouse_folders:
    
    print(f"Processing folder: {foldername}")
    
    # Perform actions on each folder here
    source_folder = os.path.join(directory, foldername, "Data")
    destination_folder = os.path.join(directory, foldername, "c4", "continuous", "Data_AP_LFP")
    
    # Create the destination folder if it doesn't exist
    os.makedirs(destination_folder, exist_ok=True)
    
    # Define file paths
    source_file = os.path.join(source_folder, "openephys.dat")
    destination_file = os.path.join(destination_folder, "continuous.dat")
       
    if os.path.exists(destination_file):
       print("Continuous file already exists, skipping.")
       continue
   
    OpenEphys.pack(source_folder, channels = channels)   
    
    
    # Move and rename the file
    if os.path.exists(source_file):
        shutil.move(source_file, destination_file)
        print(f"File moved and renamed to: {destination_file}")
    else:
        print(f"Source file does not exist: {source_file}")
        


    #also copy and paste the oebin file because that's easier
    # Define file paths
    #foldername2 = "Quimper_20230801153833 (copy)"
    #source_folder2 = os.path.join(directory, foldername2, "Extraction2Bin")
    source_folder_oebin = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ComplexSpikeToolkit/Ilse/Ilse PhD/oebin/"
    source_file_oebin = os.path.join(source_folder_oebin, "structureIK.oebin")
    destination_file_oebin = os.path.join(directory, foldername, "c4", "structureIK.oebin")
    
    # Copy and rename the file
    if os.path.exists(source_file_oebin):
        shutil.copy(source_file_oebin, destination_file_oebin)  # Copy file
        print(f"File copied to: {destination_file_oebin}")
    else:
        print(f"Source file does not exist: {source_file_oebin}")
        
    # copy the content of the Extraction2Bin folder to the c4 folder
    # Define paths
    source_folder = os.path.join(directory, foldername, "Extraction2Bin")
    destination_folder = os.path.join(directory, foldername, "c4")
    
    # Define the files to exclude
    excluded_files = {'Data4KS2.bin', 'temp_wh.dat'}  # Add all filenames to exclude
    
    def ignore_files(dir, files):
        """Custom ignore function to exclude specific files."""
        return {file for file in files if file in excluded_files}
    
    # Copy the folder, excluding the specified files
    if os.path.exists(source_folder):
        shutil.copytree(
            source_folder,
            destination_folder,
            dirs_exist_ok=True,  # Allow overwriting if destination exists (Python 3.8+)
            ignore=ignore_files
        )
        print(f"Contents of {source_folder} copied to {destination_folder}, excluding {excluded_files}")
    else:
        print(f"Source folder does not exist: {source_folder}")