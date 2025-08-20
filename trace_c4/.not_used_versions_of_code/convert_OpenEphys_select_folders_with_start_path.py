#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:02:25 2025

@author: Ilse Klinkhamer
"""
import os
import re
import shutil
import OpenEphys
from tkinter import Tk, filedialog
from pathlib import Path
import tkinter as tk

LAST_FOLDER_FILE = Path.home() / ".last_folder_convert_open_ephys.txt"  # Hidden file in home dir

def select_folder(prompt="Select a folder", start_path=None):
    # Load last folder if available, unless a start_path is provided
    if start_path:
        initial_dir = start_path
    else:
        initial_dir = LAST_FOLDER_FILE.read_text() if LAST_FOLDER_FILE.exists() else str(Path.home())

    # Open folder selection dialog
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(title=prompt, initialdir=initial_dir)

    # Save selected folder for next time
    if folder_path:
        LAST_FOLDER_FILE.write_text(folder_path)

    return folder_path if folder_path else None


def find_full_source_name(folder, source, channel):
    pattern = re.compile(rf"^{re.escape(source)}(.*?)_CH{channel}\.continuous$")
    for filename in os.listdir(folder):
        match = pattern.match(filename)
        if match:
            return source + match.group(1)
    return None  # or raise an error if preferred


"""
# Example usage
if __name__ == "__main__":
    folder = select_folder()
    if folder:
        print(f"You selected: {folder}")
    else:
        print("No folder selected.")
"""

def convertOpenEphysDataToContinuous(   channels=list(range(1,33))
                                     ,  open_ephys_data_directory=None#os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
                                     ,  kilosort_output_directory=None
                                     ,  destination_directory=None
                                     ):
    
        open_ephys_data_directory = select_folder(
            "Select the folder with raw OpenEphys data",
            start_path="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/pythonpipeline/data/ephys/continuous_data/"
        )
        
        print(f"{open_ephys_data_directory}")
        
        kilosort_output_directory = select_folder(
            "Select the phy compatible folder with kilosort input and output data",
            start_path="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/pythonpipeline/data/ephys/kilosort_output/"
        )
        
        destination_directory = select_folder(
            "Select the destination folder for c4",
            start_path="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/Ilse/ephys/"
        )

    
            
        print(f"Processing folder: {open_ephys_data_directory}")
        
        # Perform actions on each folder here
        source_folder = open_ephys_data_directory
        if "c4" not in destination_directory:
            destination_directory = os.path.join(destination_directory, "c4")
            
        destination_folder_continuous = os.path.join(destination_directory, "continuous", "Data_AP_LFP")
        
        # Define file paths
        source_file = os.path.join(source_folder, "openephys.dat")
        destination_file = os.path.join(destination_folder_continuous, "continuous.dat")

        if os.path.exists(destination_file):
           print("Continuous file already exists, skipping.")
           return

        if not os.path.exists(source_folder):
            print("Data folder not found, skipping.")
            return
           
        for filename in os.listdir(source_folder):
            if filename.endswith(".continuous"):
                match = re.match(r"(\d{3})", filename)
                if match:
                    match = str(match.group(1))  # Convert to an integer
                    break  # Stop after finding the first match
        full_source = find_full_source_name(source_folder, match, channels[0])
        
       
        OpenEphys.pack_2(folderpath=source_folder, filename="openephys.dat", source=full_source, channels = channels)
        
        # Create the destination folder if it doesn't exist
        os.makedirs(destination_folder_continuous, exist_ok=True)
        
        # Move and rename the file
        if os.path.exists(source_file):
            shutil.move(source_file, destination_file)
            print(f"File moved and renamed to: {destination_file}")
        else:
            print(f"Source file does not exist: {source_file}")
                
        #also copy and paste the oebin file because that's easier
        # Define file paths
        #foldername2 = "Quimper_20230801153833 (copy)"
        #source_folder2 = os.path.join(mouse_folder, foldername2, "Extraction2Bin")
        source_folder_oebin = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ComplexSpikeToolkit/Ilse/Ilse PhD/oebin/"
        source_file_oebin = os.path.join(source_folder_oebin, "structureIK.oebin")
        destination_file_oebin = os.path.join(destination_directory, "structureIK.oebin")
        
        # Copy and rename the file
        if os.path.exists(source_file_oebin):
            shutil.copy(source_file_oebin, destination_file_oebin)  # Copy file
            print(f"File copied to: {destination_file_oebin}")
        else:
            print(f"Source file does not exist: {source_file_oebin}")
            
        # copy the content of the Extraction2Bin folder to the c4 folder
        # Define paths
        source_folder = os.path.join(kilosort_output_directory)
        destination_folder_c4 = destination_directory
        if source_folder == destination_folder_c4:
            return
        
        # Define the files to exclude
        excluded_files = {'Data4KS2.bin', 'temp_wh.dat', 'rez.mat', 'pc_features.npy', 'template_features.npy', 'kilosort_input.bin'}  # Add all filenames to exclude
        #included_files = {'amplitudes.npy', 'cluster_group.tsv', 'params.py', 'spike_clusters.npy', 'spike_times.npy'}
        
        def ignore_files(dir, files):
            """Custom ignore function to exclude specific files."""
            return {file for file in files if file in excluded_files}
        
        # Copy the folder, excluding the specified files
        if os.path.exists(source_folder):
            shutil.copytree(
                source_folder,
                destination_folder_c4,
                dirs_exist_ok=True,  # Allow overwriting if destination exists (Python 3.8+)
                ignore=ignore_files
            )
            print(f"Contents of {source_folder} copied to {destination_folder_c4}")
        else:
            print(f"Source folder does not exist: {source_folder}")
            
            
def main():
    convertOpenEphysDataToContinuous()
        
if __name__ == "__main__":
    main()