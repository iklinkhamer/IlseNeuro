#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 22:00:35 2025

@author: Ilse Klinkhamer
"""

import os
import re
import shutil
import OpenEphys
from tkinter import Tk, filedialog
from datetime import datetime
from pathlib import Path

def select_folder(prompt="Select a folder"):
    root = Tk()
    root.withdraw()  # Hide the main window
    folder_path = filedialog.askdirectory(title=prompt)
    root.destroy()
    return folder_path

def find_full_source_name(folder, source, channel):
    pattern = re.compile(rf"^{re.escape(source)}(.*?)_CH{channel}\.continuous$")
    for filename in os.listdir(folder):
        match = pattern.match(filename)
        if match:
            return source + match.group(1)
    return None  # or raise an error if preferred

def find_session_and_mouse_name(path):
    path = Path(path).resolve()

    for parent in [path] + list(path.parents):
        name = parent.name

        # Remove trailing session number if present
        base = name.rsplit('_', 1)[0] if name.rsplit('_', 1)[-1].isdigit() else name

        # Extract digits only
        digits = ''.join(c for c in base if c.isdigit())

        # Try parsing with different datetime formats
        for fmt in ("%Y%m%d%H%M%S", "%Y%m%d%H%M", "%Y%m%d"):
            try:
                datetime.strptime(digits, fmt)
                mouse_name = parent.parent.name  # Get the parent folder name
                return name, mouse_name
            except ValueError:
                continue

    return None, None  # No session-like folder found


def find_best_record_node_folder(session_folder_open_ephys_data):
    record_node_folders = []

    # Look for subfolders matching "Record Node 10#"
    for name in os.listdir(session_folder_open_ephys_data):
        match = re.match(r"Record Node (\d+)", name)
        if match:
            node_number = int(match.group(1))
            full_path = os.path.join(session_folder_open_ephys_data, name)

            if os.path.isdir(full_path):
                # Check for files starting with number ≥ 100 and ending with _CH1.continuous
                for fname in os.listdir(full_path):
                    file_match = re.match(r"(\d+)_.*_CH1\.continuous$", fname)
                    if file_match and int(file_match.group(1)) >= 100:
                        record_node_folders.append((node_number, full_path))
                        break  # no need to check other files in this folder

    # If we found any matching folders, return the one with the highest number
    if record_node_folders:
        best_folder = max(record_node_folders, key=lambda x: x[0])[1]
        return best_folder
    else:
        return None



"""
# Example usage
if __name__ == "__main__":
    folder = select_folder()
    if folder:
        print(f"You selected: {folder}")
    else:
        print("No folder selected.")
"""

def convertOpenEphysDataToContinuous(   open_ephys_data_directory=None#os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
                                     ,  kilosort_output_directory=None
                                     ,  destination_directory=None
                                     , channels=list(range(1,33))
                                     ):
    if open_ephys_data_directory is None:
        open_ephys_data_directory=os.path.join(select_folder("Select the folder with raw OpenEphys data"))
    if kilosort_output_directory is None:
        kilosort_output_directory=os.path.join(select_folder("Select the phy compatible folder with kilosort input and output data"))
    if destination_directory is None:
        destination_directory=os.path.join(select_folder("Select the destination folder for c4"))
    
            
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
            
            
def main(open_ephys_data_directory=None#os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
                                 ,  kilosort_output_directory=None
                                 ,  destination_directory=None
                                 ):
    
    if open_ephys_data_directory is None:
        open_ephys_data_directory=os.path.join(select_folder("Select the folder with raw OpenEphys data"))
    if kilosort_output_directory is None:
        kilosort_output_directory=os.path.join(select_folder("Select the phy compatible folder with kilosort input and output data"))
    if destination_directory is None:
        destination_directory=os.path.join(select_folder("Select the destination folder for c4"))
    # Find all valid session subfolders
    subfolders_open_ephys = []
    for f in os.listdir(open_ephys_data_directory):
        full_path = os.path.join(open_ephys_data_directory, f)
        if os.path.isdir(full_path):
            session_name, mouse_name = find_session_and_mouse_name(full_path)
            if session_name is not None and mouse_name is not None:
                full_path = os.path.join(full_path)
                subfolders_open_ephys.append(full_path)
                subfolders_open_ephys.sort()
        
    subfolders_kilosort = []
    for f in os.listdir(open_ephys_data_directory):
        full_path = os.path.join(open_ephys_data_directory, f)
        if os.path.isdir(full_path):
            session_name, mouse_name = find_session_and_mouse_name(full_path)
            if session_name is not None and mouse_name is not None:
                full_path = os.path.join(full_path)
                subfolders_kilosort.append(full_path)
                subfolders_kilosort.sort()
        
    subfolders_destination = []
    for f in os.listdir(open_ephys_data_directory):
        full_path = os.path.join(open_ephys_data_directory, f)
        if os.path.isdir(full_path):
            session_name, mouse_name = find_session_and_mouse_name(full_path)
            if session_name is not None and mouse_name is not None:
                full_path = os.path.join(full_path)
                subfolders_destination.append(full_path)
                subfolders_destination.sort()
 
    if subfolders_destination and subfolders_kilosort and subfolders_open_ephys:
        for session_folder_open_ephys_data, session_folder_kilosort, session_folder_destination in zip(subfolders_open_ephys, subfolders_kilosort, subfolders_destination):
            print(f"\n--- Processing: {session_folder_open_ephys_data} ---")
    
            
            best_node_folder = find_best_record_node_folder(session_folder_open_ephys_data)
            if best_node_folder:
                best_node_folder = os.path.join(best_node_folder)
                print(f"Best Record Node folder: {best_node_folder}")
            else:
                #best_node_folder = session_folder_open_ephys_data
                print("No suitable Record Node folder found.")            
          
            convertOpenEphysDataToContinuous(best_node_folder, session_folder_kilosort, session_folder_destination)
    else:
        convertOpenEphysDataToContinuous(open_ephys_data_directory, kilosort_output_directory, destination_directory)

        
if __name__ == "__main__":
    main()

