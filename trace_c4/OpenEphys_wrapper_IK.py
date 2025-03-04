#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:14:00 2024

@author: Ilse Klinkhamer
"""

import OpenEphys
import os
import shutil  
import sys
import re
import stitch_CH_continuous_files_switch_sessions
from exceptiongroup import catch
from get_dropbox_path import get_dropbox_path


def convertOpenEphysDataToContinuous(   mouse_name
                                     ,  switch_sessions=False
                                     ,  channels=list(range(1,33))
                                     ,  directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
                                     ):

    dp_base = os.path.join(directory, mouse_name)
    if "ReserveMouse" in mouse_name:
        dp_base = dp_base.replace("MainFolder", "ReserveFolder")

    # Get all folders with mouse_name in their name and without "copy"
    mouse_folders = [
        folder for folder in os.listdir(dp_base)
        if os.path.isdir(os.path.join(dp_base, folder))
           and mouse_name in folder
           and "copy" not in folder
           and "Copy" not in folder
    ]
    mouse_folders.sort()

    if switch_sessions:

        switch_folder = os.path.join(dp_base, "SwitchSessionStitching")
        if os.path.exists(switch_folder):
            continuous_switch_file = os.path.join(switch_folder, "continuous.dat")

            if os.path.exists(continuous_switch_file):
                mouse_folders.append(switch_folder)
            else:
                switch_folder_data_file = os.path.join(switch_folder, "Data", "100_CH1.continuous")
                if not os.path.exists(switch_folder_data_file):
                    try:
                        stitch_CH_continuous_files_switch_sessions.main(mouse_name)
                        mouse_folders.append(switch_folder)
                    except:
                        print("Stitching switch-sessions failed. Skipping analysis for switch sessions")
                else:
                    mouse_folders.append(switch_folder)



                # Print the list of matching folders
    print(mouse_folders)
    
    for foldername in mouse_folders:
        
        print(f"Processing folder: {foldername}")
        
        # Perform actions on each folder here
        source_folder = os.path.join(dp_base, foldername, "Data")
        destination_folder = os.path.join(dp_base, foldername, "c4", "continuous", "Data_AP_LFP")
        
        # Create the destination folder if it doesn't exist
        os.makedirs(destination_folder, exist_ok=True)
        
        # Define file paths
        source_file = os.path.join(source_folder, "openephys.dat")
        destination_file = os.path.join(destination_folder, "continuous.dat")

        if os.path.exists(destination_file):
           print("Continuous file already exists, skipping.")
           continue

        if not os.path.exists(source_folder):
            print("Data folder not found, skipping.")
            continue
           
        for filename in os.listdir(source_folder):
            if filename.endswith(".continuous"):
                match = re.match(r"(\d{3})", filename)
                if match:
                    match = str(match.group(1))  # Convert to an integer
                    break  # Stop after finding the first match
       
        OpenEphys.pack_2(folderpath=source_folder, filename="openephys.dat", source=match, channels = channels)
        
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
        destination_file_oebin = os.path.join(dp_base, foldername, "c4", "structureIK.oebin")
        
        # Copy and rename the file
        if os.path.exists(source_file_oebin):
            shutil.copy(source_file_oebin, destination_file_oebin)  # Copy file
            print(f"File copied to: {destination_file_oebin}")
        else:
            print(f"Source file does not exist: {source_file_oebin}")
            
        # copy the content of the Extraction2Bin folder to the c4 folder
        # Define paths
        source_folder = os.path.join(dp_base, foldername, "Extraction2Bin")
        destination_folder = os.path.join(dp_base, foldername, "c4")
        
        # Define the files to exclude
        excluded_files = {'Data4KS2.bin', 'temp_wh.dat', 'rez.mat', 'pc_features.npy', 'template_features.npy'}  # Add all filenames to exclude
        
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
            
def main(   mouse_name=None
         ,  switch_sessions=True):
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
        else:
            print("Error: No mouse name provided")
            sys.exit(1)
    convertOpenEphysDataToContinuous(   mouse_name
                                     ,  switch_sessions=switch_sessions
                                     )    
        
if __name__ == "__main__":
    main()