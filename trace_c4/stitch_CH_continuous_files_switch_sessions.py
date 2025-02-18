#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:18:59 2025

@author: Ilse Klinkhamer
"""

import os
import re

def read_continuous_file(filename, skip_header=False, data_folder="Data"):
    
    """Reads a .continuous file, optionally skipping the first 1024 bytes."""
    with open(filename, 'rb') as f:
        if skip_header:
            f.seek(1024)  # Skip the 1024-byte header
        return f.read()
    
    
def get_file_names(mouse_name
                   , folders
                   , mouse_folder 
                   , channel_number=0
                   , data_folder_name="Data"
                   , channels=list(range(1,33))
                   , source='100'
                   , chprefix='CH'
                   ):   

    files = []

    """Gets file paths for the available folders (2 or 3)."""
    for f in range(len(folders)):
        session_folder=f"{mouse_name}_{folders[f]}"
        source_folder = os.path.join(mouse_folder, session_folder, data_folder_name)
        for filename in os.listdir(source_folder):
            if filename.endswith(".continuous"):
                source = re.match(r"(\d{3})", filename)
                if source:
                    source = str(source.group(1))  # Convert to an integer
                    break  # Stop after finding the first source
        filename = f"{source}_{chprefix}{channel_number}.continuous"
        files.append(os.path.join(mouse_folder, session_folder, data_folder_name, filename))
                        
                        
        #filelist = [source + '_'+chprefix + x + '.continuous' for x in map(str,channels)]
        #files = [f"{mouse_folder}/{mouse_name}{folder}/100_CH{idx + 1}.continuous" for idx, folder in enumerate(folders)]
 
            
    return files

    
    
def read_folders_from_file(mouse_name
                           , mouse_folder
                           , filename="sessionIds"
                           , switch_folder_name="SwitchSessionStitching"):
    """Reads the first 2 or 3 lines of a file as folder names."""
    
    file_path = f"{mouse_folder}/{switch_folder_name}/{filename}"

    with open(file_path, 'r') as f:
        folders = [line.strip() for line in f.readlines()[:3]]  # Read up to 3 lines, strip whitespace

    return folders  # Could return 2 or 3 folders depending on file content


def main(mouse_name="Ana4"
         , directory="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/" 
         , switch_folder_name="SwitchSessionStitching"
         , data_folder_name="Data"
         , output_folder_name="c4/continuous/Data_AP_LFP"
         , channels=list(range(1,33))
         ):    


    mouse_folder = os.path.join(directory, mouse_name)
    if "ReserveMouse" in mouse_name:
        mouse_folder = mouse_folder.replace("MainFolder", "ReserveFolder")
    destination_folder = os.path.join(mouse_folder, switch_folder_name, data_folder_name)
    os.makedirs(destination_folder, exist_ok=True)
        
    
    folders = read_folders_from_file(  mouse_name
                                     , mouse_folder
                                     , switch_folder_name=switch_folder_name
                                     )
    
    #print("Folders to check:", folders)
    for c in range(len(channels)):
        files = get_file_names(  mouse_name
                               , folders
                               , mouse_folder
                               , channel_number=channels[c]
                               , data_folder_name=data_folder_name
                               )
        
        #print("Files to merge:", files)
    
        # Read and concatenate files
        combined_data = b""  # Initialize empty binary string
        for idx, file in enumerate(files):
            skip_header = (idx > 0)  # Skip header for all except the first file
            combined_data += read_continuous_file(  file
                                                  , skip_header=skip_header
                                                  )
    
    
        # Save the merged file
        output_filename = os.path.join(destination_folder, f"100_CH{channels[c]}.continuous")
     
    
        with open(output_filename, "wb") as out_file:
            out_file.write(combined_data)
    
        print(f"Merged 100_CH{channels[c]}.continuous file created successfully!")


if __name__ == "__main__":
    main()




"""
   if file3:
       output_filename += "_CH3"
   output_filename += ".continuous"
"""