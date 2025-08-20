#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:36:03 2025

@author: Ilse Klinkhamer
"""

import numpy as np
from scipy.io import loadmat, savemat
from get_dropbox_path import get_dropbox_path
import os

def load_and_save_as_mat(   spike_file
                         ,  event_file
                         ,  destination_folder
                         ,  sampling_rate=30000):
    spike_times = np.load(spike_file) / sampling_rate  # Convert frames to seconds
    event_data = loadmat(event_file)  # Load the .mat file
    event_times = event_data['StampCS'].flatten() / sampling_rate  # Convert event times

    # Save as .mat file, ensuring data is wrapped in a dictionary
    savemat(os.path.join(destination_folder, 'spike_times_mat.mat'), {'spike_times': spike_times})
    savemat(os.path.join(destination_folder, 'stim_times_mat.mat'), {'event_times': event_times})

def main(   mouse_name="Venice"
         ,  directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
         ,  extraction_folder="Extraction2Bin"
         ,  spike_file_name="spike_times.npy"
         ,  stim_file_name="EventStamps.mat"
         ,  sampling_rate=30000
         ,  switchSessions=False
         ):
    
    mouse_folder = f"{directory}{mouse_name}"
    if "ReserveMouse" in mouse_name:
        mouse_folder = mouse_folder.replace("MainFolder", "ReserveFolder")

    if switchSessions:
        mouse_folders = [os.path.join(mouse_folder, "SwitchSessionStitching")]
    else:
        # Get all folders with mouse_name in their name and without "copy"
        mouse_folders = [
            folder for folder in os.listdir(mouse_folder)
            if os.path.isdir(os.path.join(mouse_folder, folder)) 
            and mouse_name in folder 
            and "copy" not in folder
            and "Copy" not in folder
        ]
        mouse_folders.sort()    
        
    spike_file = os.path.join(mouse_folder, extraction_folder, spike_file_name)
    event_file = os.path.join(mouse_folder, extraction_folder, stim_file_name)
    
    load_and_save_as_mat(   spike_file
                         ,  event_file
                         ,  os.path.join(mouse_folder, extraction_folder)
                         ,  sampling_rate=sampling_rate
                         )

if __name__ == "__main__":
    main()
