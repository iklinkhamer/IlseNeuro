#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 09:06:50 2025

@author: no1
"""

import os
from get_dropbox_path import get_dropbox_path

switch_sessions = True
directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
continuous_dir = os.path.join(get_dropbox_path(), "C4_conversion")
mouse_name = "ReserveMouse3"
continuous_path = os.path.join(continuous_dir, mouse_name)
data_folder = "continuous/Data_AP_LFP/"
data_file = "continuous.dat"

dp_base = os.path.join(directory, mouse_name)

if "ReserveMouse" in mouse_name:
    dp_base = dp_base.replace("MainFolder", "ReserveFolder")

# Get all folders with mouse_name in their name and without "copy"
mouse_folders = [
    folder for folder in os.listdir(continuous_path)
    if os.path.isdir(os.path.join(continuous_path, folder))
       and mouse_name in folder
       and "copy" not in folder
       and "Copy" not in folder
]
mouse_folders.sort()


if switch_sessions:
    switch_folder = os.path.join(continuous_path, "SwitchSessionStitching")
    if os.path.exists(switch_folder):
        mouse_folders.append(switch_folder)

for folder in mouse_folders:
    byte_size = os.path.getsize(os.path.join(continuous_path, folder, data_folder, data_file))
    byte_size_file_path = os.path.join(continuous_path, folder)
    print(byte_size)
    byte_size_file_save_location = os.path.join(dp_base, folder, "c4")
    if os.path.exists(byte_size_file_save_location):
        with open(f"{byte_size_file_save_location}/byte_size_continuous_file.txt", "w") as f:
            f.write(str(byte_size))    
    print(f"{mouse_name} data folder not synchronized to computer")
    with open(f"{byte_size_file_path}/byte_size_continuous_file.txt", "w") as f:
        f.write(str(byte_size))
            