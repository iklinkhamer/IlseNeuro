#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 11:32:44 2025

@author: Ilse Klinkhamer
"""
"""
import shutil
source_folder="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/ReserveFolder/ReserveMouse3/ReserveMouse3_20180708125233/c4/.NeuroPyxels/"
destination_folder="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/C4 Cache/ReserveMouse3/ReserveMouse3_20180708125233/.NeuroPyxels"
shutil.move(source_folder, destination_folder)
"""

import os
import shutil  
import sys
from get_dropbox_path import get_dropbox_path

def moveNeuroPyxelCacheFolders(   mouse_name
                                     ,  switch_sessions=True
                                     ,  directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
                                     ,cache_dir=os.path.join(os.path.dirname(get_dropbox_path().rstrip("/")), "C4 Cache")
                                     ):

    dp_base = os.path.join(directory, mouse_name)
    if "ReserveMouse" in mouse_name:
        dp_base = dp_base.replace("MainFolder", "ReserveFolder")
        
        
        
    if not os.path.exists(dp_base):
        print(f"{mouse_name} folder not synchronized to computer. skipping...")
        return

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
            switch_folder_name = "SwitchSessionStitching"                
            mouse_folders.append(switch_folder_name)    

    # Print the list of matching folders
    print(mouse_folders)
    
    for foldername in mouse_folders:
        
        print(f"Processing folder: {foldername}")
        
        # Perform actions on each folder here
        source_folder = os.path.join(dp_base, foldername, "c4",".NeuroPyxels")
        destination_folder = os.path.join(cache_dir, mouse_name, foldername)       
       

        if not os.path.exists(source_folder):
            print("Data folder not found, skipping.")
            continue          

        # Create the destination folder if it doesn't exist
        #os.makedirs(destination_folder, exist_ok=True)       
        if not os.path.exists(destination_folder):
            print("Destination folder not found, skipping")
            continue
        
        shutil.move(source_folder, destination_folder)
        print(f"Folder moved and renamed to: {destination_folder}")

                
                  
def main(   mouse_name=None
         ,  switch_sessions=True):
    
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
            mice = [mouse_name]
        else:
            try:
                mice = get_mice()
            except:
                print("Error: No mouse name provided")
                return
    else:
        mice = [mouse_name]
        
    for mouse_name in mice:
        moveNeuroPyxelCacheFolders(   mouse_name
                                         ,  switch_sessions=switch_sessions
                                         )    
    
    
def get_mice():
    """Returns a dictionary containing categorized mouse groups."""
    return ["ReserveMouse3"
            , "ReserveMouse1"
            , "ReserveMouse2"
            , "ReserveMouse4"
            , "ReserveMouse5"
            , "Dallas"
            , "Flint"
            , "Greene"
            , "Houston"
            , "Iowa"
            , "Jackson"
            , "Lincoln"
            , "Newark"
            , "Missouri"
            , "Pittsburg"
            , "Queens"
            , "Orleans"
            , "Reno"
            , "Seattle"
            , "Yosemite"
            , "Zachary"
            , "Kyiv"
            , "Istanbul"
            , "Copenhagen"
            , "Rotterdam"
            , "Tallinn"
            , "Quimper"
            , "Porto"
            , "Lisbon"
            , "Madrid"
            , "Uppsala"
            , "Venice"
            , "Willemstad"
            , "Zurich"
            , "York"
            , "Xanthi"
            , "Ana1"
            , "Ana2"
            , "Ana3"
            , "Ana4"
            , "Ana5"
            ]
        
if __name__ == "__main__":
    main()