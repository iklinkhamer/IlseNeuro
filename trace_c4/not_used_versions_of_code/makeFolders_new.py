#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 18:42:32 2025

@author: Ilse Klinkhamer
"""
import os
import shutil  
from get_dropbox_path import get_dropbox_path


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
            , "Amsterdam"]

def get_switch_mice():
    """Returns a dictionary containing categorized mouse groups."""
    return [ "Houston"
            , "Iowa"
            , "Jackson"
            , "Lincoln"            
            , "Pittsburg"
            , "Orleans"
            , "Ana1"
            , "Ana2"
            , "Ana3"
            , "Ana4"
            , "Ana5"]

def main(   mouse_name="-"
        ,  switch_sessions=True
        ,  directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
        ,  destination_folder="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/C4_conversion/"
        ):
    if mouse_name is None:
        mice = get_mice()
    else:
        mice = [mouse_name]
        
        
    switchMice = get_switch_mice()
 
    for mouse_name in mice:
        dp_base = os.path.join(directory, mouse_name)
        if "ReserveMouse" in mouse_name:
            dp_base = dp_base.replace("MainFolder", "ReserveFolder")
        if os.path.exists(dp_base):
            # Get all folders with mouse_name in their name and without "copy"
            mouse_folders = [
                folder for folder in os.listdir(dp_base)
                if os.path.isdir(os.path.join(dp_base, folder))
                   and mouse_name in folder
                   and "copy" not in folder
                   and "Copy" not in folder
            ]
            mouse_folders.sort()
        else:
            print(f"Folder {mouse_name} not found, skipping copying files to Analysis Output...")
            continue
        

        switch_folder_name = "SwitchSessionStitching"                
        mouse_folders.append(switch_folder_name)                
                
        print(mouse_folders)
        
        for foldername in mouse_folders:
            results_folder_path = os.path.join(dp_base, foldername)
          
            if os.path.exists(os.path.join(results_folder_path)):            
                new_folder = os.path.join(destination_folder,mouse_name,foldername)
                os.makedirs(new_folder, exist_ok=True)
            
            else:
                print(f"{results_folder_path} not found, skipping copying this file in folder {foldername}...")
            
        
if __name__ == "__main__":
    main()