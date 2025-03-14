#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 15:50:54 2025

@author: Ilse Klinkhamer
"""
import os
from get_dropbox_path import get_dropbox_path
from collections import defaultdict
import json


def makeAllMouseFoldersFile( directory=os.path.join(get_dropbox_path(),"C4_conversion/ZZZ_FolderCopy_IK_LeaveAlone/")
                           , save_dir=os.path.join(get_dropbox_path(), "ExperimentOutput/Ephys4Trace1/")
                           , include_switch_sessions=True):
        
    mice = get_mice()
    mice.sort()
    
    all_mouse_folders = defaultdict(list)
    
    for mouse in mice:
        
        dp_base = os.path.join(directory, mouse)        
        if os.path.exists(dp_base):
            if include_switch_sessions:                
                mouse_folders = [
                    folder for folder in os.listdir(dp_base)
                    if os.path.isdir(os.path.join(dp_base, folder)) 
                    and mouse in folder 
                    or "SwitchSessionStitching" in folder
                ]
            else:               
                mouse_folders = [
                    folder for folder in os.listdir(dp_base)
                    if os.path.isdir(os.path.join(dp_base, folder)) 
                    and mouse in folder 
                ] 
        
            mouse_folders.sort()     
            
            all_mouse_folders[mouse] = mouse_folders    
                

        with open(os.path.join(save_dir, "allMouseSessions.json"), "w") as file:
            json.dump(all_mouse_folders, file, indent=4)  # Pretty format with indentation
            
            
    return all_mouse_folders

def getAllMouseFolders(directory=os.path.join(get_dropbox_path(), "ExperimentOutput/Ephys4Trace1/")):
    with open(os.path.join(directory, "allMouseSessions.json"), "r") as all_mouse_folders_file:
        loaded_mouse_folders = json.load(all_mouse_folders_file)
    return loaded_mouse_folders

def getMouseFolders(mouse_name, directory=os.path.join(get_dropbox_path(), "ExperimentOutput/Ephys4Trace1/")):
    with open(os.path.join(directory, "allMouseSessions.json"), "r") as all_mouse_folders_file:
        loaded_mouse_folders = json.load(all_mouse_folders_file)
    mouse_folders = loaded_mouse_folders[mouse_name]
    return mouse_folders




def get_mice():
    """Returns a dictionary containing categorized mouse groups."""
    return ["ReserveMouse3", "ReserveMouse1", "ReserveMouse2", "ReserveMouse4", "ReserveMouse5", "Dallas", "Flint", "Greene"
            , "Houston", "Iowa", "Jackson", "Lincoln", "Newark", "Missouri", "Pittsburg", "Queens", "Orleans", "Reno"
            , "Seattle", "Yosemite", "Zachary", "Kyiv", "Istanbul", "Copenhagen", "Rotterdam", "Tallinn", "Quimper", "Porto"
            , "Lisbon", "Madrid", "Uppsala", "Venice", "Willemstad", "Zurich", "York", "Xanthi", "Ana1", "Ana2", "Ana3"
            , "Ana4", "Ana5"]

def main(mouseName = None):
    if mouseName == None:
        mouse_folders = getAllMouseFolders()  
    else:
        mouse_folders = getMouseFolders(mouseName)
    #makeAllMouseFoldersFile()
    return mouse_folders


if __name__ == "__main__":
    main()