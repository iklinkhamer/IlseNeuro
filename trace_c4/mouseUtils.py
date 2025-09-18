#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 15:50:54 2025

@author: Ilse Klinkhamer
"""
import os
from getBayesLabDropboxRoot import get_dropbox_path
from collections import defaultdict
import json

def makeMouseMetaDataFile(
    directory=os.path.join(get_dropbox_path(), "TraceExperiments/C4_conversion/ZZZ_FolderCopy_IK_LeaveAlone/"),
    save_dir=os.path.join(get_dropbox_path(), "TraceExperiments/ExperimentOutput/Ephys4Trace1/"),
    include_switch_sessions=True
):
    mice = get_mice()
    mice.sort()

    all_mouse_folders = defaultdict(list)
    metadata = defaultdict(dict)  # Corrected to use dict instead of list

    for mouse in mice:
        dp_base = os.path.join(directory, mouse)
        if os.path.exists(dp_base):
            if include_switch_sessions:
                mouse_folders = [
                    folder for folder in os.listdir(dp_base)
                    if os.path.isdir(os.path.join(dp_base, folder)) and (mouse in folder or "SwitchSessionStitching" in folder)
                ]
            else:
                mouse_folders = [
                    folder for folder in os.listdir(dp_base)
                    if os.path.isdir(os.path.join(dp_base, folder)) and mouse in folder
                ]

            mouse_folders.sort()

            all_mouse_folders[mouse] = mouse_folders
            metadata[mouse]["folders"] = mouse_folders
            metadata[mouse]["distType"] = next(
                (group for group, mice_list in get_mouse_groups().items() if mouse in mice_list),
                "Unknown"
            )
            metadata[mouse]["ephysSessions"] = True
            metadata[mouse]["switchOccurredDuringEphys"] = any("SwitchSessionStitching" in folder for folder in mouse_folders)


    # Save metadata
    os.makedirs(save_dir, exist_ok=True)
    with open(os.path.join(save_dir, "miceMetadata.json"), "w") as file:
        json.dump(metadata, file, indent=4)

    return all_mouse_folders

def makeAllMouseFoldersFile(
    directory=os.path.join(get_dropbox_path(), "TraceExperiments/C4_conversion/ZZZ_FolderCopy_IK_LeaveAlone/"),
    save_dir=os.path.join(get_dropbox_path(), "TraceExperiments/ExperimentOutput/Ephys4Trace1/"),
    include_switch_sessions=True
):
    mice = get_mice()
    mice.sort()

    all_mouse_folders = defaultdict(list)

    for mouse in mice:
        dp_base = os.path.join(directory, mouse)
        if os.path.exists(dp_base):
            if include_switch_sessions:
                mouse_folders = [
                    folder for folder in os.listdir(dp_base)
                    if os.path.isdir(os.path.join(dp_base, folder)) and (mouse in folder or "SwitchSessionStitching" in folder)
                ]
            else:
                mouse_folders = [
                    folder for folder in os.listdir(dp_base)
                    if os.path.isdir(os.path.join(dp_base, folder)) and mouse in folder
                ]

            mouse_folders.sort()
            all_mouse_folders[mouse] = mouse_folders

    os.makedirs(save_dir, exist_ok=True)
    with open(os.path.join(save_dir, "allMouseSessions.json"), "w") as file:
        json.dump(all_mouse_folders, file, indent=4)

    return all_mouse_folders

def getAllMouseFolders(directory=os.path.join(get_dropbox_path(), "TraceExperiments/ExperimentOutput/Ephys4Trace1/")):
    file_path = os.path.join(directory, "allMouseSessions.json")
    if not os.path.exists(file_path):
        return {}
    with open(file_path, "r") as file:
        return json.load(file)

def getMouseFolders(mouse_name, directory=os.path.join(get_dropbox_path(), "TraceExperiments/ExperimentOutput/Ephys4Trace1/")):
    all_mouse_folders = getAllMouseFolders(directory)
    return all_mouse_folders.get(mouse_name, [])

def get_mice():
    """Returns a list of all mice."""
    return [
        "ReserveMouse3", "ReserveMouse1", "ReserveMouse2", "ReserveMouse4", "ReserveMouse5", "Dallas", "Flint", "Greene",
        "Houston", "Iowa", "Jackson", "Lincoln", "Newark", "Missouri", "Pittsburg", "Queens", "Orleans", "Reno",
        "Seattle", "Yosemite", "Zachary", "Kyiv", "Istanbul", "Copenhagen", "Rotterdam", "Tallinn", "Quimper", "Porto",
        "Lisbon", "Madrid", "Uppsala", "Venice", "Willemstad", "Zurich", "York", "Xanthi", "Ana1", "Ana2", "Ana3",
        "Ana4", "Ana5"
    ]

def get_mouse_groups():
    """Returns a dictionary containing categorized mouse groups."""
    return {
        "DEUN": [
            "ReserveMouse3", "ReserveMouse1", "ReserveMouse2", "ReserveMouse4", "ReserveMouse5", "Dallas", "Flint", "Greene", "Houston", "Iowa", "Jackson",
            "Lincoln", "Newark", "Missouri", "Pittsburg", "Queens", "Orleans"
        ],
        "Wide": ["Reno", "Seattle", "Yosemite", "Zachary", "Kyiv", "Istanbul", "Copenhagen"],
        "Narrow": ["Rotterdam", "Tallinn", "Quimper", "Porto", "Lisbon", "Madrid"],
        "Bimodal": ["Uppsala", "Venice", "Willemstad", "Zurich", "York", "Xanthi"],
        "Naive": ["Ana1", "Ana2", "Ana3", "Ana4", "Ana5"]
    }

def main(mouse_name=None):
    if mouse_name:
        return getMouseFolders(mouse_name)
    makeMouseMetaDataFile()
    return getAllMouseFolders()

if __name__ == "__main__":
    main()
