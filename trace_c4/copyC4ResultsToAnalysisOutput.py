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
            , "Ana5"]

def main(   mouse_name=None
        ,  switch_sessions=True
        ,  directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
        ,  destination_folder="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/AnalysisOutput/c4 results stats/"
        ):
    if mouse_name is None:
        mice = get_mice()
    else:
        mice = [mouse_name]
 
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
        
        if switch_sessions:
            switch_folder = os.path.join(dp_base, "SwitchSessionStitching")
            if os.path.exists(switch_folder):
                mouse_folders.append(switch_folder)
                
                
            print(mouse_folders)
        
        for foldername in mouse_folders:
            c4_results_folder_path = os.path.join(dp_base, foldername, "c4", "c4_results_fpfnThreshold_0.1_confidenceRatio_1.5")
            c4_results_file = 'cluster_predicted_cell_type.tsv'
            c4_results_files = ["cluster_confidence_ratio.tsv", "cluster_model_votes.tsv", "cluster_predicted_cell_type.tsv", "cluster_pred_probability.tsv", f"{foldername}_discharge_stats.tsv"]
            
            if os.path.exists(os.path.join(c4_results_folder_path, c4_results_file)):            
                new_folder = os.path.join(destination_folder,mouse_name,foldername)
                os.makedirs(new_folder, exist_ok=True)
                
                for file in c4_results_files:
                    source_file = os.path.join(c4_results_folder_path, file)
                    destination_file = os.path.join(new_folder, file)
                    if os.path.exists(source_file):
                        shutil.copy(source_file, destination_file)  # Copy file
                        print(f"File copied to: {destination_file}")
                    else:
                        print(f"Source file does not exist: {source_file}")
            else:
                print(f"c4 results {mouse_name} not found, skipping copying this file in folder {foldername}...")
            
        
if __name__ == "__main__":
    main()