#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:20:46 2025

@author: Ilse Klinkhamer
"""

import json
import os
import pandas as pd
import matplotlib.pyplot as plt
from get_dropbox_path import get_dropbox_path

def load_mice_metadata(directory):
    """Load mouse metadata."""
    
    json_file_path = os.path.join("/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/", "miceMetadata.json")
    
    with open(json_file_path, 'r') as f:
        mice_metadata = json.load(f)
    
    mouse_names = [mouse for mouse in mice_metadata.keys() if mouse not in ["ReserveMouse1", "ReserveMouse4"]]
    
    # Remove ReserveMouse1 and ReserveMouse4 from metadata
    mice_metadata = {mouse: data for mouse, data in mice_metadata.items() if mouse not in ["ReserveMouse1", "ReserveMouse4"]}
    
    return mouse_names, mice_metadata

def count_clusters_per_mouse(directory, mouse_names, mice_metadata, switch_sessions):
    """Count number of clusters per mouse and per session from cluster_group.tsv files."""
    
    cluster_counts = {}
    session_counts = {}
    group_counts = {}
    total_clusters_all_mice_not_switch_sessions = 0
    total_clusters_all_mice_switch_sessions = 0
    
    for mouse_name in mouse_names:
        dp_base = directory.replace("MainFolder", "ReserveFolder") if "ReserveMouse" in mouse_name else directory
        mouse_path = os.path.join(dp_base, mouse_name)

        if not os.path.exists(mouse_path):
            print(f"Warning: {mouse_path} does not exist. Skipping.")
            continue

        mouse_folders = [
            folder for folder in os.listdir(mouse_path)
            if os.path.isdir(os.path.join(mouse_path, folder))
            and mouse_name in folder
            and "copy" not in folder.lower()
        ]
        if switch_sessions:
            switch_folder = os.path.join(mouse_path, "SwitchSessionStitching")
            if os.path.exists(switch_folder):
                switch_folder_name = "SwitchSessionStitching"                
                mouse_folders.append(switch_folder_name)      
        
        mouse_folders.sort()
        total_clusters = 0
        session_counts.setdefault(mouse_name, {})
        total_single_clusters = 0
        total_wide_clusters = 0
        

        # Determine if it's a "Switch" mouse
        distType = mice_metadata[mouse_name]["distType"]
        is_switch = distType in ["DEUN", "Naive"]
        
        # Read firstWideSession.txt if it's a switch mouse
        first_wide_session = None
        if is_switch:
            wide_session_file = os.path.join(get_dropbox_path(), f"AnalysisOutput/c4 results stats/{mouse_name}/firstWideSession.txt")
            if os.path.isfile(wide_session_file):
                with open(wide_session_file, 'r') as f:
                    first_wide_session = f.read().strip()
        now_wide_sessions = False
        for folder in mouse_folders:
            c4_folder_path = os.path.join(mouse_path, folder, "c4")
            cluster_file = os.path.join(c4_folder_path, "cluster_group.tsv")

            if not os.path.isfile(cluster_file):
                print(f"Error: Missing {cluster_file}. Sync it with Dropbox first!")
                continue
            
            df = pd.read_csv(cluster_file, sep='\t')
            if {"cluster_id", "group"}.issubset(df.columns):
                session_clusters = df["cluster_id"].nunique()

                if folder==first_wide_session:
                    now_wide_sessions = True
                    
                if not now_wide_sessions:                        
                    total_single_clusters += session_clusters      
                    if distType == "DEUN":
                        if "Single" not in group_counts:
                            group_counts["Single"] = 0
                        group_counts["Single"] += session_clusters
                    elif distType == "Naive":
                        if "SingleNaive" not in group_counts:
                            group_counts["SingleNaive"] = 0
                        group_counts["SingleNaive"] += session_clusters
                elif now_wide_sessions and folder != switch_folder_name:
                    total_wide_clusters += session_clusters     
                    if distType == "DEUN":
                        if "WideNovice" not in group_counts:
                            group_counts["WideNovice"] = 0
                        group_counts["WideNovice"] += session_clusters
                    elif distType == "Naive":
                        if "WideNaive" not in group_counts:
                            group_counts["WideNaive"] = 0
                        group_counts["WideNaive"] += session_clusters
                                                                      
                if folder==switch_folder_name:                        
                    total_clusters_all_mice_switch_sessions += session_clusters
                else:
                    total_clusters += session_clusters
                    if not mouse_name.startswith("Ana"):
                        total_clusters_all_mice_not_switch_sessions += session_clusters   
                    if distType not in group_counts:
                        group_counts[distType] = 0
                    group_counts[distType] += session_clusters
                    
                session_counts[mouse_name][folder] = session_clusters
               
        
        cluster_counts[mouse_name] = total_clusters
        if is_switch:
            cluster_counts[f"{mouse_name}-Single"] = total_single_clusters
            cluster_counts[f"{mouse_name}-Wide"] = total_wide_clusters
    cluster_counts["Total"] = total_clusters_all_mice_not_switch_sessions
    cluster_counts["Total Switch Sessions"] = total_clusters_all_mice_switch_sessions
    
    return cluster_counts, session_counts, group_counts

def save_results(directory, cluster_counts, session_counts, group_counts):
    """Save cluster count results to TSV files and generate bar plots."""
    
    cluster_output_file = os.path.join(directory, "cluster_summary.tsv")
    session_output_file = os.path.join(directory, "session_summary.tsv")
    group_output_file = os.path.join(directory, "group_summary.tsv")
    plot_output_file = os.path.join(directory, "cluster_counts.png")
    
    # Save cluster counts
    with open(cluster_output_file, 'w') as f:
        f.write("Mouse\tClusterCount\n")
        for mouse, count in cluster_counts.items():
            f.write(f"{mouse}\t{count}\n")
    
    print(f"Cluster counts saved to {cluster_output_file}")
    
    # Save session counts
    with open(session_output_file, 'w') as f:
        f.write("Mouse\tSession\tClusterCount\n")
        for mouse, sessions in session_counts.items():
            for session, count in sessions.items():
                f.write(f"{mouse}\t{session}\t{count}\n")
    
    print(f"Session counts saved to {session_output_file}")
    
    # Save session counts
    with open(group_output_file, 'w') as f:
        f.write("Group\tClusterCount\n")
        for group, count in group_counts.items():
            f.write(f"{group}\t{count}\n")
    
    print(f"Group counts saved to {group_output_file}")
    
    bar_plot_data = cluster_counts
    bar_plot_data = {mouse: data for mouse, data in bar_plot_data.items() if mouse not in ["Total", "Total Switch Sessions"]}
    # Plot bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(bar_plot_data.keys(), bar_plot_data.values(), color='skyblue')
    plt.xlabel("Mouse")
    plt.ylabel("Cluster Count")
    plt.title("Cluster Counts per Mouse")
    plt.xticks(rotation=90)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.savefig(plot_output_file)
    plt.show()
    
    print(f"Cluster count plot saved to {plot_output_file}")

def main():
    directory = os.path.join(get_dropbox_path(), "ExperimentOutput/Ephys4Trace1/MainFolder/")
    mouse_names, mice_metadata = load_mice_metadata(directory)
    cluster_counts, session_counts, group_counts = count_clusters_per_mouse(directory, mouse_names, mice_metadata, switch_sessions=True)
    save_results(directory, cluster_counts, session_counts, group_counts)

if __name__ == "__main__":
    main()



            