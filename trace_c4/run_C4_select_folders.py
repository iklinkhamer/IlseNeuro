#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:45:07 2025

@author: Ilse Klinkhamer
"""

import os
import re
import shutil
import OpenEphys
from tkinter import Tk, filedialog
import sys
import npyx
from get_dropbox_path import get_dropbox_path
from os import path
from pathlib import Path
from datetime import datetime
import pandas as pd
import tkinter as tk
import numpy as np

"""
def select_folder(prompt="Select a folder"):
    root = Tk()
    root.withdraw()  # Hide the main window
    folder_path = filedialog.askdirectory(title=prompt)
    root.destroy()
    return folder_path
"""

LAST_FOLDER_FILE = Path.home() / ".last_folder.txt"  # Hidden file in home dir

def select_folder(prompt="Select a folder"):
    # Load last folder if available
    initial_dir = LAST_FOLDER_FILE.read_text() if LAST_FOLDER_FILE.exists() else str(Path.home())

    # Open folder selection dialog
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(title=prompt, initialdir=initial_dir)

    # Save selected folder for next time
    if folder_path:
        LAST_FOLDER_FILE.write_text(folder_path)

    return folder_path if folder_path else None

def find_session_and_mouse_name(path):
    path = Path(path).resolve()

    for parent in [path] + list(path.parents):
        name = parent.name

        # Remove trailing session number if present
        base = name.rsplit('_', 1)[0] if name.rsplit('_', 1)[-1].isdigit() else name

        # Extract digits only
        digits = ''.join(c for c in base if c.isdigit())

        # Try parsing with different datetime formats
        for fmt in ("%Y%m%d%H%M%S", "%Y%m%d%H%M", "%Y%m%d"):
            try:
                datetime.strptime(digits, fmt)
                mouse_name = parent.parent.name  # Get the parent folder name
                return name, mouse_name
            except ValueError:
                continue

    return None, None  # No session-like folder found

def filter_units_by_confidence_ratio(phy_compatible_c4_folder, contamination_ratio=0.99, threshold_old=0, threshold_filter_new=1.5):
    base_results_folder = os.path.join(phy_compatible_c4_folder, f"c4_results_fpfnThreshold_{contamination_ratio}_confidenceRatio_{threshold_old}")
        #'/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/Ilse/ephys/Ilo_Cbx/I_Cbx_2025-02-25_10-33-00_7/c4/c4_results_fpfnThreshold_0.1_confidenceRatio_0'
    # Find the cluster confidence file in base_results_folder
    confidence_file = os.path.join(base_results_folder, 'cluster_confidence_ratio.tsv')
    if not os.path.exists(confidence_file):
        print(f"{confidence_file} does not exist.")
        return
        #raise FileNotFoundError(f"{confidence_file} does not exist.")

    # Load confidence ratio data
    conf_df = pd.read_csv(confidence_file, sep='\t')

    # Filter units with confidence_ratio > threshold
    filtered_df = conf_df[conf_df['confidence_ratio'] > threshold_filter_new]

    # Prepare new results folder path
    folder_name = os.path.basename(base_results_folder)
    parent_folder = os.path.dirname(base_results_folder)
    new_folder_name = folder_name.replace(f'confidenceRatio_{threshold_old}', f'confidenceRatio_{threshold_filter_new}')
    new_results_folder = os.path.join(parent_folder, new_folder_name)
    os.makedirs(new_results_folder, exist_ok=True)

    # Prepare cell_type_classification folders
    old_classif_folder = os.path.join(base_results_folder, 'cell_type_classification')
    new_classif_folder = os.path.join(new_results_folder, 'cell_type_classification')
    os.makedirs(new_classif_folder, exist_ok=True)

    for cluster_id in filtered_df['cluster_id']:
        eps_file = os.path.join(old_classif_folder, f'unit_{cluster_id}_cell_type_predictions.eps')
        pdf_file = os.path.join(old_classif_folder, f'unit_{cluster_id}_cell_type_predictions.pdf')

        # Copy the files if they exist
        if os.path.exists(eps_file):
            shutil.copy2(eps_file, new_classif_folder)
        if os.path.exists(pdf_file):
            shutil.copy2(pdf_file, new_classif_folder)

    # Define the files to filter
    files = [
        ("cluster_predicted_cell_type.tsv", "predicted_cell_type"),
        ("cluster_confidence_ratio.tsv", "confidence_ratio"),
        ("cluster_model_votes.tsv", "model_votes"),
        ("cluster_pred_probability.tsv", "pred_probability")
    ]

    # Filter each file using the same cluster IDs
    for filename, column in files:
        file_path = os.path.join(base_results_folder, filename)
        if not os.path.exists(file_path):
            print(f"Warning: {filename} not found, skipping.")
            continue

        df = pd.read_csv(file_path, sep='\t')
        filtered_df_file = df[df['cluster_id'].isin(filtered_df['cluster_id'])]
        filtered_file_path = os.path.join(new_results_folder, filename)
        filtered_df_file.to_csv(filtered_file_path, sep='\t', index=False)

    print(f"Filtered results saved in: {new_results_folder}")

def run_c4( phy_compatible_c4_folder
          , classify_again=True
          , contamination_ratio=0.1
          , confidence_ratio_threshold=0
          , cache_dir = None
          , data_dir = None
          , filter_spikes = False
          , phy_c4_folder = "c4"
          ):

    from npyx.c4.predict_cell_types import run_cell_types_classifier

    if cache_dir is None:
        cache_dir = os.path.join(os.path.dirname(get_dropbox_path().rstrip("/")), "C4 Cache")

    if data_dir is None:
        data_dir = select_folder("Select the folder with the continuous folder for raw binary data in it")

    print(f"Processing folder: {phy_compatible_c4_folder}")
    save_path = os.path.join(phy_compatible_c4_folder, f"c4_results_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}")
    results_file_path = os.path.join(save_path, "cluster_predicted_cell_type.tsv")

    if not os.path.exists(phy_compatible_c4_folder):
        print(f"Folder {phy_compatible_c4_folder} does not exist, skipping c4 analysis.")
        return  # Skip to the next iteration if the folder doesn't exist
    if not os.path.exists(f"{phy_compatible_c4_folder}/params.py"):
        print(f"params.py folder not found in folder {phy_compatible_c4_folder}, skipping c4 analysis...")
        return
    if not classify_again and os.path.isfile(results_file_path):
        print(f"Session {phy_compatible_c4_folder} has already been classified and classify again is false, skipping c4 analysis.")
        return

    os.makedirs(save_path, exist_ok=True)
    session_folder, mouse_name = find_session_and_mouse_name(phy_compatible_c4_folder)
    cache_path = os.path.join(cache_dir, mouse_name, session_folder, ".NeuroPyxels")
    os.makedirs(cache_path, exist_ok=True)

    if not os.path.exists(os.path.join(data_dir, "continuous")):
        data_dir = os.path.join(data_dir, phy_c4_folder)
        if not os.path.exists(os.path.join(data_dir, "continuous")):
            print("continuous folder not found in data path folder, so please check.")

    run_cell_types_classifier(phy_compatible_c4_folder, quality='all', parallel=False, fp_threshold=contamination_ratio,
                              fn_threshold=contamination_ratio, threshold=confidence_ratio_threshold,
                              save_path=save_path, cache_path=cache_path, dat_path=data_dir, filter_spikes=filter_spikes)


def create_cluster_group_file(phy_compatible_c4_folder):
    ks_label_file = os.path.join(phy_compatible_c4_folder, "cluster_KSLabel.tsv")
    spike_clusters_file = os.path.join(phy_compatible_c4_folder, "spike_clusters.npy")
    output_file = os.path.join(phy_compatible_c4_folder, "cluster_group.tsv")

    if not os.path.exists(ks_label_file):
        return
        #raise FileNotFoundError(f"{ks_label_file} not found.")
    if not os.path.exists(spike_clusters_file):
        return
        #raise FileNotFoundError(f"{spike_clusters_file} not found.")

    # Load KSLabel info
    ks_df = pd.read_csv(ks_label_file, sep='\t')

    # Load spike_clusters and get unique cluster IDs
    spike_clusters = np.load(spike_clusters_file)
    unique_spike_clusters = pd.DataFrame({'cluster_id': np.unique(spike_clusters)})

    # Combine all cluster IDs from both sources
    all_cluster_ids = pd.DataFrame({'cluster_id': pd.concat([
        ks_df['cluster_id'], unique_spike_clusters['cluster_id']
    ]).drop_duplicates().sort_values().reset_index(drop=True)})

    # Merge and fill missing KSLabels with 'unsorted'
    merged_df = all_cluster_ids.merge(ks_df, on='cluster_id', how='left')
    merged_df['KSLabel'] = merged_df['KSLabel'].fillna('unsorted')

    # Save to cluster_group.tsv
    merged_df.to_csv(output_file, sep='\t', index=False)

    print(f"cluster_group.tsv written to: {output_file}")


def main(phy_compatible_c4_folder=os.path.join(select_folder("Select the phy compatible c4 input folder")),
         classify_again=True,
         contamination_ratio=0.99,
         confidence_ratio_threshold_c4_run=0,
         confidence_ratio_threshold_results_filter=1.5):
    # Find all valid session subfolders
    subfolders = []
    for f in os.listdir(phy_compatible_c4_folder):
        full_path = os.path.join(phy_compatible_c4_folder, f)
        if os.path.isdir(full_path):
            session_name, mouse_name = find_session_and_mouse_name(full_path)
            if session_name is not None and mouse_name is not None:
                if "c4" not in full_path:
                    full_path = os.path.join(full_path, "c4")
                subfolders.append(full_path)
                subfolders.sort()

    if not subfolders:
        subfolders = phy_compatible_c4_folder

    for session_folder in subfolders:
        print(f"\n--- Processing: {session_folder} ---")

        if not os.path.exists(os.path.join(session_folder, "cluster_group.tsv")):
            create_cluster_group_file(session_folder)

        #try:
        run_c4(session_folder, classify_again, contamination_ratio, confidence_ratio_threshold_c4_run, data_dir=session_folder, filter_spikes=False)
        #except:
        #    print("Crashed, likely due to: ValueError: No units were found with the provided parameter choices after quality checks.")

        filter_units_by_confidence_ratio(session_folder,
                                         contamination_ratio=contamination_ratio,
                                         threshold_old=confidence_ratio_threshold_c4_run,
                                         threshold_filter_new=confidence_ratio_threshold_results_filter)

if __name__ == "__main__":
    main()
