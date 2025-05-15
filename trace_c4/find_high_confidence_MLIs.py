#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 15:11:15 2025

@author: Ilse Klinkhamer
"""

import pandas as pd
from pathlib import Path

def get_mouse_groups():
    """Returns a list of mouse names."""
    return [
        "ReserveMouse3", "Dallas", "Flint", "Greene", "Houston", "Iowa", "Jackson",
        "Lincoln", "Newark", "Missouri", "Pittsburg", "Queens", "Orleans",       
        "Reno", "Seattle", "Yosemite", "Zachary", "Kyiv", "Istanbul", "Copenhagen",
        "Rotterdam", "Tallinn", "Quimper", "Porto", "Lisbon", "Madrid",
        "Uppsala", "Venice", "Willemstad", "Zurich", "York", "Xanthi",
        "Ana1", "Ana2", "Ana3", "Ana4", "Ana5"
    ]

confidence_thres = 0.8
all_cells_with_high_confidence = []
max_confidence = 0.0  # to track the highest confidence score
cell_type = "PkC_ss"

base_folder = Path("/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/AnalysisOutput/c4 results stats/")  # change to your actual base path

for mouse in get_mouse_groups():
    mouse_folder = base_folder / mouse
    if not mouse_folder.exists():
        continue

    subfolders = [f for f in mouse_folder.iterdir() if f.is_dir() and mouse in f.name]

    for subfolder in subfolders:
        try:
            types_file = subfolder / "cluster_predicted_cell_type.tsv"
            prob_file = subfolder / "cluster_pred_probability.tsv"
    
            df_types = pd.read_csv(types_file, sep="\t")
            df_probs = pd.read_csv(prob_file, sep="\t")
    
            df = pd.concat([
                df_types["predicted_cell_type"],
                df_probs["pred_probability"],
                df_types["cluster_id"]
            ], axis=1)
    
            high_conf_cells = df[(df["predicted_cell_type"] == cell_type) & (df["pred_probability"] > confidence_thres)]
            all_cells = df[(df["predicted_cell_type"] == cell_type) & (df["pred_probability"] > 0)]
    
            if not high_conf_cells.empty:
                all_cells_with_high_confidence.append((mouse, subfolder.name, high_conf_cells))
    
            # Update max confidence
            current_max = all_cells["pred_probability"].max()
            if current_max > max_confidence:
                max_confidence = current_max
        except Exception as e:
            print(f"Failed to process {subfolder} for mouse {mouse}: {e}")
        

print(f"Highest MLI confidence found: {max_confidence:.4f}")


#Optional: flatten and combine
combined_df = pd.concat([df.assign(mouse=mouse, folder=folder) for mouse, folder, df in all_cells_with_high_confidence])
a = 1
