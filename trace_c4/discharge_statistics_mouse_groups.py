#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:37:44 2025

@author: Ilse Klinkhamer
"""

import os
import pandas as pd
from collections import defaultdict
import getting_mouse_groups
import discharge_statistics

def compute_discharge_statistics_by_mouse_type(directory, save_directory):
    """
    Computes and saves discharge statistics for each mouse type.

    Arguments:
        - directory: str, base path where mouse data is stored.
        - save_directory: str, path to save the compiled statistics.

    Returns:
        - discharge_stats: dict, structured as {mouse_type: {mouse_name: data}}.
    """
    mouse_groups = getting_mouse_groups.main()
    discharge_stats = defaultdict(dict)  # Store per mouse type

    os.makedirs(save_directory, exist_ok=True)  # Ensure save directory exists

    for mouse_type, mice in mouse_groups.items():
        print(f"Processing category: {mouse_type}")
        mouse_type_stats = []

        for mouse in mice:
            mouse_path = os.path.join(directory, mouse)
            if not os.path.exists(mouse_path):
                print(f"Skipping {mouse}, data directory not found.")
                continue

            print(f"Computing statistics for {mouse} ({mouse_type})...")
            stats_file = os.path.join(mouse_path, f"{mouse}_overall_discharge_stats.tsv")

            # Compute stats if missing
            if not os.path.exists(stats_file):
                discharge_statistics.main(mouse)

            # Load computed stats
            if os.path.exists(stats_file):
                df = pd.read_csv(stats_file, sep="\t")
                discharge_stats[mouse_type][mouse] = df
                mouse_type_stats.append(df)

        # Save combined stats for this mouse type
        if mouse_type_stats:
            combined_df = pd.concat(mouse_type_stats, ignore_index=True)
            save_path = os.path.join(save_directory, f"{mouse_type}_discharge_stats.tsv")
            combined_df.to_csv(save_path, sep="\t", index=False)
            print(f"Saved {mouse_type} statistics to {save_path}")

    return discharge_stats

# Example usage
base_directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/"
save_directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/test"

stats_by_type = compute_discharge_statistics_by_mouse_type(base_directory, save_directory)

# Access statistics for a specific mouse type
print(stats_by_type["Switch"])  # Example: print data for "Switch" mice

