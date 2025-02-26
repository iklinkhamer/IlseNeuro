#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:10:03 2025

@author: Ilse Klinkhamer
"""



import os
import sys
from collections import defaultdict
from npyx.spk_t import trn, isi, inst_cv2, mean_firing_rate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from get_dropbox_path import get_dropbox_path

def save_boxplots(all_sessions_data, save_path, title, show_outliers=False):
    """Generate and save boxplots with scatter overlay for ISI, CV2, and firing rate.
    
    Arguments:
        - all_sessions_data: dict, structured spike metrics
        - save_path: str, where to save the plot
        - title: str, plot title
        - show_outliers: bool, if True, boxplot shows outliers, but scatter never does
    """
    metrics = ["Mean Firing Rate", "Mean CV", "Mean CV2", "Median ISI"]
    fixed_order = ["PkC_ss", "PkC_cs", "MLI", "GoC", "MFB"]  # Fixed order (if present)
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    
    color_map = {
        "PkC_ss": "blue",
        "PkC_cs": "black",
        "MLI": "pink",
        "GoC": "green",
        "MFB": "red"
    }

    for i, metric in enumerate(metrics):
        # Extract data, keeping only present cell types
        data = {cell_type: [entry[i+1] for entry in values] for cell_type, values in all_sessions_data.items()}
        df = pd.DataFrame({k: pd.Series(v) for k, v in data.items() if len(v) > 0})

        # Keep only present cell types but in the correct order
        present_cell_types = [ct for ct in fixed_order if ct in df.columns]
        df = df[present_cell_types]  # Reorder columns, drop missing ones

        # Identify outliers using IQR
        outliers_dict = {}
        filtered_data = {}

        for cell_type in present_cell_types:
            y = df[cell_type].dropna()

            # Compute IQR bounds
            Q1 = y.quantile(0.25)
            Q3 = y.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR

            # Identify and store outliers
            outliers = y[(y < lower_bound) | (y > upper_bound)]
            outliers_dict[cell_type] = outliers.tolist()

            # Keep only non-outlier data
            filtered_data[cell_type] = y[(y >= lower_bound) & (y <= upper_bound)]

        # Print detected outliers
        for cell_type, outliers in outliers_dict.items():
            if outliers:
                print(f"Outliers in {metric} ({cell_type}): {outliers}")

        # Boxplot with only present cell types
        box_data = [df[col].dropna() if show_outliers else filtered_data.get(col, []) for col in present_cell_types]
        box = axes[i].boxplot(box_data, patch_artist=True, labels=present_cell_types, showfliers=show_outliers)

        # Modify boxplot colors
        for patch, cell_type in zip(box["boxes"], present_cell_types):
            patch.set_facecolor(color_map.get(cell_type, "gray"))
            patch.set_alpha(0.5)

        # Modify median line to be thick and black
        for median in box["medians"]:
            median.set_color("black")
            median.set_linewidth(2.5)

        # Overlay scatter points (only non-outliers)
        for j, cell_type in enumerate(present_cell_types):
            y = filtered_data[cell_type]
            x = np.random.normal(j + 1, 0.04, size=len(y))  # Add jitter
            axes[i].scatter(x, y, facecolors=color_map.get(cell_type, "gray"), edgecolors="none", s=30, alpha=1)

        axes[i].set_title(metric)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
    print(f"Boxplots saved to {save_path}")




def get_discharge_statistics(mouse_name, switch_sessions=False, contamination_ratio=0.1, confidence_ratio_threshold=2, directory=os.path.join(get_dropbox_path(),"/ExperimentOutput/Ephys4Trace1/MainFolder/"), skip_without_continuous=True):
    
    dp_base = os.path.join(directory, mouse_name)
    if "ReserveMouse" in mouse_name:
        dp_base = dp_base.replace("MainFolder", "ReserveFolder")

    mouse_folders = [
        folder for folder in os.listdir(dp_base)
        if os.path.isdir(os.path.join(dp_base, folder)) and mouse_name in folder and "copy" not in folder.lower()
    ]
    mouse_folders.sort()

    if switch_sessions:
        switch_folder = os.path.join(dp_base, "SwitchSessionStitching")
        if os.path.exists(switch_folder):
            mouse_folders.append(switch_folder)

    phy_folder = "c4"
    all_sessions_data = defaultdict(list)  # Store data by cell type across all sessions
    all_sessions_data_stats = []
    
    for sess in mouse_folders:
        print(f"Processing session: {sess}")
        dp = os.path.join(dp_base, sess, phy_folder)
        save_path = os.path.join(dp, f"c4_results_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}")
        os.makedirs(save_path, exist_ok=True)

        unit_file = os.path.join(save_path, "cluster_predicted_cell_type.tsv")
        if not os.path.exists(unit_file):
            continue  # Skip if the file doesn't exist

        df_units = pd.read_csv(unit_file, sep="\t", usecols=["cluster_id", "predicted_cell_type"])
        df_units = df_units.dropna()

        session_data = defaultdict(list)  # Store per-session data by cell type
        session_data_stats = []  # Store per-session data
        for _, row in df_units.iterrows():
            cluster_id, cell_type = int(row["cluster_id"]), row["predicted_cell_type"]
            t = trn(dp, cluster_id)
            ISIs = isi(dp, cluster_id)
            instant_cv2 = inst_cv2(t)
            av_firing_rate = mean_firing_rate(t)
            cv = compute_cv(t)

            data_entry = [cluster_id, av_firing_rate, np.mean(cv), np.mean(instant_cv2), np.median(ISIs)]
            session_data[cell_type].append(data_entry)
            all_sessions_data[cell_type].append(data_entry)
            data_entry = [sess, cluster_id, cell_type, av_firing_rate, np.mean(cv), np.mean(instant_cv2), np.median(ISIs)]
            session_data_stats.append(data_entry)
            all_sessions_data_stats.append(data_entry)

        # Save per-session TSV
        df_session = pd.DataFrame(session_data_stats, columns=["Session", "Cluster ID", "Cell Type", "Mean Firing Rate", "Mean CV", "Mean CV2", "Median ISI"])
        session_tsv = os.path.join(save_path, f"{sess}_discharge_stats.tsv")
        df_session.to_csv(session_tsv, sep="\t", index=False)
        print(f"Session statistics saved to {session_tsv}")
        
        # Save session boxplots
        session_plot = os.path.join(save_path, f"{sess}_discharge_boxplots.png")
        save_boxplots(session_data, session_plot, f"Discharge Statistics - {sess}")

    # Save overall statistics
    df_overall = pd.DataFrame(all_sessions_data_stats, columns=["Session", "Cluster ID", "Cell Type", "Mean Firing Rate", "Mean CV", "Mean CV2", "Median ISI"])
    overall_tsv = os.path.join(dp_base, f"{mouse_name}_overall_discharge_stats.tsv")
    df_overall.to_csv(overall_tsv, sep="\t", index=False)
    print(f"Overall statistics saved to {overall_tsv}")

    # Save combined boxplots
    overall_plot = os.path.join(dp_base, f"{mouse_name}_overall_discharge_boxplots.png")
    save_boxplots(all_sessions_data, overall_plot, "Overall Discharge Statistics")
    return all_sessions_data_stats
    
def compute_cv(t):
    """
    Compute the coefficient of variation (CV) of interspike intervals.

    Arguments:
        - t: (nspikes,) np.array, spike times in any unit

    Returns:
        - cv: float, coefficient of variation of ISIs
    """
    if len(t) < 2:
        return np.nan  # CV is undefined for fewer than 2 spikes

    isis = np.diff(t)  # Compute interspike intervals
    return np.std(isis) / np.mean(isis)  # CV formula
    
def main(mouse_name=None, switch_sessions=True, contamination_ratio=0.1, confidence_ratio_threshold=1.5):
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
        else:
            print("Error: No mouse name provided")
            sys.exit(1)
    discharge_statistics = get_discharge_statistics(mouse_name, switch_sessions, contamination_ratio, confidence_ratio_threshold)
    return discharge_statistics

if __name__ == "__main__":
    main()