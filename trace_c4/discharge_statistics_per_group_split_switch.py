#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 12:58:45 2025

@author: Ilse Klinkhamer
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from get_dropbox_path import get_dropbox_path
from plot_utils_IK import c4_colors_rgb, lighten, normalize_RGB_dict
from collections import defaultdict
from scipy.stats import zscore

def get_mouse_groups():
    """Returns a dictionary containing categorized mouse groups."""
    return {
        "Switch": [
            "ReserveMouse3", "Dallas", "Flint", "Greene", "Houston", "Iowa", "Jackson",
            "Lincoln", "Newark", "Missouri", "Pittsburg", "Queens", "Orleans"
        ],
        "WideExperts": ["Reno", "Seattle", "Yosemite", "Zachary", "Kyiv", "Istanbul", "Copenhagen"],
        "Narrow": ["Rotterdam", "Tallinn", "Quimper", "Porto", "Lisbon", "Madrid"],
        "Bimodal": ["Uppsala", "Venice", "Willemstad", "Zurich", "York", "Xanthi"],
        "Naive": ["Ana1", "Ana2", "Ana3", "Ana4", "Ana5"]
    }

def read_tsv(file_path):
    """Reads a TSV file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', index_col=0)

def save_boxplots(all_sessions_data, save_path, title, show_outliers=False):
    """Generate and save a single boxplot figure with all metrics and individual data points.

    Arguments:
        - all_sessions_data: dict, structured spike metrics
        - save_path: str, where to save the plot
        - title: str, plot title
        - show_outliers: bool, if True, boxplot shows outliers, but scatter never does
    """
    metrics = ["Mean Firing Rate", "Mean CV", "Mean CV2", "Median ISI"]
    fixed_order = ["PkC_ss", "PkC_cs", "MLI", "GoC", "MFB"]  # Fixed order (if present)
    
    fig, ax = plt.subplots(figsize=(12, 6))  # Single axis for all plots
    
    color_map = c4_colors_rgb() 
    color_map['PkC_cs'] = lighten(color_map['PkC_cs'])
    color_map = normalize_RGB_dict(color_map)

    all_box_data = []
    positions = []
    color_list = []

    scatter_x = []
    scatter_y = []
    scatter_colors = []

    metric_positions = []  # To store x positions for metric labels

    # Process each metric
    for i, metric in enumerate(metrics):
        data = {cell_type: [entry[i+1] for entry in values] for cell_type, values in all_sessions_data.items()}
        df = pd.DataFrame({k: pd.Series(v) for k, v in data.items() if len(v) > 0})

        present_cell_types = [ct for ct in fixed_order if ct in df.columns]
        df = df[present_cell_types]  # Keep only present cell types in order

        # Identify outliers using IQR
        filtered_data = {}
        for cell_type in present_cell_types:
            y = df[cell_type].dropna()
            Q1, Q3 = y.quantile(0.25), y.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound, upper_bound = Q1 - 1.5 * IQR, Q3 + 1.5 * IQR
            filtered_data[cell_type] = y[(y >= lower_bound) & (y <= upper_bound)]

        # Prepare boxplot data
        metric_x_positions = []
        for j, cell_type in enumerate(present_cell_types):
            x_pos = i * 1.5 + j * 0.25  # Adjust spacing so groups stay together
            metric_x_positions.append(x_pos)
            all_box_data.append(filtered_data[cell_type])
            positions.append(x_pos)
            color_list.append(color_map.get(cell_type, "gray"))

            # Scatter data points (closer to center)
            y_values = filtered_data[cell_type]
            x_values = np.random.normal(x_pos, 0.005, size=len(y_values))  # Keep jitter minimal

            scatter_x.extend(x_values)
            scatter_y.extend(y_values)
            scatter_colors.extend([color_map.get(cell_type, "gray")] * len(y_values))

        # Pick the middle boxplot position for x-axis labels
        if metric_x_positions:
            middle_pos = metric_x_positions[len(metric_x_positions) // 2]
            metric_positions.append(middle_pos)

    # Create boxplot
    box = ax.boxplot(all_box_data, patch_artist=True, positions=positions, widths=0.1, showfliers=show_outliers)

    # Set colors for boxplots
    for patch, color in zip(box["boxes"], color_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)

    for median in box["medians"]:
        median.set_color("black")
        median.set_linewidth(2.5)

    # Scatter plot overlay (points are now more centered)
    ax.scatter(scatter_x, scatter_y, color=scatter_colors, edgecolors="none", s=30, alpha=0.9)

    # Add legend for cell types
    legend_handles = [plt.Line2D([0], [0], color=color_map[cell], lw=4, label=cell) for cell in fixed_order if cell in color_map]
    ax.legend(handles=legend_handles, title="Cell Types")

    # Set x-axis labels under the middle boxplot for each metric
    ax.set_xticks(metric_positions)
    ax.set_xticklabels(metrics)

    ax.set_ylim(-3, 6)  # Y-axis limits
    ax.set_ylabel("Feature value (z-score)")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
    print(f"Boxplots saved to {save_path}")


def main(contamination_ratio=0.1, confidence_ratio_threshold=1.5, fig_output_types=[".png", ".eps"]):
    dropbox_path = get_dropbox_path()
    print(dropbox_path)
    #folder_inside_dropbox = "ExperimentOutput/Ephys4Trace1/MainFolder/"
    
    
    

    mice_groups = get_mouse_groups()
    
    mice_to_analyze = {group: [] for group in mice_groups}
    mouse_files = {group: {} for group in mice_groups}  # Use a dictionary instead of a list
    
    # Check file existence
    for group, mice in mice_groups.items():
        for mouse in mice:
            counts_file_new = f"_cell_type_counts_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv"
            counts_file_old = "_overall_discharge_stats.tsv"
            results_folder = os.path.join(dropbox_path, "AnalysisOutput/c4 results stats/")
            file_path_new = os.path.join(results_folder, mouse, f"{mouse}{counts_file_new}")
            file_path_old = os.path.join(results_folder, mouse, f"{mouse}{counts_file_old}")
            
            # Select the correct file if it exists
            if os.path.isfile(file_path_new):
                file_path = file_path_new
            elif os.path.isfile(file_path_old):
                file_path = file_path_old
            else:
                continue  # Skip if neither file exists
    
            # Determine if it's a "Switch" mouse
            is_switch = mouse in mice_groups["Switch"] or mice_groups["Naive"]
            
            # Read firstWideSession.txt if it's a switch mouse
            first_wide_session = None
            if is_switch:
                wide_session_file = os.path.join(results_folder, mouse, "firstWideSession.txt")
                if os.path.isfile(wide_session_file):
                    with open(wide_session_file, 'r') as f:
                        first_wide_session = f.read().strip()
    
            # Store mouse and file path
            mice_to_analyze[group].append(mouse)
            mouse_files[group][mouse] = {"file": file_path, "first_wide_session": first_wide_session}

    
    total_discharge_stats_per_cell_type = {group: [] for group in mice_groups}
    total_discharge_stats_per_mouse = {}
    
    for group, mice in mice_to_analyze.items():
        for mouse in mice:
            file_info = mouse_files[group][mouse]
            file_path = file_info["file"]
            first_wide_session = file_info["first_wide_session"]  # This is None for Wide, Narrow and Bimodal
            
            
            cell_type_counts = read_tsv(file_path)
            if cell_type_counts.empty:
                continue
            
            # If the mouse is a "Switch" mouse, split its sessions
            if first_wide_session:
                sessions = list(cell_type_counts.index.unique())
                index = sessions.index(first_wide_session)  # Find its position
                last_single_session = sessions[index - 1] if index > 0 else None  # Get the previous one
                switch_index = sessions.index("SwitchSessionStitching") if "SwitchSessionStitching" in sessions else None
                if switch_index:
                    last_wide_session = sessions[switch_index - 1] if switch_index and switch_index > 0 else None                    
                    wide_sessions = cell_type_counts.loc[first_wide_session:last_wide_session].iloc[:] # Wide sessions start **from** firstWideSession
                else:
                    wide_sessions = cell_type_counts.loc[first_wide_session:].iloc[:] # Wide sessions start **from** firstWideSession
                single_sessions = cell_type_counts.loc[:last_single_session].iloc[:]  # Exclude firstWideSession
                    
                
                if mouse in mice_groups["Switch"]:
                    # Store separately as "Single" and "WideNovices"
                    total_discharge_stats_per_cell_type.setdefault("Single", []).append(single_sessions)
                    total_discharge_stats_per_cell_type.setdefault("WideNovices", []).append(wide_sessions)
                    
                    total_discharge_stats_per_mouse[mouse + "_Single"] = single_sessions
                    total_discharge_stats_per_mouse[mouse + "_WideNovice"] = wide_sessions
                elif mouse in mice_groups["Naive"]:
                    # Store separately as "Single" and "WideNovices"
                    total_discharge_stats_per_cell_type.setdefault("SingleNaive", []).append(single_sessions)
                    total_discharge_stats_per_cell_type.setdefault("WideNaive", []).append(wide_sessions)
                    
                    total_discharge_stats_per_mouse[mouse + "_SingleNaive"] = single_sessions
                    total_discharge_stats_per_mouse[mouse + "_WideNaive"] = wide_sessions
            elif not first_wide_session and mouse in mice_groups["Switch"]:
                total_discharge_stats_per_cell_type.setdefault("Single", []).append(cell_type_counts.iloc[:, :])
                total_discharge_stats_per_mouse[mouse + "_Single"] = cell_type_counts.iloc[:, :]           
            elif not first_wide_session and mouse in mice_groups["Naive"]:
                total_discharge_stats_per_cell_type.setdefault("SingleNaive", []).append(cell_type_counts.iloc[:, :])
                total_discharge_stats_per_mouse[mouse + "_SingleNaive"] = cell_type_counts.iloc[:, :]           
            else:
                # Normal case: just store as before
                total_discharge_stats_per_cell_type.setdefault(group, []).append(cell_type_counts.iloc[:, :])
                total_discharge_stats_per_mouse[mouse] = cell_type_counts.iloc[:, :]

    """
    # Ensure only non-empty groups are processed
    results_df = pd.DataFrame({
        group: pd.concat(counts).groupby(level=0)
        for group, counts in total_discharge_stats_per_cell_type.items() if counts
    })

    # Replace NaN values with zero
    results_df = results_df.fillna(0)

    # Save to file in results_folder
    output_file = os.path.join(results_folder, f"all_cell_type_totals_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv")
    results_df.to_csv(output_file, sep='\t')
    print(f"Saved cell type totals to {output_file}")
    """
    #groups = total_discharge_stats_per_cell_type.keys()
    all_mice_data = defaultdict(list)
    all_mice_stats = []
    for group in total_discharge_stats_per_cell_type.keys():
        if group in ["WideExperts", "WideNovices"]:

            for m in range(len(total_discharge_stats_per_cell_type[group])):
                overall_stats = total_discharge_stats_per_cell_type[group][m]
                
                df = overall_stats
                for session, row in df.iterrows():
                    if not session=="SwitchSessionStitching":
                        cell_type = row["Cell Type"]
                        data_entry = [row["Cluster ID"], row["Mean Firing Rate"], row["Mean CV"], row["Mean CV2"], row["Median ISI"]]
                        all_mice_data[cell_type].append(data_entry)
                        all_mice_stats.append(row.tolist())  # Store for saving
               
        
    # Convert to numpy and apply Z-score
    for cell_type in all_mice_data:
        all_mice_data[cell_type] = np.array(all_mice_data[cell_type])
        all_mice_data[cell_type][:, 1:] = zscore(all_mice_data[cell_type][:, 1:], axis=0, nan_policy='omit')

    # Save the summary plot

    for ext in fig_output_types:
        save_path_plot = os.path.join(results_folder, f"Wide_Experts_and_Novices_overall_discharge_boxplots_zscore{ext}")
        save_boxplots(all_mice_data, save_path_plot, "Overall Discharge Statistics for Wide Experts and Novices together")

    print(f"Saved bar chart to {save_path_plot}")
            
            
    for group in total_discharge_stats_per_cell_type.keys():
        if total_discharge_stats_per_cell_type[group]:
            all_mice_data = defaultdict(list)
            all_mice_stats = []
            for m in range(len(total_discharge_stats_per_cell_type[group])):
                overall_stats = total_discharge_stats_per_cell_type[group][m]
                
                df = overall_stats
                for session, row in df.iterrows():
                    if not session=="SwitchSessionStitching":
                        cell_type = row["Cell Type"]
                        data_entry = [row["Cluster ID"], row["Mean Firing Rate"], row["Mean CV"], row["Mean CV2"], row["Median ISI"]]
                        all_mice_data[cell_type].append(data_entry)
                        all_mice_stats.append(row.tolist())  # Store for saving
               
        
            # Convert to numpy and apply Z-score
            for cell_type in all_mice_data:
                all_mice_data[cell_type] = np.array(all_mice_data[cell_type])
                all_mice_data[cell_type][:, 1:] = zscore(all_mice_data[cell_type][:, 1:], axis=0, nan_policy='omit')
    
        # Save the summary plot
    
            for ext in fig_output_types:
                save_path_plot = os.path.join(results_folder, f"{group}_overall_discharge_boxplots_zscore{ext}")
                save_boxplots(all_mice_data, save_path_plot, f"Overall Discharge Statistics for group: {group}")
        
            print(f"Saved bar chart to {save_path_plot}")
            
            
    

    #per_mouse_chart_path = os.path.join(results_folder, f"all_per_mouse_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
    #plot_per_mouse_group_bar_chart(total_discharge_stats_per_mouse, per_mouse_chart_path)
    #print(f"Saved per-mouse-group bar chart to {per_mouse_chart_path}")
    
if __name__ == "__main__":
    main()