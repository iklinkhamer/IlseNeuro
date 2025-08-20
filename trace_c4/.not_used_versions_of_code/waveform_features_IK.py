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
from npyx.feat import waveform_features_json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from get_dropbox_path import get_dropbox_path
from plot_utils_IK import c4_colors_rgb, lighten, normalize_RGB_dict
from mouseUtils import getMouseFolders
from scipy.stats import zscore

def save_boxplots(all_sessions_data, save_path, title, show_outliers=False):
    """Generate and save a single boxplot figure with all metrics and individual data points.

    Arguments:
        - all_sessions_data: dict, structured spike metrics
        - save_path: str, where to save the plot
        - title: str, plot title
        - show_outliers: bool, if True, boxplot shows outliers, but scatter never does
    """
    metrics = ["Peak half-width", "Peak to trough ratio", "Waveform width"]
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



def get_discharge_statistics(mouse_name = "Iowa"
                             , switch_sessions=False
                             , contamination_ratio=0.1
                             , confidence_ratio_threshold=2
                             , directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")            
                             , save_directory=os.path.join(get_dropbox_path(), "AnalysisOutput/c4 results stats/")
                             , dat_dir=os.path.join(get_dropbox_path(), "C4_conversion")
                             , again=False
                             , fig_output_types=[".png", ".eps"]
                             ):
    dp_base = os.path.join(directory, mouse_name)
    if "ReserveMouse" in mouse_name:
        dp_base = dp_base.replace("MainFolder", "ReserveFolder")
    save_base = os.path.join(save_directory,mouse_name)
    """
    mouse_folders = [
        folder for folder in os.listdir(dp_base)
        if os.path.isdir(os.path.join(dp_base, folder)) and mouse_name in folder and "copy" not in folder.lower()
    ]
    mouse_folders.sort()
    """
    mouse_folders = getMouseFolders(mouse_name)

    if switch_sessions:
        switch_folder = os.path.join(dp_base, "SwitchSessionStitching")
        if os.path.exists(switch_folder):
            switch_folder_name = "SwitchSessionStitching"                
            mouse_folders.append(switch_folder_name)             

    phy_folder = "c4"
    all_sessions_data_plot = defaultdict(list)  # Store data by cell type across all sessions
    all_sessions_data_stats = []
    update_final_fig=True
    i = 1
    for sess in mouse_folders:
      #if i == 1:
       #    i = 2
        #   continue
        
        if not again and os.path.isfile(os.path.join(save_directory, mouse_name, sess, f"{sess}_wvf_feats.tsv")) and os.path.isfile(os.path.join(save_directory, mouse_name, sess, f"{sess}_wvf_feats_boxplots.png")) and os.path.isfile(os.path.join(save_directory, mouse_name, f"{mouse_name}_overall_wvf_feats_boxplots.png")) and os.path.isfile(os.path.join(save_directory, mouse_name, f"{mouse_name}_overall_wvf_feats.tsv")):
            print(f"Session {sess} already done, skipping this session...")
            update_final_fig=False
            continue
        
        print(f"Processing session: {sess}")
        
        save_path = os.path.join(save_base, sess)
        os.makedirs(save_path, exist_ok=True)
        
        dat_path = os.path.join(dat_dir, mouse_name, sess)
           
               
        dp = os.path.join(dp_base, sess, phy_folder)
        specific_c4_results_folder = f"c4_results_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}"

        unit_file = os.path.join(dp, specific_c4_results_folder, "cluster_predicted_cell_type.tsv")
        if not os.path.exists(unit_file):
            continue  # Skip if the file doesn't exist      
                
            
        df_units = pd.read_csv(unit_file, sep="\t", usecols=["cluster_id", "predicted_cell_type"])
        df_units = df_units.dropna()

        session_data = defaultdict(list)  # Store per-session data by cell type
        session_data_stats = []  # Store per-session data
        for _, row in df_units.iterrows():
            cluster_id, cell_type = int(row["cluster_id"]), row["predicted_cell_type"]
            [wvf_feats,waveform_features_IK] = waveform_features_json(dp, cluster_id, plot_debug=False, dat_dir=dat_path)

            data_entry = [cluster_id, waveform_features_IK["pos_half_width"], waveform_features_IK["peak_to_trough_ratio"], waveform_features_IK["wvf_width"]]            
            session_data[cell_type].append(data_entry)
            if not sess == "SwitchSessionStitching":
                all_sessions_data_plot[cell_type].append(data_entry)
            data_entry = [sess, cluster_id, cell_type, waveform_features_IK["pos_half_width"], waveform_features_IK["peak_to_trough_ratio"], waveform_features_IK["wvf_width"]]
            session_data_stats.append(data_entry)
            all_sessions_data_stats.append(data_entry)

        # Save per-session TSV
        df_session = pd.DataFrame(session_data_stats, columns=["Session", "Cluster ID", "Cell Type", "Peak half-width", "Peak to trough ratio", "Waveform width"])
        session_tsv = os.path.join(save_path, f"{sess}_wvf_feats.tsv")
        df_session.to_csv(session_tsv, sep="\t", index=False)
        print(f"Session wvf_feats saved to {session_tsv}")
        
        # Z-score each metric for each cell type
        for cell_type in session_data:
            session_data[cell_type] = np.array(session_data[cell_type])
            session_data[cell_type][:, 1:] = zscore(session_data[cell_type][:, 1:], axis=0, nan_policy='omit')
            
            
        if any(len(v) > 0 for v in session_data.values()):
            # Save session boxplots
            for ext in fig_output_types:
                session_plot = os.path.join(save_path, f"{sess}_wvf_feats_boxplots_zscore{ext}")
                save_boxplots(session_data, session_plot, f"Waveform Features - {sess}")
        else:
            print(f"No neurons found in session {sess}, not making boxplot figure...")
    if update_final_fig:
        # Save overall statistics
        df_overall = pd.DataFrame(all_sessions_data_stats, columns=["Session", "Cluster ID", "Cell Type", "Peak half-width", "Peak to trough ratio", "Waveform width"])
        overall_tsv = os.path.join(save_base, f"{mouse_name}_overall_wvf_feats.tsv")
        df_overall.to_csv(overall_tsv, sep="\t", index=False)
        print(f"Overall wvf_feats saved to {overall_tsv}")
               

        # Z-score each metric for each cell type
        for cell_type in all_sessions_data_plot:
            all_sessions_data_plot[cell_type] = np.array(all_sessions_data_plot[cell_type])
            all_sessions_data_plot[cell_type][:, 1:] = zscore(all_sessions_data_plot[cell_type][:, 1:], axis=0, nan_policy='omit')


        if all_sessions_data_plot:
            # Save combined boxplots
            for ext in fig_output_types:
                overall_plot = os.path.join(save_base, f"{mouse_name}_overall_wvf_feats_boxplots_zscore{ext}")
                save_boxplots(all_sessions_data_plot, overall_plot, "Overall Waveform Features")
        
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

def summarize_discharge_statistics(mice, save_directory, fig_output_types=[".png", ".eps"]):
    """Aggregate and plot discharge statistics for all mice together and save as a TSV."""
    
    all_mice_data = defaultdict(list)
    all_mice_stats = []

    # Exclude "Ana" mice
    filtered_mice = [mouse for mouse in mice if not mouse.startswith("Ana")]

    for mouse in filtered_mice:
        overall_stats_path = os.path.join(save_directory, mouse, f"{mouse}_overall_discharge_stats.tsv")
        if os.path.isfile(overall_stats_path):
            df = pd.read_csv(overall_stats_path, sep="\t")
            for _, row in df.iterrows():
                if not row["Session"]=="SwitchSessionStitching":
                    cell_type = row["Cell Type"]
                    data_entry = [row["Cluster ID"], row["Mean Firing Rate"], row["Mean CV"], row["Mean CV2"], row["Median ISI"]]
                    all_mice_data[cell_type].append(data_entry)
                    all_mice_stats.append(row.tolist())  # Store for saving
        else:
            print(f"Skipping {mouse}, no overall discharge file found.")

    # Convert to numpy and apply Z-score
    for cell_type in all_mice_data:
        all_mice_data[cell_type] = np.array(all_mice_data[cell_type])
        all_mice_data[cell_type][:, 1:] = zscore(all_mice_data[cell_type][:, 1:], axis=0, nan_policy='omit')

    # Save the summary plot
    for ext in fig_output_types:
        save_path_plot = os.path.join(save_directory, f"AllMice_overall_discharge_boxplots_zscore{ext}")
        save_boxplots(all_mice_data, save_path_plot, "Overall Discharge Statistics for All Mice (Excluding Ana Mice)")

    # Save summary statistics as a TSV file
    if all_mice_stats:
        df_summary = pd.DataFrame(all_mice_stats, columns=["Session", "Cluster ID", "Cell Type", "Mean Firing Rate", "Mean CV", "Mean CV2", "Median ISI"])
        save_path_tsv = os.path.join(save_directory, "AllMice_overall_discharge_stats.tsv")
        df_summary.to_csv(save_path_tsv, sep="\t", index=False)
        print(f"Summary statistics saved to {save_path_tsv}")

    
def main(mouse_name="ReserveMouse3"
         , switch_sessions=True
         , contamination_ratio=0.1
         , confidence_ratio_threshold=1.5
         , directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
         , save_directory=os.path.join(get_dropbox_path(), "AnalysisOutput/c4 results stats/")
         , again=True
         ):
    
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
            mice = [mouse_name]
        else:
            try:
                mice = get_mice()
            except:
                print("Error: No mouse name provided")
                return
    else:
        mice = [mouse_name]
        
    #summarize_discharge_statistics(mice, save_directory)
    
    for mouse_name in mice:

        if "ReserveMouse" in mouse_name:
            dp_base = directory.replace("MainFolder", "ReserveFolder")
        else:
            dp_base = directory
        discharge_statistics = []
        if os.path.exists(os.path.join(dp_base,mouse_name)):
  
            discharge_statistics = get_discharge_statistics(mouse_name, switch_sessions, contamination_ratio, confidence_ratio_threshold, again=again)
            
        else:
            print(f"Data for {mouse_name} not found. Probably not synced to this computer. Skipping...")
            continue
        
    

                #sys.exit(1)
    #discharge_statistics = get_discharge_statistics(mouse_name, switch_sessions, contamination_ratio, confidence_ratio_threshold)
    return discharge_statistics


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
            ]

if __name__ == "__main__":
    main()