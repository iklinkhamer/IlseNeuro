# -*- coding: utf-8 -*-
"""
Created on Wed May 21 15:04:17 2025

@author: Ilse Klinkhamer
"""
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from get_dropbox_path_2 import get_dropbox_path
from collections import defaultdict
from plot_utils_IK import c4_colors_rgb, lighten, normalize_RGB_dict


def load_mice_metadata(directory):
    """Load mouse metadata."""
    
    json_file_path = os.path.join(get_dropbox_path(), "ContextMouseExperiments\Ilse", "MiceMetadata_CSG.json")
    
    with open(json_file_path, 'r') as f:
        mice_metadata = json.load(f)
    
    mouse_names = [mouse for mouse in mice_metadata.keys() if mouse not in ["ReserveMouse1", "ReserveMouse4"]]
    
    # Remove ReserveMouse1 and ReserveMouse4 from metadata
    mice_metadata = {mouse: data for mouse, data in mice_metadata.items() if mouse not in ["ReserveMouse1", "ReserveMouse4"]}
    
    return mouse_names, mice_metadata

def get_mouse_groups():
    """Returns a dictionary containing categorized mouse groups."""
    return {
        "CSG": [
            "Geneva", "Georgetown_Cbx", "Georgetown_Ctx", "Helsinki", "Ilo_Cbx", "Ilo_mPFC", "Kourou_Cbx"
            ]
    }

def read_tsv(file_path):
    """Reads a TSV file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', index_col=0)

def plot_pie_charts(total_counts_per_cell_type, save_path):
    """Plot and save pie charts for all mouse groups in one figure."""
    num_groups = total_counts_per_cell_type.shape[1]
    fig, axes = plt.subplots(1, num_groups, figsize=(num_groups * 5, 5))
    
    if num_groups == 1:
        axes = [axes]
    """
    color_map = {
        'PkC_cs': 'grey',
        'MLI': 'pink',
        'MFB': 'red',
        'GoC': 'green',
        'PkC_ss': 'blue'
    }
    """
    color_map = c4_colors_rgb() 
    color_map['PkC_cs'] = lighten(color_map['PkC_cs'])
    color_map = normalize_RGB_dict(color_map)
    
    
    for ax, (group, counts) in zip(axes, total_counts_per_cell_type.items()):
        labels = [f"{cell_type} ({count})" for cell_type, count in counts.items()]
        colors = [color_map.get(cell_type, 'gray') for cell_type in counts.index]
        ax.pie(counts, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)
        ax.set_title(group)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_grouped_bar_chart(total_counts_per_cell_type, save_path):
    """Plot and save a grouped bar chart for all mouse groups."""
    df = pd.DataFrame(np.transpose(total_counts_per_cell_type))
    color_map = c4_colors_rgb() 
    color_map['PkC_cs'] = lighten(color_map['PkC_cs'])
    color_map = normalize_RGB_dict(color_map)
    
    
    df_test = pd.DataFrame(total_counts_per_cell_type)
    
    df.plot(kind='bar', figsize=(12, 6), color=[color_map.get(cell, 'gray') for cell in df_test.index])
    plt.xlabel("Mouse Group")
    plt.ylabel("Neuron Count")
    plt.title("Neuron Counts per Cell Type Across Mouse Groups")
    plt.xticks(rotation=45)
    plt.legend(title="Cell Type")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_per_mouse_group_bar_chart(total_counts_per_mouse, save_path):
    """Plot and save a bar chart with mouse groups on x-axis and bars for different cell types per mouse group."""
    df = pd.DataFrame(total_counts_per_mouse).T.fillna(0).astype(int)
    color_map = c4_colors_rgb() 
    color_map['PkC_cs'] = lighten(color_map['PkC_cs'])
    color_map = normalize_RGB_dict(color_map)

    num_mice = len(df.index)
    fig_width = max(12, num_mice * 0.5)  # Scale width dynamically
    fig_height = 8  # Slightly taller for clarity

    df.plot(kind='bar', figsize=(fig_width, fig_height), width=0.7, color=[color_map.get(cell, 'gray') for cell in df.columns])

    plt.xlabel("Mouse Name")
    plt.ylabel("Neuron Count")
    plt.title("Neuron Counts per Cell Type per Mouse Group")
    
    plt.xticks(ticks=range(num_mice), labels=df.index, rotation=60, ha='right', fontsize=10)  # More rotation
    plt.yticks(fontsize=10)
    
    plt.legend(title="Cell Type", fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()  # Prevent labels from getting cut off
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def main(contamination_ratio=0.1, confidence_ratio_threshold=1.5, in_analyzed_ephys_and_c4=False):
    dropbox_path = get_dropbox_path()
    print(dropbox_path)
    #folder_inside_dropbox = "ExperimentOutput/Ephys4Trace1/MainFolder/"
    
    counts_folder = os.path.join(dropbox_path, f"ContextMouseExperiments/Ilse/AnalysisOutput/Cell_type_counts/fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}/")
    #os.makedirs(counts_folder, exist_ok=True)

    mice_groups = get_mouse_groups()
    
    mice_to_analyze = {group: [] for group in mice_groups}
    mouse_files = {group: {} for group in mice_groups}  # Use a dictionary instead of a list
    
    mouse_names, mice_metadata = load_mice_metadata(get_dropbox_path())
    
    if in_analyzed_ephys_and_c4:    

        # Read the TSV files into DataFrames
        neurons_c4_prior = pd.read_csv(os.path.join(get_dropbox_path(), 'AnalysisOutput/NeuronCounts/c4/neuron_prior_labels.tsv'), sep='\t')
        neurons_c4 = pd.read_csv(os.path.join(get_dropbox_path(), 'AnalysisOutput/NeuronCounts/c4/neurons.tsv'), sep='\t')
        neurons_analyzed_ephys = pd.read_csv(os.path.join(get_dropbox_path(), 'AnalysisOutput/NeuronCounts/analyzedEphys/neurons.tsv'), sep='\t')
        
        # Perform an inner join on the 'NeuronId' column
        mouse_df = pd.merge(neurons_analyzed_ephys, neurons_c4_prior, how='inner', on='NeuronId')
        #mouse_df['cluster_id']=""
        #mouse_df['predicted_cell_type']=""
        
        
        
        
        for mouse in mouse_names:
            folder_names = mice_metadata[mouse]["folders"]
            for folder in folder_names:
                if os.path.exists(os.path.join(get_dropbox_path(), f"AnalysisOutput/c4 results stats/{mouse}/{folder}/cluster_predicted_cell_type.tsv")):
                    c4_results = pd.read_csv(os.path.join(get_dropbox_path(), f"AnalysisOutput/c4 results stats/{mouse}/{folder}/cluster_predicted_cell_type.tsv"), sep='\t')
                    c4_results['NeuronId'] = f"TRACE_{folder}_" + c4_results['cluster_id'].astype(str)
                    #merged_df = pd.merge(mouse_df, c4_results, how='outer', on='NeuronId')
                    if not 'mouse_c4_df' in locals():
                        mouse_c4_df = c4_results 
                    else:
                        mouse_c4_df = pd.concat([mouse_c4_df, c4_results], ignore_index=True)

        mouse_df = pd.merge(mouse_df, mouse_c4_df, how='inner', on='NeuronId')
        
        priors = np.unique(mouse_df["Prior"])
        total_counts_per_cell_type = {prior: {} for prior in priors}
        cell_types = np.unique(mouse_df["predicted_cell_type"])
        cell_counts_group = defaultdict()
        
        # Get wide novices once outside the loop
        wide_novices = [mouse for mouse, metadata in mice_metadata.items() if metadata.get("distType") == "DEUN"]
        
        for prior in priors:
            for cell_type in cell_types:
                # General count for this prior and cell type, excluding "Ana" mice
                count = ((mouse_df["Prior"] == prior) & 
                         (mouse_df["predicted_cell_type"] == cell_type) & 
                         (~mouse_df["Mouse"].str.contains("Ana", case=False, na=False))).sum()
        
                if prior == "Wide":  # Special handling for "Wide" prior
                    count_experts = ((mouse_df["Prior"] == prior) & 
                                     (mouse_df["predicted_cell_type"] == cell_type) & 
                                     (~mouse_df["Mouse"].str.contains("Ana", case=False, na=False)) & 
                                     (~mouse_df["Mouse"].isin(wide_novices))).sum()
        
                    count_novices = ((mouse_df["Prior"] == prior) & 
                                     (mouse_df["predicted_cell_type"] == cell_type) & 
                                     (~mouse_df["Mouse"].str.contains("Ana", case=False, na=False)) & 
                                     (mouse_df["Mouse"].isin(wide_novices))).sum()
        
                    counts = [count, count_experts, count_novices]
                    prior_labels = ["Wide", "WideExperts", "WideNovices"]
                else:
                    counts = [count]
                    prior_labels = [prior]
        
                for p, prior_label in enumerate(prior_labels):
                    # Check if the prior label exists in the dictionary and ensure it's initialized
                    if prior_label not in total_counts_per_cell_type:
                        total_counts_per_cell_type[prior_label] = {}
        
                    # Initialize the cell_type if it doesn't exist
                    if cell_type not in total_counts_per_cell_type[prior_label]:
                        total_counts_per_cell_type[prior_label][cell_type] = 0
        
                    # Update the counts correctly
                    total_counts_per_cell_type[prior_label][cell_type] += counts[p]
        
        # Debugging: check the final result to see if anything was overwritten incorrectly
        print(total_counts_per_cell_type)

        # Convert total_counts_per_cell_type to a DataFrame
        results_df = pd.DataFrame(total_counts_per_cell_type)  # Transpose to get cell_types as index and priors as columns        
   
        # Display the resulting DataFrame
        print(results_df)
        
        # Replace NaN values with zero
        results_df = results_df.fillna(0).astype(int)
        
        # Generate and save plots
        
        GoCs= 0
        MFBs=0
        PkC_sss=0
        PkC_css=0
        MLIs=0
     
        for group, counts in results_df.items():
            labels = [f"{cell_type} ({count})" for cell_type, count in counts.items()]
            if not group in ["SingleNaive", "WideNaive", "Wide"]:
                GoCs+=counts.GoC
                MFBs+=counts.MFB
                PkC_css+=counts.PkC_cs
                PkC_sss+=counts.PkC_ss
                MLIs+=counts.MLI            
            else:
                a = 1
        results_all_mice = {
            'Index': ['GoC', 'MFB', 'MLI', 'PkC_cs', 'PkC_ss'],
            'AllMice': [GoCs, MFBs, MLIs, PkC_css, PkC_sss]
            }
        results_all_mice_df = pd.DataFrame(results_all_mice)
        results_all_mice_df = results_all_mice_df.set_index('Index')
        a =1
        
        """
        pie_chart_path = os.path.join(counts_folder, f"aaa_new_overall_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
            
        plot_pie_charts(results_all_mice_df, pie_chart_path)
            
        pie_chart_path = os.path.join(counts_folder, f"aaa_new_overall_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
          
        plot_pie_charts(results_all_mice_df, pie_chart_path)
        
        # Save to file in counts_folder
        output_file = os.path.join(counts_folder, f"aaa_new_overall_all_cell_type_totals_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv")
        results_all_mice_df.to_csv(output_file, sep='\t')
        print(f"Saved cell type totals to {output_file}")
        

        # Save to file in counts_folder
        output_file = os.path.join(counts_folder, f"all_new_cell_type_totals_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv")
        results_df.to_csv(output_file, sep='\t')
        print(f"Saved cell type totals to {output_file}")
        
        # Generate and save plots
        pie_chart_path = os.path.join(counts_folder, f"all_new_groups_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
        bar_chart_path = os.path.join(counts_folder, f"all_new_groups_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
        
        plot_pie_charts(results_df, pie_chart_path)
        plot_grouped_bar_chart(results_df, bar_chart_path)
        
        pie_chart_path = os.path.join(counts_folder, f"all_new_groups_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
        bar_chart_path = os.path.join(counts_folder, f"all_new_groups_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
        
        plot_pie_charts(results_df, pie_chart_path)
        plot_grouped_bar_chart(results_df, bar_chart_path)
        
        print(f"Saved pie chart to {pie_chart_path}")
        print(f"Saved bar chart to {bar_chart_path}")
        """

    else:
            
        # Check file existence
        for group, mice in mice_groups.items():
            for mouse in mice:
                counts_file_new = f"_cell_type_counts_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv"
                counts_file_old = "_cell_type_counts.tsv"
                
                file_path_new = os.path.join(counts_folder, f"{mouse}{counts_file_new}")
                file_path_old = os.path.join(counts_folder, f"{mouse}{counts_file_old}")
                
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
                    wide_session_file = os.path.join(dropbox_path, f"AnalysisOutput/c4 results stats/{mouse}/firstWideSession.txt")
                    if os.path.isfile(wide_session_file):
                        with open(wide_session_file, 'r') as f:
                            first_wide_session = f.read().strip()
        
                # Store mouse and file path
                mice_to_analyze[group].append(mouse)
                mouse_files[group][mouse] = {"file": file_path, "first_wide_session": first_wide_session}
    
        
        total_counts_per_cell_type = {group: [] for group in mice_groups}
        total_counts_per_mouse = {}
        
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
                    single_sessions = cell_type_counts.loc[:first_wide_session].iloc[:-1]  # Exclude firstWideSession
                    wide_sessions = cell_type_counts.loc[first_wide_session:].iloc[:-2] # Wide sessions start **from** firstWideSession
                    
                    if mouse in mice_groups["Switch"]:
                        # Store separately as "Single" and "WideNovices"
                        total_counts_per_cell_type.setdefault("Single", []).append(single_sessions.sum().astype(int))
                        total_counts_per_cell_type.setdefault("WideNovices", []).append(wide_sessions.sum().astype(int))
                        
                        total_counts_per_mouse[mouse + "_Single"] = single_sessions.sum().astype(int)
                        total_counts_per_mouse[mouse + "_WideNovice"] = wide_sessions.sum().astype(int)
                    elif mouse in mice_groups["Naive"]:
                        # Store separately as "Single" and "WideNovices"
                        total_counts_per_cell_type.setdefault("SingleNaive", []).append(single_sessions.sum().astype(int))
                        total_counts_per_cell_type.setdefault("WideNaive", []).append(wide_sessions.sum().astype(int))
                        
                        total_counts_per_mouse[mouse + "_SingleNaive"] = single_sessions.sum().astype(int)
                        total_counts_per_mouse[mouse + "_WideNaive"] = wide_sessions.sum().astype(int)
                elif not first_wide_session and mouse in mice_groups["Switch"]:
                    total_counts_per_cell_type.setdefault("Single", []).append(cell_type_counts.iloc[-2, :].astype(int))
                    total_counts_per_mouse[mouse + "_Single"] = cell_type_counts.iloc[-2, :].astype(int)           
                elif not first_wide_session and mouse in mice_groups["Naive"]:
                    total_counts_per_cell_type.setdefault("SingleNaive", []).append(cell_type_counts.iloc[-2, :].astype(int))
                    total_counts_per_mouse[mouse + "_SingleNaive"] = cell_type_counts.iloc[-2, :].astype(int)           
                else:
                    # Normal case: just store as before
                    total_counts_per_cell_type.setdefault(group, []).append(cell_type_counts.iloc[-2, :].astype(int))
                    total_counts_per_mouse[mouse] = cell_type_counts.iloc[-2, :].astype(int)
    
        
        # Ensure only non-empty groups are processed
        results_df = pd.DataFrame({
            group: pd.concat(counts).groupby(level=0).sum().astype(int)
            for group, counts in total_counts_per_cell_type.items() if counts
        })
    
        # Replace NaN values with zero
        results_df = results_df.fillna(0).astype(int)
        
        # Generate and save plots
       
        GoCs= 0
        MFBs=0
        PkC_sss=0
        PkC_css=0
        MLIs=0
     
        for group, counts in results_df.items():
            labels = [f"{cell_type} ({count})" for cell_type, count in counts.items()]
            if not group in ["SingleNaive", "WideNaive"]:
                GoCs+=counts.GoC
                MFBs+=counts.MFB
                PkC_css+=counts.PkC_cs
                PkC_sss+=counts.PkC_ss
                MLIs+=counts.MLI            
            else:
                a = 1
        results_all_mice = {
            'Index': ['GoC', 'MFB', 'MLI', 'PkC_cs', 'PkC_ss'],
            'AllMice': [GoCs, MFBs, MLIs, PkC_css, PkC_sss]
            }
        results_all_mice_df = pd.DataFrame(results_all_mice)
        results_all_mice_df = results_all_mice_df.set_index('Index')
        a = 1
        
        """
        
        pie_chart_path = os.path.join(counts_folder, f"aaa_overall_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
            
        plot_pie_charts(results_all_mice_df, pie_chart_path)
            
        pie_chart_path = os.path.join(counts_folder, f"aaa_overall_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
          
        plot_pie_charts(results_all_mice_df, pie_chart_path)
        
        # Save to file in counts_folder
        output_file = os.path.join(counts_folder, f"aaa_overall_all_cell_type_totals_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv")
        results_all_mice_df.to_csv(output_file, sep='\t')
        print(f"Saved cell type totals to {output_file}")
        
    
        # Save to file in counts_folder
        output_file = os.path.join(counts_folder, f"all_cell_type_totals_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv")
        results_df.to_csv(output_file, sep='\t')
        print(f"Saved cell type totals to {output_file}")
        
        # Generate and save plots
        pie_chart_path = os.path.join(counts_folder, f"all_groups_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
        bar_chart_path = os.path.join(counts_folder, f"all_groups_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
        
        plot_pie_charts(results_df, pie_chart_path)
        plot_grouped_bar_chart(results_df, bar_chart_path)
        
        pie_chart_path = os.path.join(counts_folder, f"all_groups_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
        bar_chart_path = os.path.join(counts_folder, f"all_groups_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
        
        plot_pie_charts(results_df, pie_chart_path)
        plot_grouped_bar_chart(results_df, bar_chart_path)
        
        print(f"Saved pie chart to {pie_chart_path}")
        print(f"Saved bar chart to {bar_chart_path}")
    
        per_mouse_chart_path = os.path.join(counts_folder, f"all_per_mouse_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
        plot_per_mouse_group_bar_chart(total_counts_per_mouse, per_mouse_chart_path)
        print(f"Saved per-mouse-group bar chart to {per_mouse_chart_path}")
        per_mouse_chart_path = os.path.join(counts_folder, f"all_per_mouse_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
        plot_per_mouse_group_bar_chart(total_counts_per_mouse, per_mouse_chart_path)
        
        """
    
    
if __name__ == "__main__":
    main()

