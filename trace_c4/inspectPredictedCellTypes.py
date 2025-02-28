#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:19:06 2024

@author: Ilse Klinkhamer
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from get_dropbox_path import get_dropbox_path
# Load the .tsv file
def load_tsv(file_path):
    """Load the TSV file and return it as a DataFrame."""
    return pd.read_csv(file_path, sep='\t')

# Compute cell type counts and fractions
def compute_fractions(data, cell_type_column):
    """Compute counts and fractions of cell types in the given column."""
    cell_type_counts = data[cell_type_column].value_counts()
    total = cell_type_counts.sum()
    fractions = cell_type_counts / total
    return cell_type_counts, fractions

# Plot and save a pie chart
def plot_pie_chart(cell_type_counts, fractions, save_path, title="Cell Type Distribution"):
    """Plot and save a pie chart of the cell type fractions."""
    labels = [f"{cell_type} ({count})" for cell_type, count in cell_type_counts.items()]
    plt.figure(figsize=(8, 8))

    # Define fixed colors for cell types
    color_map = {
        'PkC_cs': 'grey',
        'MLI': 'pink',
        'MFB': 'red',
        'GoC': 'green',
        'PkC_ss': 'blue'
    }
    colors = [color_map.get(cell_type, 'gray') for cell_type in cell_type_counts.index]

    plt.pie(
        fractions, labels=labels, colors=colors,
        autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100. * sum(cell_type_counts)))})", startangle=140
    )

    plt.title(title)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

# Plot and save a bar chart
def plot_bar_chart(cell_type_counts, save_path, title="Neuron Counts per Cell Type"):
    """Plot and save a bar chart showing neuron counts per cell type."""
    plt.figure(figsize=(8, 6))

    # Define fixed colors for cell types
    color_map = {
        'PkC_cs': 'grey',
        'MLI': 'pink',
        'MFB': 'red',
        'GoC': 'green',
        'PkC_ss': 'blue'
    }
    colors = [color_map.get(cell_type, 'gray') for cell_type in cell_type_counts.index]

    plt.bar(cell_type_counts.index, cell_type_counts.values, color=colors)
    plt.xlabel("Cell Type")
    plt.ylabel("Neuron Count")
    plt.title(title)
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

# Plot total neuron count per session
def plot_total_neurons_per_session(session_counts, save_path):
    """Plot and save a bar chart of total neuron count per session."""
    plt.figure(figsize=(10, 6))
    
    plt.bar(session_counts.keys(), session_counts.values(), color='skyblue')
    plt.xlabel("Session")
    plt.ylabel("Total Neurons")
    plt.title("Total Neurons Recorded per Session")
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
# New: Plot cell type counts across sessions
def plot_cell_type_counts_per_session(session_data, save_path):
    """Plot and save a grouped bar chart of cell type counts per session."""
    cell_types = ['PkC_cs', 'MLI', 'MFB', 'GoC', 'PkC_ss']
    sessions = list(session_data.keys())

    # Extract cell type counts per session
    counts = {cell_type: [session_data[session].get(cell_type, 0) for session in sessions] for cell_type in cell_types}

    x = np.arange(len(sessions))  # Session indices
    width = 0.15  # Width of each bar
    fig, ax = plt.subplots(figsize=(12, 6))

    color_map = {'PkC_cs': 'grey', 'MLI': 'pink', 'MFB': 'red', 'GoC': 'green', 'PkC_ss': 'blue'}

    # Plot bars for each cell type
    for i, cell_type in enumerate(cell_types):
        ax.bar(x + i * width, counts[cell_type], width, label=cell_type, color=color_map[cell_type])

    ax.set_xlabel("Session")
    ax.set_ylabel("Neuron Count")
    ax.set_title("Neuron Counts per Cell Type Across Sessions")
    ax.set_xticks(x + width)
    ax.set_xticklabels(sessions, rotation=45)
    ax.legend()
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
def save_counts_to_tsv(cell_type_counts_per_session, total_cell_type_counts, total_fractions, save_path):
    """Save cell type counts and fractions to a TSV file."""
    with open(save_path, 'w') as f:
        # Write header
        f.write("Session\t" + "\t".join(total_cell_type_counts.index) + "\n")
        
        # Write per-session counts
        for session, counts in cell_type_counts_per_session.items():
            row_counts = [str(counts.get(cell_type, 0)) for cell_type in total_cell_type_counts.index]
            f.write(session + "\t" + "\t".join(row_counts) + "\n")
        
        # Write total counts
        f.write("Total Counts\t" + "\t".join(map(str, total_cell_type_counts.values)) + "\n")
        
        # Write total fractions
        f.write("Total Fractions\t" + "\t".join(map(lambda x: f"{x:.4f}", total_fractions.values)) + "\n")



# Main function
def main(mouse_name=None,
         general_results=False,
         contamination_ratio=0.1,
         confidence_ratio_threshold=1.5,
         directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/"),
         save_dir=os.path.join(get_dropbox_path(), "AnalysisOutput/Cell_type_counts/"),
         switch_sessions=False):
    
    
    
    if not general_results:
        save_dir = os.path.join(save_dir,f"fnfpThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}")
        os.makedirs(save_dir, exist_ok=True)
    
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
        else:
            print("Error: No mouse name provided")
            sys.exit(1)
    
    dp_base = os.path.join(directory, mouse_name)
    if "ReserveMouse" in mouse_name:
        dp_base = dp_base.replace("MainFolder", "ReserveFolder")

    if switch_sessions:
        mouse_folders = [os.path.join(dp_base, "SwitchSessionStitching")]
    else:
        mouse_folders = [
            folder for folder in os.listdir(dp_base)
            if os.path.isdir(os.path.join(dp_base, folder)) 
            and mouse_name in folder 
            and "copy" not in folder
            and "Copy" not in folder
        ]
        mouse_folders.sort()   

    # Create directories for saving figures
    os.makedirs(save_dir, exist_ok=True)
    save_dir_pie_charts = os.path.join(save_dir, "pie_charts")
    os.makedirs(save_dir_pie_charts, exist_ok=True)
    save_dir_bar_charts = os.path.join(save_dir, "bar_charts")
    os.makedirs(save_dir_bar_charts, exist_ok=True)

    # Aggregate data across all sessions
    total_cell_type_counts = pd.Series(dtype=int)
    total_neurons_per_session = {}
    cell_type_counts_per_session = {}

    try:
        for folder in mouse_folders:
            if general_results:
                folder_path = os.path.join(dp_base, folder, "c4")
            else:
                folder_path = os.path.join(dp_base, folder, "c4", "c4_results_fpfnThreshold_0.1_confidenceRatio_1.5")
                print(f"Results folder: {folder_path}")
            
            if not os.path.exists(folder_path):
                print("Results folder not found, skipping.")
                continue
            
            tsv_files = [f for f in os.listdir(folder_path) if f == 'cluster_predicted_cell_type.tsv']
            confidence_ratio_files = [f for f in os.listdir(folder_path) if f == 'cluster_confidence_ratio.tsv']
            
            for tsv_file, confidence_ratio_file in zip(tsv_files, confidence_ratio_files):

                file_path = os.path.join(folder_path, tsv_file)
                print(f"Processing file: {file_path}")

                # Load data
                data = load_tsv(file_path)
                confidence_ratios = load_tsv(os.path.join(folder_path, confidence_ratio_file))

                # Ensure necessary columns exist
                required_columns = ['cluster_id', 'predicted_cell_type']
                if not all(col in data.columns for col in required_columns):
                    print(f"Skipping file {file_path}: Missing required columns.")
                    continue
                
                # Filter data where confidence ratio for each cluster_id is above the threshold
                if 'cluster_id' in confidence_ratios.columns and 'confidence_ratio' in confidence_ratios.columns:
                    valid_clusters = confidence_ratios[confidence_ratios['confidence_ratio'] > confidence_ratio_threshold]['cluster_id']
                    data = data[data['cluster_id'].isin(valid_clusters)]

                # Compute counts and fractions
                cell_type_counts, fractions = compute_fractions(data, 'predicted_cell_type')

                # Filter for specific cell types
                relevant_cell_types = ['PkC_cs', 'MLI', 'MFB', 'GoC', 'PkC_ss']
                filtered_counts = cell_type_counts[cell_type_counts.index.isin(relevant_cell_types)]
                filtered_fractions = fractions[fractions.index.isin(relevant_cell_types)]
                
                # Store session data
                cell_type_counts_per_session[folder] = filtered_counts.to_dict()
                total_neurons_per_session[folder] = cell_type_counts.sum()
                


                # Update total counts
                total_cell_type_counts = total_cell_type_counts.add(filtered_counts, fill_value=0).astype(int)
                

                # Save paths
                save_path_pie_charts = os.path.join(save_dir_pie_charts, f"{folder}_cell_type_distribution_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
                save_path_bar_charts = os.path.join(save_dir_bar_charts, f"{folder}_cell_type_counts_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")

                # Generate and save plots
                plot_pie_chart(filtered_counts, filtered_fractions, save_path_pie_charts, title=f"Cell Type Distribution for {folder}")
                plot_bar_chart(filtered_counts, save_path_bar_charts, title=f"Neuron Counts per Cell Type in {folder}")

        # Generate overall charts for all sessions combined
        if not total_cell_type_counts.empty:
            total_fractions = total_cell_type_counts / total_cell_type_counts.sum()

            overall_pie_chart_path = os.path.join(save_dir_pie_charts, f"{mouse_name}_overall_cell_type_distribution_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
            overall_bar_chart_path = os.path.join(save_dir_bar_charts, f"{mouse_name}_overall_cell_type_counts_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
            total_neurons_chart_path = os.path.join(save_dir_bar_charts, f"{mouse_name}_total_neurons_per_session_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")

            plot_pie_chart(total_cell_type_counts, total_fractions, overall_pie_chart_path, title=f"Overall Cell Type Distribution for {mouse_name}")
            plot_bar_chart(total_cell_type_counts, overall_bar_chart_path, title=f"Overall Neuron Counts per Cell Type in {mouse_name}")
            plot_total_neurons_per_session(total_neurons_per_session, total_neurons_chart_path)
            
            # Generate overall plots            
            plot_cell_type_counts_per_session(cell_type_counts_per_session, os.path.join(save_dir_bar_charts, f"{mouse_name}_cell_types_per_session_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png"))


            # Inside main function, after generating overall charts
            save_counts_to_tsv(
                cell_type_counts_per_session, 
                total_cell_type_counts, 
                total_fractions, 
                os.path.join(save_dir, f"{mouse_name}_cell_type_counts_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.tsv")
            )
        
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()



