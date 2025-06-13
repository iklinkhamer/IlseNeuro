#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 11:38:10 2025

@author: Ilse Klinkhamer
"""

import os
import re
from tkinter import filedialog, messagebox
from pathlib import Path
import tkinter as tk
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from get_dropbox_path_2 import get_dropbox_path
from collections import defaultdict
from plot_utils_IK import c4_colors_rgb, lighten, normalize_RGB_dict

LAST_FOLDER_FILE = Path.home() / ".last_folder_cell_type_counts.txt"  # Hidden file in home dir

def select_file_or_folder(prompt="Select a file or folder", start_path=None):
    # Determine the initial directory
    initial_dir = str(start_path) if start_path else (
        LAST_FOLDER_FILE.read_text() if LAST_FOLDER_FILE.exists() else str(Path.home())
    )

    root = tk.Tk()
    root.withdraw()

    # Ask the user what they want to select
    choice = messagebox.askquestion(
        "Selection Type", "Do you want to select a **file**? Click 'Yes'.\nTo select a **folder**, click 'No'."
    )

    if choice == 'yes':
        selected_path = filedialog.askopenfilename(title=f"{prompt}file", initialdir=initial_dir)
    else:
        selected_path = filedialog.askdirectory(title=f"{prompt}folder", initialdir=initial_dir)

    # Save path for next time
    if selected_path:
        path_obj = Path(selected_path)
        save_path = path_obj.parent if path_obj.is_file() else path_obj
        LAST_FOLDER_FILE.write_text(str(save_path))

    return selected_path if selected_path else None



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



def find_full_source_name(folder, source, channel):
    pattern = re.compile(rf"^{re.escape(source)}(.*?)_CH{channel}\.continuous$")
    for filename in os.listdir(folder):
        match = pattern.match(filename)
        if match:
            return source + match.group(1)
    return None  

def read_tsv(file_path):
    """Reads a TSV file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', index_col=0)

def plot_pie_charts(total_counts_per_cell_type, save_path):
    """Plot and save pie charts for all mouse groups in one figure."""
    num_groups = total_counts_per_cell_type.shape[0 if isinstance(total_counts_per_cell_type, pd.Series) else 1]
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
        ax.pie(counts.values, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)
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


def countCellTypesSession(results_file=None, session=None):
    if not results_file:
        results_file = select_file_or_folder(
            "Select the cell type ",
            start_path="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/Ilse/ephys/"
        )       
                      
    print(f"Processing folder: {results_file}")  
       
    session, mouse_name = find_session_and_mouse_name(results_file)   
    if os.path.exists(results_file):
        c4_results = pd.read_csv(results_file, sep='\t')
        c4_results['NeuronId'] = f"{session}_N" + c4_results['cluster_id'].astype(str)
    else:
        return
    
    c4_counts = c4_results['predicted_cell_type'].value_counts()
    counts_df = c4_counts.reset_index()
    counts_df.columns = ['cell_type', 'count']
    #counts_df = counts_df.T
    
    return c4_results, counts_df, c4_counts
            
            
def main(results_file_or_folder=None,contamination_ratio=0.1,
         confidence_ratio_threshold=1.5, select_folder=False):
    dropbox_path = get_dropbox_path()
    if not results_file_or_folder and select_folder:
        results_file_or_folder = select_file_or_folder(
            "Select the cell type",
            start_path="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/Ilse/ephys/"
        )  
    elif not results_file_or_folder:
        results_file_or_folder = os.path.join(dropbox_path, "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/Ilse/ephys/")
        

    print(dropbox_path)    
    counts_folder = os.path.join(dropbox_path, f"ContextMouseExperiments/Ilse/AnalysisOutput/Cell_type_counts/fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}/")
    os.makedirs(counts_folder, exist_ok=True)
    if os.path.isfile(results_file_or_folder):
        results_file = results_file_or_folder    
        session_name, mouse_name = find_session_and_mouse_name(results_file)
        c4_results, cell_type_counts, c4_counts = countCellTypesSession(results_file)
        pie_chart_path = os.path.join(counts_folder, f"{session_name}_pie_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
        plot_pie_charts(c4_counts, pie_chart_path)
        a = 1
    else:        
        mouse_directory = results_file_or_folder
        # Find all valid session subfolders
        subfolders = []
        sessions = []
        for f in os.listdir(mouse_directory):
            full_path = os.path.join(mouse_directory, f)
            if os.path.isdir(full_path):
                session_name, mouse_name = find_session_and_mouse_name(full_path)
                if not session_name or not mouse_name:
                    for f_sub in os.listdir(full_path):
                        full_sub_path = os.path.join(full_path, f_sub)
                        if os.path.isdir(full_sub_path):
                            session_name, mouse_name = find_session_and_mouse_name(full_sub_path)
                            if session_name is not None and mouse_name is not None:                    
                                subfolders.append(full_path)
                                sessions.append(session_name)
                else:                 
                    subfolders.append(full_path)
                    sessions.append(session_name)
        sessions.sort()
        subfolders.sort()    
        # List to collect neuron counts per session       
        session_counts = []
        
        for session_folder, session in zip(subfolders, sessions):
            print(f"\n--- Processing: {session_folder} ---")
            
            c4_results_folder = os.path.join(
                session_folder,
                "c4",
                f"c4_results_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}"
            )
            results_file = os.path.join(c4_results_folder, "cluster_predicted_cell_type.tsv")
            
            c4_results, cell_type_counts = countCellTypesSession(results_file, session)
            
            cell_type_names = ['GoC', 'MFB', 'MLI', 'PkC_cs', 'PkC_ss']
        
            if cell_type_counts is None:
                print(f"Warning: No data for session {session}. Adding zeros.")
                count_dict = {cell_type: 0 for cell_type in cell_type_names}
            else:
                if 'mouse_c4_df' not in locals():
                    mouse_c4_df = cell_type_counts
                else:
                    mouse_c4_df = pd.concat([mouse_c4_df, cell_type_counts], ignore_index=True)
                
                count_dict = {
                    cell_type: (cell_type_counts['predicted_cell_type'] == cell_type).sum()
                    for cell_type in cell_type_names
                }
        
            count_dict['session'] = session
            session_counts.append(count_dict)
        
        # Save summary
        summary_df = pd.DataFrame(session_counts)
        summary_df.to_csv(os.path.join(counts_folder, "neuron_counts_per_session.tsv"), sep='\t', index=False)
        
        # Store full data only if at least one session had real data
        if 'mouse_c4_df' in locals():
            mouse_df = mouse_c4_df
        else:
            mouse_df = pd.DataFrame()  # or None, depending on your downstream code
        
        # Summarize counts per cell type
        cell_type_names = ['GoC', 'MFB', 'MLI', 'PkC_cs', 'PkC_ss']
        count_dict = {cell_type: (cell_type_counts['predicted_cell_type'] == cell_type).sum()
                      for cell_type in cell_type_names}
        count_dict['session'] = session
        session_counts.append(count_dict)
        
        # Convert counts to DataFrame
        summary_df = pd.DataFrame(session_counts)
        
        # Save to TSV
        summary_file = "neuron_counts_per_session.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        
        # Store the full mouse data
        mouse_df = mouse_c4_df
        
        total_counts_per_cell_type = {}
        total_counts_per_mouse = {}
        cell_types = np.unique(mouse_df["predicted_cell_type"])
        cell_counts_group = defaultdict()
        
        for cell_type in cell_types:
            # General count for this cell type
            counts = (mouse_df["predicted_cell_type"] == cell_type).sum()
            # Initialize the cell_type if it doesn't exist
            if cell_type not in total_counts_per_cell_type:
                total_counts_per_cell_type[cell_type] = 0
            # Update the counts correctly
            total_counts_per_cell_type[cell_type] += counts                    
            #total_counts_per_mouse[mouse] = cell_type_counts.iloc[-2, :].astype(int)
            
    total_counts_per_cell_type = {k: int(v) for k, v in total_counts_per_cell_type.items()}

    # Debugging: check the final result to see if anything was overwritten incorrectly
    print(total_counts_per_cell_type)
    
    # Convert total_counts_per_cell_type to a DataFrame
    results_df = pd.DataFrame(total_counts_per_cell_type, index=[0])  # Transpose to get cell_types as index and priors as columns        
    
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
    
    for index, row in results_df.iterrows():
        labels = [f"{cell_type} ({count})" for cell_type, count in row.items()]
        GoCs += row['GoC']
        MFBs += row['MFB']
        PkC_css += row['PkC_cs']
        PkC_sss += row['PkC_ss']
        MLIs += row['MLI']
        
       
    results_all_mice = {
        'Index': ['GoC', 'MFB', 'MLI', 'PkC_cs', 'PkC_ss'],
        'AllMice': [GoCs, MFBs, MLIs, PkC_css, PkC_sss]
        }
    
    #results_all_mice = {k: int(v) for k, v in results_all_mice.items()}
    
    
    
    results_all_mice_df = pd.DataFrame(results_all_mice)
    results_all_mice_df = results_all_mice_df.set_index('Index')    

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
    per_mouse_chart_path = os.path.join(counts_folder, f"all_per_mouse_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.png")
    plot_per_mouse_group_bar_chart(total_counts_per_mouse, per_mouse_chart_path)
    print(f"Saved per-mouse-group bar chart to {per_mouse_chart_path}")
    per_mouse_chart_path = os.path.join(counts_folder, f"all_per_mouse_bar_chart_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}.eps")
    plot_per_mouse_group_bar_chart(total_counts_per_mouse, per_mouse_chart_path)
    """

                
if __name__ == "__main__":
    main()