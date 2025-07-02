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
    """
    Prompts the user to select a file or folder. Remembers the last accessed directory.
    """
    init_dir = str(start_path or (LAST_FOLDER_FILE.read_text() if LAST_FOLDER_FILE.exists() else Path.home()))
    tk.Tk().withdraw()

    select = filedialog.askopenfilename if messagebox.askyesno("Selection Type", "Select **file**? (Yes)\nFolder? (No)") else filedialog.askdirectory
    path = Path(select(title=prompt, initialdir=init_dir))

    if path.exists():
        LAST_FOLDER_FILE.write_text(str(path.parent if path.is_file() else path))
        return str(path)


def extract_session_and_mouse_from_path(path):
    """
    Infers the session folder name and mouse name by scanning upward through the directory tree.
    Assumes session folders contain a date or datetime string in formats like YYYYMMDD or YYYYMMDDHHMMSS.

    Parameters:
        path (str or Path): A file or folder path somewhere within the session directory structure.

    Returns:
        tuple[str, str] or (None, None): A tuple with (session_name, mouse_name), where:
            - session_name is the name of the folder containing the date/datetime.
            - mouse_name is the parent folder of that session.
            Returns (None, None) if no valid session date is found.
    """
    path = Path(path).resolve()
    formats = ("%Y%m%d%H%M%S", "%Y%m%d%H%M", "%Y%m%d")
    for p in [path] + list(path.parents):
        name = p.name
        base = name.rsplit('_', 1)[0] if name.rsplit('_', 1)[-1].isdigit() else name
        digits = ''.join(filter(str.isdigit, base))
        for f in formats:
            try:
                datetime.strptime(digits, f)
                return name, p.parent.name
            except ValueError:
                pass
    return None, None

def path_contains_datetime(path):
    """
    Returns True if any folder in the path or its parents contains
    a datetime string matching one of the specified formats.
    Otherwise returns False.
    """
    path = Path(path).resolve()
    formats = ("%Y%m%d%H%M%S", "%Y%m%d%H%M", "%Y%m%d")
    pattern = re.compile(r'\d{8}(\d{2}(\d{2})?)?')  # matches 8, 10, or 14 digits

    for p in [path] + list(path.parents):
        name = p.name
        digit_groups = [m.group(0) for m in pattern.finditer(name)]
        for digits in digit_groups:
            for fmt in formats:
                if len(digits) == len(fmt.replace('%', '')):
                    try:
                        datetime.strptime(digits, fmt)
                        return True
                    except ValueError:
                        continue
    return False


def read_tsv(file_path):
    """Reads a TSV file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', index_col=0)

def plot_pie_charts(counts_per_group, save_paths):
    """Plot and save pie charts of cell type counts per group."""
    save_paths = [save_paths] if isinstance(save_paths, str) else save_paths    
    is_series = isinstance(counts_per_group, pd.Series)
    num = counts_per_group.shape[0 if is_series else 1]
    fig, axes = plt.subplots(1, num, figsize=(num * 5, 5))
    axes = [axes] if num == 1 else axes

    cmap = normalize_RGB_dict({**c4_colors_rgb(), 'PkC_cs': lighten(c4_colors_rgb()['PkC_cs'])})

    for ax, (group, counts) in zip(axes, counts_per_group.items()):
        labels = [f"{ct} ({n})" for ct, n in counts.items()]
        colors = [cmap.get(ct, 'gray') for ct in counts.index]
        ax.pie(counts.values, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)
        ax.set_title(group)

    plt.tight_layout()
    for path in save_paths:
        plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_grouped_bar_chart(counts, save_paths):
    """Plot and save a grouped bar chart of neuron counts per cell type across mice."""
    save_paths = [save_paths] if isinstance(save_paths, str) else save_paths
    df = pd.DataFrame(counts).T

    cmap = normalize_RGB_dict({**c4_colors_rgb(), 'PkC_cs': lighten(c4_colors_rgb()['PkC_cs'])})
    ax = df.plot(kind='bar', figsize=(12, 6), color=[cmap.get(col, 'gray') for col in df.columns])

    for c in ax.containers:
        ax.bar_label(c, label_type='edge')

    totals = df.sum(axis=1)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels([f"{m}\n({int(t)})" for m, t in zip(df.index, totals)], rotation=45)

    ax.set(xlabel="Mouse", ylabel="Neuron Count", title="Neuron Counts per Cell Type Across Mice")
    ax.legend(title="Cell Type")
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    for path in save_paths:
        plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()

    
def plot_per_mouse_group_bar_chart(counts, save_paths):
    """Plot and save a bar chart of cell type counts per mouse group."""
    save_paths = [save_paths] if isinstance(save_paths, str) else save_paths
    df = pd.DataFrame(counts).T.fillna(0).astype(int)

    cmap = normalize_RGB_dict({**c4_colors_rgb(), 'PkC_cs': lighten(c4_colors_rgb()['PkC_cs'])})
    fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.5), 8))
    df.plot(kind='bar', width=0.7, ax=ax, color=[cmap.get(c, 'gray') for c in df.columns])

    for bars in ax.containers:
        ax.bar_label(bars, label_type='edge')

    ax.set(xlabel="Mouse Name", ylabel="Neuron Count", title="Neuron Counts per Cell Type per Mouse Group")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df.index, rotation=60, ha='right', fontsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.legend(title="Cell Type", fontsize=10)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    for path in save_paths:
        plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def count_cell_types_session(results_file=None, session=None):
    """Load and count predicted cell types from a result file."""
    if not results_file:
        results_file = select_file_or_folder(
            "Select the cell type",
            start_path=os.path.join(get_dropbox_path(), "ContextMouseExperiments", "Ilse", "ephys")
        )

    print(f"Processing folder: {results_file}")
    
    if not os.path.exists(results_file):
        return None, None, None

    #session, mouse_name = extract_session_and_mouse_from_path(results_file)
    df = pd.read_csv(results_file, sep='\t')
    df['NeuronId'] = f"{session}_N" + df['cluster_id'].astype(str)

    counts = df['predicted_cell_type'].value_counts()
    counts_df = counts.rename_axis('cell_type').reset_index(name='count')

    return df, counts_df, counts
            
            
def main(results_file_or_folder=None, contamination_ratio=0.1, confidence_ratio_threshold=1.5, select_folder=False):
    dropbox_path = get_dropbox_path()
    base_path = os.path.join(dropbox_path, "ContextMouseExperiments", "Ilse", "ephys")
    if not results_file_or_folder:
        results_file_or_folder = select_file_or_folder("Select the cell type", start_path=base_path) if select_folder else base_path

    counts_folder = os.path.join(
        dropbox_path, "ContextMouseExperiments", "Ilse", "AnalysisOutput", "Cell_type_counts",
        f"fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}"
    )
    os.makedirs(counts_folder, exist_ok=True)

    if os.path.isfile(results_file_or_folder):
        session_name, _ = extract_session_and_mouse_from_path(results_file_or_folder)
        _, _, c4_counts = count_cell_types_session(results_file_or_folder)
        pie_path = os.path.join(
            counts_folder,
            f"total_cell_type_distribution_all_mice_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}.png"
        )
        plot_pie_charts(c4_counts, pie_path)
        return

    sessions, subfolders, mouse = [], [], []
    results_file_or_folder = Path(results_file_or_folder)
    directory = os.listdir(results_file_or_folder)
    directory.sort()
    for f in directory:
        full_path = results_file_or_folder / f
        if os.path.isdir(full_path):
            if path_contains_datetime(full_path):
                sessions.append(f)
                subfolders.append(full_path)     
                mouse.append(f.name)
            else:
                for f_sub in os.listdir(full_path):
                    full_sub_path = full_path/ f_sub
                    if os.path.isdir(full_sub_path):
                        subfolders.append(full_sub_path)
                        sessions.append(f_sub)
                        mouse.append(f)
                        
    session_counts, cell_type_names = [], ['GoC', 'MFB', 'MLI', 'PkC_cs', 'PkC_ss']
    mouse_c4_df = pd.DataFrame()
    for session_folder, session, mouse_name in zip(subfolders, sessions, mouse):
        print(f"\n--- Processing: {session_folder} ---")
        results_file = os.path.join(session_folder, "c4", f"c4_results_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}", "cluster_predicted_cell_type.tsv")
        _, cell_type_counts, _ = count_cell_types_session(results_file, session)

        if cell_type_counts is not None:
            mouse_c4_df = pd.concat([mouse_c4_df, cell_type_counts], ignore_index=True)
            counts = {ct: int(cell_type_counts.loc[cell_type_counts['cell_type'] == ct, 'count'].sum()) for ct in cell_type_names}
        else:
            print(f"Warning: No data for session {session}. Adding zeros.")
            counts = {ct: 0 for ct in cell_type_names}

        counts['session'] = session
        counts['mouse'] = mouse_name
        session_counts.append(counts)
        

    summary_df = pd.DataFrame(session_counts)
    summary_df.to_csv(os.path.join(counts_folder, f"neuron_counts_per_session_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}.tsv"), sep='\t', index=False)

    #results_df = mouse_c4_df.groupby("cell_type").size().reindex(cell_type_names, fill_value=0).to_frame().T
    #results_df = results_df.fillna(0).astype(int)
    #results_all_mice_df = results_df.copy()
    #results_all_mice_df.index = ['AllMice']
    #results_all_mice_df = results_all_mice_df.T
    
    total_counts = summary_df[cell_type_names].sum()
    total_counts_df = total_counts.to_frame().T
    total_counts_df.index = ['AllMice']


    save_types = [".png", ".eps"]
    for ext in save_types:
        pie_path = os.path.join(
            counts_folder,
            f"total_cell_type_distribution_all_mice_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}{ext}"
        )
        bar_path = os.path.join(
            counts_folder,
            f"total_cell_type_counts_all_mice_bar_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}{ext}"
        )
        plot_pie_charts(total_counts_df.T, pie_path)
        plot_grouped_bar_chart(total_counts_df.T, bar_path)

    total_counts_df.to_csv(os.path.join(
        counts_folder,
        f"cell_type_total_counts_all_mice_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}.tsv"
    ), sep='\t')

    #summary_df["base_name"] = summary_df["session"].str.extract(r'^([^0-9]*)')[0].str.rstrip('_')
    grouped_df = summary_df.groupby("mouse")[cell_type_names].sum().T
    filtered_grouped_df = grouped_df.loc[:, (grouped_df != 0).any(axis=0)]
    for ext in save_types:
        pie_path = os.path.join(
            counts_folder,
            f"cell_type_distribution_by_mouse_group_pie_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}{ext}"
        )
        bar_path = os.path.join(
            counts_folder,
            f"cell_type_counts_by_mouse_group_bar_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}{ext}"
        )
        
        plot_grouped_bar_chart(grouped_df, bar_path)
        plot_pie_charts(filtered_grouped_df, pie_path)

    print(f"Saved plots and results to {counts_folder}")
    
    grouped_df.to_csv(os.path.join(
        counts_folder,
        f"cell_type_totals_mice_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}.tsv"
    ), sep='\t')

                
if __name__ == "__main__":
    main()
    
    
"""
sessions, subfolders = [], []
for f in os.listdir(results_file_or_folder):
    full_path = os.path.join(results_file_or_folder, f)
    if os.path.isdir(full_path):
        session_name, mouse_name = extract_session_and_mouse_from_path(full_path)
        if session_name and mouse_name:
            subfolders.append(full_path)
            sessions.append(session_name)
        else:
            for f_sub in os.listdir(full_path):
                full_sub_path = os.path.join(full_path, f_sub)
                if os.path.isdir(full_sub_path):
                    session_name, mouse_name = extract_session_and_mouse_from_path(full_sub_path)
                    if session_name and mouse_name:
                        subfolders.append(full_sub_path)
                        sessions.append(session_name)
                        
                        
results_df.to_csv(os.path.join(
        counts_folder,
        f"cell_type_total_counts_combined_sessions_fpfn_{contamination_ratio}_conf_{confidence_ratio_threshold}.tsv"
    ), sep='\t')
"""