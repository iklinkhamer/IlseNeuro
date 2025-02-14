#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:19:06 2024

@author: Ilse Klinkhamer
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

# Load the .tsv file
def load_tsv(file_path):
    """Load the TSV file and return it as a DataFrame."""
    return pd.read_csv(file_path, sep='\t')

# Process the data to compute cell type fractions
def compute_fractions(data, cell_type_column):
    """Compute fractions of cell types in the given column."""
    cell_type_counts = data[cell_type_column].value_counts()
    total = cell_type_counts.sum()
    fractions = cell_type_counts / total
    return cell_type_counts, fractions

# Plot the fractions as a pie chart
def plot_pie_chart(cell_type_counts, fractions, title="Cell Type Distribution"):
    """Plot a pie chart of the cell type fractions."""
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

    wedges, texts, autotexts = plt.pie(
        fractions, labels=labels, colors=colors,
        autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100. * sum(cell_type_counts)))})", startangle=140
    )
    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_fontsize(10)
    plt.title(title)
    plt.show()

# Main function
def main():
    mouse_name = input("Enter the mouse name: ")
    directory = f"/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/{mouse_name}/"

    try:
        # Get all folders for the mouse
        mouse_folders = [
            folder for folder in os.listdir(directory)
            if os.path.isdir(os.path.join(directory, folder)) 
            and mouse_name in folder 
            and "copy" not in folder
        ]
        mouse_folders.sort()

        for folder in mouse_folders:
            folder_path = os.path.join(directory, folder, "c4")
            tsv_files = [f for f in os.listdir(folder_path) if f == 'cluster_predicted_cell_type.tsv']

            for tsv_file in tsv_files:
                file_path = os.path.join(folder_path, tsv_file)
                print(f"Processing file: {file_path}")

                # Load data
                data = load_tsv(file_path)

                # Ensure necessary columns exist
                required_columns = ['cluster_id', 'predicted_cell_type']
                if not all(col in data.columns for col in required_columns):
                    print(f"Skipping file {file_path}: Missing required columns.")
                    continue

                # Compute counts and fractions
                cell_type_counts, fractions = compute_fractions(data, 'predicted_cell_type')

                # Filter for specific cell types
                relevant_cell_types = ['PkC_cs', 'MLI', 'MFB', 'GoC', 'PkC_ss']
                filtered_counts = cell_type_counts[cell_type_counts.index.isin(relevant_cell_types)]
                filtered_fractions = fractions[fractions.index.isin(relevant_cell_types)]

                # Plot pie chart
                plot_pie_chart(filtered_counts, filtered_fractions, title=f"Cell Type Distribution for {folder}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
