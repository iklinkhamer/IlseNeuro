#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:57:24 2025

@author: Ilse Klinkhamer
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_mouse_groups():
    """Returns a dictionary containing categorized mouse groups."""
    return {
        "Switch": [
            "ReserveMouse3", "Dallas", "Flint", "Greene", "Houston", "Iowa", "Jackson",
            "Lincoln", "Newark", "Missouri?", "Pittsburg", "Queens?", "Orleans"
        ],
        "WideExperts": ["Reno", "Seattle", "Yosemite", "Zachary", "Kyiv", "Istanbul", "Copenhagen"],
        "Narrow": ["Rotterdam", "Tallinn", "Quimper", "Porto", "Lisbon", "Madrid"],
        "Bimodal": ["Uppsala", "Venice", "Willemstad", "Zurich", "York", "Xanthi"],
        "Naive": ["Ana1", "Ana2", "Ana3", "Ana4", "Ana5"]
    }

def read_tsv(file_path):
    """Reads a TSV file into a DataFrame."""
    return pd.read_csv(file_path, sep='\t', index_col=0)

def plot_pie_charts(total_counts_per_cell_type, save_path):
    """Plot and save pie charts for all mouse groups in one figure."""
    num_groups = len(total_counts_per_cell_type)
    fig, axes = plt.subplots(1, num_groups, figsize=(num_groups * 5, 5))
    
    if num_groups == 1:
        axes = [axes]
    
    color_map = {
        'PkC_cs': 'grey',
        'MLI': 'pink',
        'MFB': 'red',
        'GoC': 'green',
        'PkC_ss': 'blue'
    }
    
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
    color_map = {'PkC_cs': 'grey', 'MLI': 'pink', 'MFB': 'red', 'GoC': 'green', 'PkC_ss': 'blue'}
    
    
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

def main():
    counts_folder = "/home/no1/Documents/c4_neurons_temp_output/cell_type_classification_results/fnfpThreshold_0.1_confidenceRatio_1.5/"
    counts_file = "_cell_type_counts.tsv"
    mice_groups = get_mouse_groups()
    
    mice_to_analyze = {group: [] for group in mice_groups}
    
    # Check file existence
    for group, mice in mice_groups.items():
        for mouse in mice:
            file_path = os.path.join(counts_folder, f"{mouse}{counts_file}")
            if os.path.isfile(file_path):
                mice_to_analyze[group].append(mouse)
    
    total_counts_per_cell_type = {group: [] for group in mice_groups}
    
    # Read and aggregate data
    for group, mice in mice_to_analyze.items():
        for mouse in mice:
            file_path = os.path.join(counts_folder, f"{mouse}{counts_file}")
            cell_type_counts = read_tsv(file_path)
            
            if cell_type_counts.empty:
                continue
            
            total_counts_per_cell_type[group].append(cell_type_counts.iloc[:, 0:5].sum().astype(int))
    
    # Ensure only non-empty groups are processed
    results_df = pd.DataFrame({
        group: pd.concat(counts).groupby(level=0).sum().astype(int)
        for group, counts in total_counts_per_cell_type.items() if counts
    })
    
    # Save to file in counts_folder
    output_file = os.path.join(counts_folder, "cell_type_totals.tsv")
    results_df.to_csv(output_file, sep='\t')
    print(f"Saved cell type totals to {output_file}")
    
    # Generate and save plots
    pie_chart_path = os.path.join(counts_folder, "all_groups_pie_chart.png")
    bar_chart_path = os.path.join(counts_folder, "all_groups_bar_chart.png")
    
    plot_pie_charts(results_df, pie_chart_path)
    plot_grouped_bar_chart(results_df, bar_chart_path)
    
    print(f"Saved pie chart to {pie_chart_path}")
    print(f"Saved bar chart to {bar_chart_path}")
    
if __name__ == "__main__":
    main()





    