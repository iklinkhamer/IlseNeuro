#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:49:12 2025

@author: Ilse Klinkhamer
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.io import loadmat
from pathlib import Path
from get_dropbox_path import get_dropbox_path

def plot_trial_spikes(trial, spike_times, y=None, ax=None, annotation_fields=None, spike_height=None):
    if y is None:
        y = trial.index
    if ax is None:
        fig, ax = plt.subplots()
    if annotation_fields is None:
        annotation_fields = []
    if spike_height is None:
        spike_height = default_spike_height()
    
    if len(spike_times) == 0:
        return []
    
    t_vals_matrix = np.array([spike_times, spike_times, np.full_like(spike_times, np.nan)])
    t_vals = t_vals_matrix.flatten(order='F')
    y_vals = np.tile([y, y + spike_height, np.nan], len(spike_times))
    
    h, = ax.plot(1e3 * t_vals, y_vals, color='black', linewidth=0.5)
    
    for field in annotation_fields:
        values = np.full_like(t_vals, str(getattr(trial, field)), dtype=object)
        ax.text(1e3 * spike_times[0], y, values[0], fontsize=8, verticalalignment='bottom')
    
    return h

def default_spike_height():
    return 1.0  # Set this to an appropriate default value

def load_and_align_spikes(spike_file, event_file, cluster_file, session_file, pre_event_ms=750, post_event_ms=750, sampling_rate=30000):
    spike_times = np.load(spike_file) / sampling_rate  # Convert frames to seconds
    event_data = loadmat(event_file)  # Load the .mat file
    event_times = event_data['StampCS'].flatten() / sampling_rate  # Convert event times
    #session_info = loadmat(session_file)
    #trial_types = session_info['sess']['type'].flatten()
    
    
    cluster_data = np.load(cluster_file)  # Load cluster IDs
    unique_clusters = np.unique(cluster_data)  # Get unique cluster IDs
    unique_clusters = np.sort(unique_clusters)
    
    aligned_spikes = {cluster: [] for cluster in unique_clusters}
    
    for cluster in unique_clusters:
        cluster_indices = np.where(cluster_data == cluster)[0]  # Find spike indices for this cluster
        cluster_spikes = spike_times[cluster_indices]  # Get spike times
        
        for event_time in event_times:
            start_time = event_time - (pre_event_ms / 1000)
            end_time = event_time + (post_event_ms / 1000)
            aligned_spikes[cluster].append(cluster_spikes[(cluster_spikes >= start_time) & (cluster_spikes <= end_time)] - event_time)

    return aligned_spikes  # Ensure this is a dictionary

def main():
    mouse_session = "Venice_20240523111730"
    directory=os.path.join(get_dropbox_path(),"/ExperimentOutput/Ephys4Trace1/MainFolder/", "Venice/Venice_20240523111730/Extraction2Bin/")
    spike_file = f'{directory}spike_times.npy'
    event_file = f'{directory}EventStamps.mat'
    cluster_file = f'{directory}spike_clusters.npy'
    short_dir = Path(directory).parent
    session_file = f'{short_dir}/EphysSession_{mouse_session}.mat'
    save_dir = '/home/no1/Documents/SpikePlots'
    os.makedirs(save_dir, exist_ok=True)
    
    aligned_spikes = load_and_align_spikes(spike_file, event_file, cluster_file, session_file)
    
    session_info = loadmat(session_file)
    trial_types = session_info['sess']['type'].flatten()
    trial_types_unique = np.unique(trial_types)
    indices = np.lexsort((np.arange(len(trial_types)), trial_types))
    sorted_trials = trial_types[indices]

    
    for cluster, spikes_per_trial in aligned_spikes.items():
        fig, ax = plt.subplots()
        for i, spikes in enumerate(spikes_per_trial):
            plot_trial_spikes(trial=type('Trial', (object,), {'index': np.where(indices == i)[0][0]})(), spike_times=spikes, y=np.where(indices == i)[0][0], ax=ax)
        
        ax.set_xlim([-300, 800])
        ax.set_xticks([0, 200, 275, 350, 425, 500])
        ax.set_xticklabels(["0", "200", "275", "350", "425", "500"])
        for x in [0, 200, 275, 350, 425, 500]:
            ax.axvline(x, color='black', linestyle='dashed', linewidth=0.75)
        
        ax.set_ylim(bottom=0)  # Ensure no extra space at the bottom of y-axis
        ax.spines['left'].set_bounds(0, ax.get_ylim()[1])
        ax.spines['bottom'].set_bounds(-300, 800)
        
        plt.xlabel('Time (ms)')
        plt.ylabel('Trial Index')
        plt.title(f'Aligned Spike Times - {mouse_session}_{cluster}')
        plt.savefig(os.path.join(save_dir, f'Spike_raster_python_phy_{mouse_session}_{cluster}.png'))
        plt.close(fig)
    
if __name__ == "__main__":
    main()








