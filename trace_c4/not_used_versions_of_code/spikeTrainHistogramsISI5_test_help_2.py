#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:38:56 2025

@author: Ilse Klinkhamer
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
from utils_IK import load_data
import gc  # Garbage collection module
from matplotlib.ticker import ScalarFormatter

def autocorrelation(x):
    result = np.correlate(x, x, mode='full')
    return result   #[result.size // 2:]

sampling_rate = 30000
max_lag_ms = 1000  # maximum lag in milliseconds

# Compute the maximum lag in samples
max_lag_samples = int(max_lag_ms / 1000 * sampling_rate)  # 1000 ms = 1 second


# Example Usage
directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909142859/c4/"
#directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909145528/c4/"
spike_times_file_path = f"{directory}spike_times.npy"
spike_clusters_file_path = f"{directory}spike_clusters.npy"
isi_values_path = f"{directory}isi_values/spike_differences.npy"

spike_times = load_data(spike_times_file_path)
spike_times = spike_times / sampling_rate
spike_clusters = load_data(spike_clusters_file_path)
#isi_values = load_data(isi_values_path)
bin_edges = np.geomspace(1, 1000, num=50)  # Log-spaced bins

# Precompute unique clusters
unique_clusters = np.unique(spike_clusters)
n_clusters = len(unique_clusters)

# Initialize a list to store differences, but process and save them incrementally
spike_diffs = [[] for _ in range(n_clusters)]

start_time = time.time()
# Efficiently process pairwise differences and save to disk
for c in tqdm(unique_clusters, desc="Processing spikes"):
    # Find the indices for the current cluster
    cluster_indices = np.where(spike_clusters == c)[0]
    
    # Extract the spike times for the current cluster
    spike_times_cluster = spike_times[cluster_indices]
    
    # # Calculate pairwise differences using broadcasting
    # diffs_all = spike_times_cluster[:, None] - spike_times_cluster  # Shape: (n_spikes, n_spikes)
    
    # # Use np.tril_indices to get only the lower triangle (without the diagonal)
    # lower_triangle_indices = np.tril_indices(len(diffs_all), -1)
    # diffs_all = diffs_all[lower_triangle_indices]
    
    # plt.hist(diffs_all, bins=50, color='b', alpha=0.7)
    # plt.xlabel("Time (s)")
    # plt.ylabel("Spike Count")
    # plt.title("Spike Time Histogram")
    # plt.show()
    
    diffs = np.diff(spike_times_cluster.flatten())
    
    plt.hist(diffs, bins=50, color='b', alpha=0.7)
    plt.xlabel("Time (s)")
    plt.ylabel("Spike Count")
    plt.title("Spike Time Histogram")
    plt.show()
    
    y_max = 20
    diffs_ms = diffs*1000
    diffs_ms = diffs_ms[diffs_ms <= y_max]
    plt.hist(diffs_ms, bins=50, color='b', alpha=0.7)
    plt.xlim(0, y_max)
    plt.xlabel("Time (ms)")
    plt.ylabel("Spike Count")
    plt.title("Spike Time Histogram")
    plt.show()
    
    # Calculate autocorrelation for the current cluster
    auto_corr = autocorrelation(diffs.flatten())
    #auto_corr = autocorrelation(spike_times_cluster.flatten())
    
    # Normalize the autocorrelation by bin width
    # The bin width is the time difference between each consecutive lag, which corresponds to 1/sampling_rate in seconds
    bin_width_seconds = 1 / sampling_rate  # In seconds
    
    # Normalize by dividing by the bin width (convert to spikes per second)
    auto_corr_normalized = auto_corr / bin_width_seconds
    # Compute the range of time lags (keeping both negative and positive shifts)
    num_samples = len(auto_corr)
    time_lags_ms = np.arange(-num_samples // 2, num_samples // 2) * (1000 / sampling_rate)

   
    plt.plot(time_lags_ms, -auto_corr_normalized)
    plt.xlim(-50, 50)  # Ensure proper mirroring in x-axis
    plt.xscale('symlog')
    
    plt.title(f"Cluster: {c} Autocorrelation (mirrored, up to ±1000 ms)")
    plt.xlabel("Lag (ms)")
    plt.ylabel("Spikes / s")  # Updated ylabel to spikes per second
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.show()
    
  
    # Explicitly delete variables no longer needed and collect garbage to free memory
    del diffs, auto_corr, spike_times_cluster
    gc.collect()  # Garbage collection
    
# End timing
end_time = time.time()
elapsed_time = end_time - start_time

# Print elapsed time for autocorrelation calculation
print(f"Cluster {c} - Elapsed time: {elapsed_time:.4f} seconds")
  
