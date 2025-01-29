#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:52:59 2025

@author: Ilse Klinkhamer
"""

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from scipy.interpolate import make_interp_spline

def load_data(spike_times_file_path):
    """
    Load spike times from a file.

    Parameters:
    - spike_times_file_path (str): Path to the file containing spike times.

    Returns:
    - spike_times (np.array): Array of spike times.
    """
    # Determine the file extension
    _, file_extension = os.path.splitext(spike_times_file_path)

    if file_extension == ".txt":
        spike_times = np.loadtxt(spike_times_file_path)  # Load from a text file
    elif file_extension == ".csv":
        spike_times = np.genfromtxt(spike_times_file_path, delimiter=",")  # Load from CSV
    elif file_extension == ".npy":
        spike_times = np.load(spike_times_file_path)  # Load from NumPy binary format
    else:
        raise ValueError("Unsupported file format. Please use .txt, .csv, or .npy files.")

    return spike_times

def compute_autocorrelogram(spike_times, spike_clusters, bin_width=1.0, max_lag=50, save_diffs_path=None):
    """
    Compute and plot the autocorrelogram for spike times.

    Parameters:
    - spike_times (list or np.array): Array of spike times (in ms or s).
    - bin_width (float): Width of bins for the histogram (same units as spike_times).
    - max_lag (float): Maximum lag to consider for autocorrelogram (same units as spike_times).

    Returns:
    - lags (np.array): Lag times for the autocorrelogram.
    - autocorr (np.array): Counts for each lag bin.
    """
    spike_diffs = []

    # Check if saved spike differences exist
    if save_diffs_path and os.path.exists(save_diffs_path):
        print(f"Loading previously saved spike differences from: {save_diffs_path}")
        spike_diffs = np.load(save_diffs_path)
    else:
        # Calculate pairwise differences between spike times
        for c in tqdm(range(len(spike_clusters)), desc="Processing spikes"):
            for i in range(len(spike_times)):
                spike_times_cluster = spike_times[spike_clusters==c]
                diffs = spike_times_cluster - spike_times_cluster[i]
                spike_diffs[:,c].extend(diffs[np.abs(diffs) <= max_lag])  # Keep only diffs within max_lag

        # Exclude zero-lag (spike with itself)
        spike_diffs = np.array(spike_diffs)
        spike_diffs = spike_diffs[spike_diffs != 0]

        # Save the spike differences if a path is provided
        if save_diffs_path:
            np.save(save_diffs_path, spike_diffs)
            print(f"Spike differences saved to: {save_diffs_path}")

    # Create histogram
    bins = np.arange(-max_lag, max_lag + bin_width, bin_width)
    autocorr, edges = np.histogram(spike_diffs, bins=bins)

    # Normalize to account for the number of spikes
    autocorr = autocorr / len(spike_times)

    # Prepare lag times
    lags = edges[:-1] + bin_width / 2

    # Plot the autocorrelogram
    plt.figure(figsize=(8, 4))
    plt.bar(lags, autocorr, width=bin_width, color='gray', edgecolor='black', alpha=0.7, label="Histogram")

    # Create a smooth line plot
    lags_mid = lags[:-1] + (lags[1] - lags[0]) / 2
    spline = make_interp_spline(lags, autocorr, k=3)
    lags_smooth = np.linspace(lags[0], lags[-1], 500)
    autocorr_smooth = spline(lags_smooth)
    plt.plot(lags_smooth, autocorr_smooth, color='blue', label="Smoothed Line")

    plt.xlabel("Lag (ms or s)")
    plt.ylabel("Normalized Count")
    plt.title("Autocorrelogram of Spike Times")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.show(block=True)

    return lags, autocorr

# Example Usage
if __name__ == "__main__":
    # Specify the path to the spike times file
    directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909142859/c4/"  # Update with your file path

    spike_times_file_path = f"{directory}spike_times.npy"
    spike_clusters_file_path = f"{directory}spike_clusters.npy"
    try:
        # Load spike times
        spike_times = load_data(spike_times_file_path)
        spike_clusters = load_data(spike_clusters_file_path)

        # Path to save the spike differences
        save_diffs_path = f"{directory}isi_values/spike_differences_clusters.npy"

        # Create autocorrelogram and save spike differences
        compute_autocorrelogram(spike_times, spike_clusters, bin_width=5, max_lag=100, save_diffs_path=save_diffs_path)

    except Exception as e:
            print(f"Error: {e}")