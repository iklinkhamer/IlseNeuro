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


def load_data(spike_times_file_path):
    """
    Load spike times from a file.

    Parameters:
    - spike_times_file_path (str): Path to the file containing spike times.

    Returns:
    - spike_times (np.array): Array of spike times.
    """
    _, file_extension = os.path.splitext(spike_times_file_path)

    if file_extension == ".txt":
        spike_times = np.loadtxt(spike_times_file_path)
    elif file_extension == ".csv":
        spike_times = np.genfromtxt(spike_times_file_path, delimiter=",")
    elif file_extension == ".npy":
        spike_times = np.load(spike_times_file_path)
    else:
        raise ValueError("Unsupported file format. Please use .txt, .csv, or .npy files.")

    return spike_times


def compute_autocorrelogram(spike_times, spike_clusters, bin_width=1.0, max_lag=50):
    """
    Compute and plot the autocorrelogram for spike times using np.corrcoef.

    Parameters:
    - spike_times (np.array): Array of spike times.
    - spike_clusters (np.array): Cluster labels for each spike.
    - bin_width (float): Width of bins for the histogram.
    - max_lag (float): Maximum lag to consider for autocorrelogram.
    """
    unique_clusters = np.unique(spike_clusters)

    for c in unique_clusters:
        spike_times_cluster = spike_times[spike_clusters == c]
        isi = np.diff(spike_times_cluster)

        # Compute correlation coefficients for different lags
        autocorr_values = []
        lags = np.arange(-max_lag, max_lag + bin_width, bin_width)

        for lag in lags:
            shifted_isi = np.roll(isi, lag)
            corr = \
            np.corrcoef(isi[:-lag] if lag > 0 else isi[lag:], shifted_isi[:-lag] if lag > 0 else shifted_isi[lag:])[
                0, 1]
            autocorr_values.append(corr)

        # Plot autocorrelation
        plt.figure(figsize=(8, 4))
        plt.plot(lags, autocorr_values, marker='o', linestyle='-', color='blue', label="Autocorrelation")

        plt.xlabel("Lag (ms or s)")
        plt.ylabel("Correlation Coefficient")
        plt.title(f"Autocorrelogram for Cluster {c}")
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend()
        plt.show(block=True)


if __name__ == "__main__":
    directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909142859/c4/"

    spike_times_file_path = f"{directory}spike_times.npy"
    spike_clusters_file_path = f"{directory}spike_clusters.npy"

    try:
        spike_times = load_data(spike_times_file_path)
        spike_clusters = load_data(spike_clusters_file_path)

        compute_autocorrelogram(spike_times, spike_clusters, bin_width=5, max_lag=100)
    except Exception as e:
        print(f"Error: {e}")
