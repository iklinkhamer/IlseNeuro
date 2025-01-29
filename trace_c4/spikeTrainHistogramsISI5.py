#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:52:59 2025

@author: Ilse Klinkhamer
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import uniform_filter1d


def load_spike_times(file_path):
    return np.loadtxt(file_path)  # Adjust based on your file format (e.g., CSV, MAT, etc.)


def compute_3d_acg(spike_times, bin_edges, smooth_window=250):
    # Compute interspike intervals (ISIs)
    isi = np.diff(spike_times)
    inst_rate = 1 / isi  # Instantaneous firing rate

    # Smooth firing rate with a boxcar filter (250 ms width)
    smooth_inst_rate = uniform_filter1d(inst_rate, size=smooth_window)

    # Stratify firing rates into 10 deciles
    deciles = np.percentile(smooth_inst_rate, np.linspace(0, 100, 11))
    spike_bins = np.digitize(smooth_inst_rate, deciles) - 1  # Assign each spike to a decile bin

    # Compute 2D-ACGs for each decile
    time_lags = (bin_edges[:-1] + bin_edges[1:]) / 2
    acgs = np.zeros((10, len(time_lags)))

    for d in range(10):
        spike_indices = np.where(spike_bins == d)[0]
        for i in spike_indices:
            spike_diffs = spike_times[i] - spike_times
            hist, _ = np.histogram(spike_diffs, bins=bin_edges)
            acgs[d] += hist

    # Normalize by bin width (spikes/s)
    acgs /= np.diff(bin_edges)

    return acgs, time_lags


def plot_3d_acg(acgs, time_lags):
    plt.figure(figsize=(10, 6))
    sns.heatmap(acgs, xticklabels=np.round(time_lags, 2), cmap='viridis')
    plt.xlabel('Time from spike (ms)')
    plt.ylabel('Firing Rate Decile')
    plt.title('3D Autocorrelogram')
    plt.show()


# Example Usage
file_path = "spike_times.txt"  # Replace with actual file path
spike_times = load_spike_times(file_path)
bin_edges = np.geomspace(1, 1000, num=50)  # Log-spaced bins
acgs, time_lags = compute_3d_acg(spike_times, bin_edges)
plot_3d_acg(acgs, time_lags)
