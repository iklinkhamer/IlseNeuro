#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:52:59 2025

@author: Ilse Klinkhamer
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.signal import correlate

plot_histogram = False
plot_boxplot = False


def get_isi(data_path, isi_save_path, sampling_rate, threshold):
    if os.path.exists(isi_save_path):
        print("Loading precomputed ISIs...")
        isi = np.loadtxt(isi_save_path, delimiter=",")
    else:
        print("Computing ISIs...")
        # Load raw data
        data = np.fromfile(data_path, dtype='int16')
        
        # Detect spikes
        spike_indices = np.where(data > threshold)[0]
        spike_times = spike_indices / sampling_rate
        
        # Compute ISIs
        isi = np.diff(spike_times)
        
        # Save ISIs
        np.savetxt(isi_save_path, isi, delimiter=",")
        print(f"ISIs saved to {isi_save_path}")
    
    return isi

# Example usage
directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909142859/c4/"
data_path = f"{directory}continuous/Data_AP_LFP/continuous.dat"
isi_save_path = f"{directory}isi_values/isi_values.csv"
sampling_rate = 30000  # Adjust to your setup

# Load the data once for both threshold calculation and ISI computation
data = np.fromfile(data_path, dtype='int16')
threshold = 5 * np.std(data)  # Calculate threshold based on standard deviation of the raw data

# Get ISI values
isi = get_isi(data_path, isi_save_path, sampling_rate, threshold)

# Count ISIs greater than and less than or equal to 0.005
greater_than_005 = np.sum(isi > 0.005)  # Count ISIs greater than 0.005
less_than_or_equal_005 = np.sum(isi <= 0.005)  # Count ISIs less than or equal to 0.005

print(f"Number of ISIs > 0.005: {greater_than_005}")
print(f"Number of ISIs <= 0.005: {less_than_or_equal_005}")

# Define a reasonable range for the histogram
min_isi = np.min(isi)
max_isi = np.max(isi)
bin_width = 0.001  # Set a small bin width (e.g., 1 ms)
bins = np.arange(min_isi, max_isi + bin_width, bin_width)  # Dynamically define bins

# Define the ISI range and bin count
max_display_range = 0.002  # Focus on ISIs up to 5 ms
num_bins = 200  # Reasonable number of bins

# Filter ISI values for the desired range
filtered_isi = isi[isi <= max_display_range]

# Use np.histogram for faster binning
hist_counts, bin_edges = np.histogram(filtered_isi, bins=num_bins, range=(0, max_display_range))

if plot_histogram:
    # Plot the histogram
    plt.figure(figsize=(8, 6))
    plt.bar(bin_edges[:-1], hist_counts, width=np.diff(bin_edges), color='blue', alpha=0.7, edgecolor='black')
    plt.xlabel('Interspike Interval (s)')
    plt.ylabel('Count')
    plt.title('Histogram of Interspike Intervals (ISIs)')
    plt.xlim(0, max_display_range)  # Limit x-axis to focus on ISIs up to 5 ms
    plt.tight_layout()
    plt.show()

if plot_boxplot:
    # Boxplot to compare datasets
    plt.figure(figsize=(8, 6))
    plt.boxplot([isi, isi], tick_labels=["Data", "Same Data"])
    plt.ylabel("Interspike Interval (s)")
    plt.title("Boxplot of ISI Spread")
    #plt.ylim(0, 5)  # Focus on a specific range
    plt.show()



def autocorrelation_from_isi(isi, max_lag, sampling_rate):
    """
    Calculate the autocorrelation directly from ISI data using np.correlate.
    :param isi: Interspike intervals (ISI) in seconds.
    :param max_lag: Maximum lag for which to compute the autocorrelation.
    :param sampling_rate: The sampling rate in Hz (samples per second).
    :return: Normalized autocorrelation array.
    """
    # Reconstruct the spike train from ISI values
    # Cumulatively sum ISIs to get spike times
    spike_times = np.cumsum(isi)
    
    # Convert spike times to indices based on the sampling rate
    spike_indices = np.floor(spike_times * sampling_rate).astype(int)
    
    # Length of the spike train
    total_samples = spike_indices[-1] + 1  # Total samples until the last spike
    
    # Create an array for the spike train (binary)
    spike_train = np.zeros(total_samples)
    spike_train[spike_indices] = 1  # Mark spikes at their respective indices

    # Compute autocorrelation using np.correlate
    corr = correlate(spike_train, spike_train, mode='full')
    
    # Since we only want the positive lags, we select the second half of the result
    corr = corr[corr.size // 2:]
    
    # Normalize the autocorrelation
    corr /= np.max(corr)
    
    # Convert autocorrelation to spikes per second (Hz)
    corr *= sampling_rate
    
    # Limit to the specified max_lag
    corr = corr[:max_lag]
    
    return corr

# Example usage
# Assuming 'isi' is your array of ISI values (in seconds)
# isi = np.array([...])  # Replace with your actual ISI values

# Set sampling rate and maximum lag
sampling_rate = 30000  # Hz (samples per second)
max_lag = 100  # Maximum lag in samples

# Calculate the autocorrelation from ISI values using np.correlate
auto_corr = autocorrelation_from_isi(isi, max_lag, sampling_rate)

# Plot the autocorrelogram
lags_ms = np.arange(0, max_lag) * (1000 / sampling_rate)  # Convert lag to milliseconds
plt.figure(figsize=(10, 6))
plt.plot(lags_ms, auto_corr, color='blue')
plt.title('Autocorrelogram from ISI Data (Spikes per Second)')
plt.xlabel('Lag (in milliseconds)')
plt.ylabel('Autocorrelation (Spikes per Second)')
plt.grid(True)
plt.show()




"""
import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.signal import correlate

plot_histogam = False
plot_boxplot = False


def get_isi(data_path, isi_save_path, sampling_rate, threshold):
    if os.path.exists(isi_save_path):
        print("Loading precomputed ISIs...")
        isi = np.loadtxt(isi_save_path, delimiter=",")
    else:
        print("Computing ISIs...")
        # Load raw data
        data = np.fromfile(data_path, dtype='int16')
        
        # Detect spikes
        spike_indices = np.where(data > threshold)[0]
        spike_times = spike_indices / sampling_rate
        
        # Compute ISIs
        isi = np.diff(spike_times)
        
        # Save ISIs
        np.savetxt(isi_save_path, isi, delimiter=",")
        print(f"ISIs saved to {isi_save_path}")
    
    return isi

# Example usage
directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909142859/c4/"
data_path = f"{directory}continuous/Data_AP_LFP/continuous.dat"
isi_save_path = f"{directory}isi_values/isi_values.csv"
sampling_rate = 30000  # Adjust to your setup
threshold = 5 * np.std(np.fromfile(data_path, dtype='int16'))  # Example threshold

isi = get_isi(data_path, isi_save_path, sampling_rate, threshold)

# Assuming 'isi' is your array of Interspike Intervals
greater_than_005 = np.sum(isi > 0.005)  # Count ISIs greater than 0.005
less_than_or_equal_005 = np.sum(isi <= 0.005)  # Count ISIs less than or equal to 0.005

print(f"Number of ISIs > 0.005: {greater_than_005}")
print(f"Number of ISIs <= 0.005: {less_than_or_equal_005}")

# Define a reasonable range for the histogram
min_isi = np.min(isi)
max_isi = np.max(isi)
bin_width = 0.001  # Set a small bin width (e.g., 1 ms)
bins = np.arange(min_isi, max_isi + bin_width, bin_width)  # Dynamically define bins

# Define the ISI range and bin count
max_display_range = 0.002  # Focus on ISIs up to 5 ms
num_bins = 200  # Reasonable number of bins

# Filter ISI values for the desired range
filtered_isi = isi[isi <= max_display_range]

# Use np.histogram for faster binning
hist_counts, bin_edges = np.histogram(filtered_isi, bins=num_bins, range=(0, max_display_range))

if plot_histogam:
    # Plot the histogram
    plt.figure(figsize=(8, 6))
    plt.bar(bin_edges[:-1], hist_counts, width=np.diff(bin_edges), color='blue', alpha=0.7, edgecolor='black')
    plt.xlabel('Interspike Interval (s)')
    plt.ylabel('Count')
    plt.title('Histogram of Interspike Intervals (ISIs)')
    plt.xlim(0, max_display_range)  # Limit x-axis to focus on ISIs up to 5 ms
    plt.tight_layout()
    plt.show()

if plot_boxplot:
    # Boxplot to compare datasets
    plt.figure(figsize=(8, 6))
    plt.boxplot([isi, isi], tick_labels=["Data", "Same Data"])
    plt.ylabel("Interspike Interval (s)")
    plt.title("Boxplot of ISI Spread")
    #plt.ylim(0, 5)  # Focus on a specific range
    plt.show()


# Calculate the autocorrelation function with a progress bar
def autocorrelation(x, max_lag=100):
    #Compute the autocorrelation of a 1D signal with a progress bar.
    corr = np.zeros(max_lag)  # Initialize the array for the autocorrelation
    N = len(x)
    for lag in tqdm(range(1, max_lag), desc="Calculating Autocorrelation"):
        # Compute the autocorrelation for each lag
        corr[lag] = np.sum(x[:N-lag] * x[lag:])  # Ensure no overlap by limiting the range
    corr /= np.max(corr)  # Normalize to max correlation value (for visualization)
    return corr

# Calculate autocorrelation of ISIs (limiting to first 100 lags)
isi_autocorr = autocorrelation(isi, max_lag=1000)

# Plot the autocorrelogram
plt.figure(figsize=(10, 6))
plt.plot(isi_autocorr, color='blue')
plt.title('Autocorrelogram of ISIs')
plt.xlabel('Lag')
plt.ylabel('Normalized Autocorrelation')
plt.grid(True)
plt.show()



def autocorrelation(x, max_lag=100, sampling_rate=30000):
    #Compute the autocorrelation of a 1D signal with a progress bar.
    N = len(x)
    max_time = np.sum(x)  # Total time from all ISIs
    total_samples = int(np.floor(max_time * sampling_rate))  # Total samples in spike train
    
    # Limit total_samples to a practical size (e.g., first 10 minutes of data)
    max_samples_limit = 10 * 60 * sampling_rate  # Limit to 10 minutes
    total_samples = min(total_samples, max_samples_limit)  # Apply the limit
    
    # Scale the ISIs into spike trains
    spike_train = np.zeros(total_samples)  # Create a binary spike train array of sufficient size
    spike_indices = np.floor(np.cumsum(x) * sampling_rate).astype(int)  # Calculate spike indices
    
    # Ensure that the spike indices do not exceed the total size
    spike_indices = spike_indices[spike_indices < total_samples]  # Clip to valid range
    spike_train[spike_indices] = 1  # Mark spikes at times given by ISI

    # Compute autocorrelation using np.correlate for efficiency
    corr = correlate(spike_train, spike_train, mode='full')[total_samples-1:total_samples+max_lag-1]
    
    # Normalize and convert to spikes per second (Hz)
    corr /= np.max(corr)  # Normalize to the max value
    corr *= sampling_rate  # Convert to spikes per second (Hz)

    return corr

# Set the max_lag value before calling the autocorrelation function
max_lag = 100  # or any value you want to set

# Ensure the autocorrelation function is called and assigned
isi_autocorr = autocorrelation(isi, max_lag=max_lag, sampling_rate=sampling_rate)

# Convert lag from samples to milliseconds
lags_ms = np.arange(0, max_lag) * (1000 / sampling_rate)  # Lag in milliseconds, start from 0

# Plot the autocorrelogram in spikes per second (Hz)
plt.figure(figsize=(10, 6))
plt.plot(lags_ms, isi_autocorr[:max_lag], color='blue')  # Plot with lag in ms, use the correct slice
plt.title('Autocorrelogram of ISIs (Spikes per Second)')
plt.xlabel('Lag (in milliseconds)')
plt.ylabel('Autocorrelation (Spikes per Second)')
plt.grid(True)
plt.show()




"""







"""
# Parameters for the .dat file
sampling_rate = 30000  # Sampling rate in Hz (adjust as per your system)
dtype = 'int16'        # Data type of the .dat file (usually int16 or int32)
n_channels = 1         # Number of channels in the .dat file

# Load the continuous.dat file
data_path = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909140005/c4/continuous/Data_AP_LFP/continuous.dat"
data = np.fromfile(data_path, dtype=dtype).reshape(-1, n_channels)


# Set a threshold (adjust based on your data)
threshold = 5 * np.std(data)  # 5x standard deviation

# Find spike indices (samples where the signal crosses the threshold)
spike_indices = np.where(data > threshold)[0]

# Convert spike indices to spike times in seconds
spike_times = spike_indices / sampling_rate

# Calculate ISIs (differences between consecutive spike times)
isi = np.diff(spike_times)

import matplotlib.pyplot as plt

# Plot the ISI histogram
plt.figure(figsize=(8, 6))
#plt.hist(isi, bins=50, color='blue', alpha=0.7, edgecolor='black')
plt.hist(isi)
plt.xlim((0, 50))
plt.xlabel('Interspike Interval (s)')
plt.ylabel('Count')
plt.title('Histogram of Interspike Intervals (ISIs)')
plt.show()
"""

