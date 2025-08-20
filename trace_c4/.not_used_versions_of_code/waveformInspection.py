#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:17:24 2025

@author: Ilse Klinkhamer
"""


import numpy as np
import matplotlib.pyplot as plt

def plot_random_waveforms(spike_waveforms, num_samples=10):
    """
    Plots a random selection of spike waveforms from the provided dataset.

    Parameters:
    - spike_waveforms: numpy array of shape (num_neurons, num_timepoints, num_channels)
    - num_samples: Number of waveforms to randomly select and plot
    """
    if spike_waveforms.size == 0:
        print("No waveforms available.")
        return

    # Print original shape
    print(f"Original shape of spike_waveforms: {spike_waveforms.shape}")

    # Use only the first channel
    single_channel_waveforms = spike_waveforms[..., 0]  # Shape: (num_neurons, num_timepoints)

    # Flatten across neurons: (num_neurons, num_timepoints) → (num_neurons * num_timepoints,)
    reshaped_waveforms = single_channel_waveforms.reshape(-1, single_channel_waveforms.shape[-1])

    num_available = reshaped_waveforms.shape[0]
    num_samples = min(num_samples, num_available)  # Avoid oversampling

    # Randomly sample waveforms
    random_indices = np.random.choice(num_available, num_samples, replace=False)
    selected_waveforms = reshaped_waveforms[random_indices]

    # Plot
    plt.figure(figsize=(10, 5))
    for i, waveform in enumerate(selected_waveforms):
        plt.plot(waveform, alpha=0.6, label=f"Waveform {i+1}" if i < 5 else "_nolegend_")  # Limit legend clutter

    plt.xlabel("Time (samples)")
    plt.ylabel("Amplitude")
    plt.title("Randomly Selected Spike Waveforms")
    plt.legend()
    plt.show(block=True)







import numpy as np
import matplotlib.pyplot as plt

def plot_random_waveforms(spike_waveforms, num_samples=10):
    """
    Plots a random selection of spike waveforms from the provided dataset.

    Parameters:
    - spike_waveforms: numpy array of shape (num_neurons * num_timepoints, num_timepoints, num_channels)
    - num_samples: Number of waveforms to randomly select and plot
    """
    if spike_waveforms.size == 0:
        print("No waveforms available.")
        return

    # Concatenate first two axes (neurons and timepoints) → shape (5000, 120, 32)
    concatenated_waveforms = spike_waveforms.reshape(-1, spike_waveforms.shape[2], spike_waveforms.shape[3])
    print(f"New shape after concatenation: {concatenated_waveforms.shape}")  # Should be (5000, 120, 32)

    # Use only the first channel
    single_channel_waveforms = concatenated_waveforms[..., 0]  # Shape: (5000, 120)

    num_available = single_channel_waveforms.shape[0]
    num_samples = min(num_samples, num_available)  # Avoid oversampling

    # Randomly sample waveforms
    random_indices = np.random.choice(num_available, num_samples, replace=False)
    selected_waveforms = single_channel_waveforms[random_indices]

    # Plot
    plt.figure(figsize=(10, 5))
    for i, waveform in enumerate(selected_waveforms):
        plt.plot(waveform, alpha=0.6, label=f"Waveform {i+1}" if i < 5 else "_nolegend_")  # Limit legend clutter

    plt.xlabel("Time (samples)")
    plt.ylabel("Amplitude")
    plt.title("Randomly Selected Spike Waveforms")
    plt.legend()
    plt.show()
    
    
    
    ####################################
    
    import numpy as np
import matplotlib.pyplot as plt

def plot_random_waveforms(spike_waveforms, num_samples=10):
    """
    Plots a random selection of spike waveforms from the provided dataset.

    Parameters:
    - spike_waveforms: numpy array of shape (num_neurons, num_timepoints, num_channels)
    - num_samples: Number of waveforms to randomly select and plot
    """
    if spike_waveforms.size == 0:
        print("No waveforms available.")
        return

    # Print original shape
    print(f"Original shape of spike_waveforms: {spike_waveforms.shape}")

    # Use only the first channel
    single_channel_waveforms = spike_waveforms[..., 0]  # Shape: (num_neurons, num_timepoints)

    # Flatten across neurons: (num_neurons, num_timepoints) → (num_neurons * num_timepoints,)
    reshaped_waveforms = single_channel_waveforms.reshape(-1, single_channel_waveforms.shape[-1])

    num_available = reshaped_waveforms.shape[0]
    num_samples = min(num_samples, num_available)  # Avoid oversampling

    # Randomly sample waveforms
    random_indices = np.random.choice(num_available, num_samples, replace=False)
    selected_waveforms = reshaped_waveforms[random_indices]

    # Plot
    plt.figure(figsize=(10, 5))
    for i, waveform in enumerate(selected_waveforms):
        plt.plot(waveform, alpha=0.6, label=f"Waveform {i+1}" if i < 5 else "_nolegend_")  # Limit legend clutter

    plt.xlabel("Time (samples)")
    plt.ylabel("Amplitude")
    plt.title("Randomly Selected Spike Waveforms")
    plt.legend()
    plt.show(block=True)
    
    
    #######################
    
    
    import numpy as np
import matplotlib.pyplot as plt

def plot_random_waveforms(spike_waveforms, num_samples=100):
    """
    Plots a random selection of spike waveforms from the provided dataset in individual subplots.

    Parameters:
    - spike_waveforms: numpy array of shape (num_neurons, num_timepoints, num_channels)
    - num_samples: Number of waveforms to randomly select and plot
    """
    if spike_waveforms.size == 0:
        print("No waveforms available.")
        return

    # Print original shape
    print(f"Original shape of spike_waveforms: {spike_waveforms.shape}")

    # Use only the first channel
    single_channel_waveforms = spike_waveforms[..., 0]  # Shape: (num_neurons, num_timepoints)

    # Flatten across neurons: (num_neurons, num_timepoints) → (num_neurons * num_timepoints, timepoints)
    reshaped_waveforms = single_channel_waveforms.reshape(-1, single_channel_waveforms.shape[-1])

    num_available = reshaped_waveforms.shape[0]
    num_samples = min(num_samples, num_available)  # Avoid oversampling

    # Randomly sample waveforms
    random_indices = np.random.choice(num_available, num_samples, replace=False)
    selected_waveforms = reshaped_waveforms[random_indices]

    # Set up subplots (10x10 grid)
    fig, axes = plt.subplots(10, 10, figsize=(15, 15))
    fig.suptitle("Randomly Selected Spike Waveforms", fontsize=16)

    for i, ax in enumerate(axes.flatten()):
        if i < num_samples:
            ax.plot(selected_waveforms[i], alpha=0.6)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.axis("off")  # Hide extra subplots

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)  # Adjust title positioning
    plt.show(block=True)

    
    
#########


# Good code!!!

import numpy as np  #IK change: added
import matplotlib.pyplot as plt #IK change: added

def plot_random_waveforms(spike_waveforms, num_samples=100): #IK change: added
    """
    Plots a random selection of spike waveforms from the provided dataset in individual subplots.

    Parameters:
    - spike_waveforms: numpy array of shape (num_neurons, num_timepoints, num_channels)
    - num_samples: Number of waveforms to randomly select and plot
    """
    if spike_waveforms.size == 0:
        print("No waveforms available.")
        return

    # Print original shape
    print(f"Original shape of spike_waveforms: {spike_waveforms.shape}")

    # Use only the first channel
    single_channel_waveforms = spike_waveforms[..., 0]  # Shape: (num_neurons, num_timepoints)

    # Flatten across neurons: (num_neurons, num_timepoints) → (num_neurons * num_timepoints, timepoints)
    reshaped_waveforms = single_channel_waveforms.reshape(-1, single_channel_waveforms.shape[-1])

    num_available = reshaped_waveforms.shape[0]
    num_samples = min(num_samples, num_available)  # Avoid oversampling

    # Randomly sample waveforms
    random_indices = np.random.choice(num_available, num_samples, replace=False)
    selected_waveforms = reshaped_waveforms[random_indices]

    # Set up subplots (10x10 grid)
    fig, axes = plt.subplots(10, 10, figsize=(15, 15))
    fig.suptitle("Randomly Selected Spike Waveforms", fontsize=16)

    for i, ax in enumerate(axes.flatten()):
        if i < num_samples:
            ax.plot(selected_waveforms[i], alpha=0.6)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.axis("off")  # Hide extra subplots

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)  # Adjust title positioning
    plt.show(block=True)