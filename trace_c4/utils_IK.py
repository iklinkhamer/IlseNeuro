#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:00:44 2025

@author: Ilse Klinkhamer
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import json

def save_json(data, filename):
    """Save a dictionary as a JSON file."""
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)

def load_json(filename):
    """Load a JSON file into a dictionary."""
    with open(filename, 'r') as f:
        return json.load(f)

def ensure_dir(directory):
    """Create a directory if it does not exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def plot_signal(signal, title='Signal', xlabel='Time', ylabel='Amplitude'):
    """Plot a given signal using Matplotlib."""
    plt.figure(figsize=(10, 4))
    plt.plot(signal)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.show()

def normalize_array(arr):
    """Normalize a NumPy array to the range [0, 1]."""
    return (arr - np.min(arr)) / (np.max(arr) - np.min(arr))

def moving_average(data, window_size):
    """Compute the moving average of a 1D array."""
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

