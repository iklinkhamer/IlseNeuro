#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 15:50:45 2025

@author: no1
"""

import os

#source_dir = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/pythonpipeline/data/ephys/kilosort_output/Kourou_Cbx"
#dest_dir = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/Ilse/ephys/Kourou_Cbx"

# Make sure the destination directory exists
os.makedirs(dest_dir, exist_ok=True)

# List all subdirectories in source
for folder_name in os.listdir(source_dir):
    src_path = os.path.join(source_dir, folder_name)
    if os.path.isdir(src_path):
        dest_path = os.path.join(dest_dir, folder_name)
        os.makedirs(dest_path, exist_ok=True)
        print(f"Created: {dest_path}")
