#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 16:57:45 2025

@author: Ilse Klinkhamer
"""

import OpenEphys_wrapper_IK
import time
from get_dropbox_path import get_dropbox_path
import os

mouse_names = ["Georgetown_Cbx"]
directory = os.path.join(os.path.dirname(os.path.dirname(get_dropbox_path())), "ContextMouseExperiments/Ilse/ephys")

"""
for m in range(len(mouse_names)):
    mouse_name = mouse_names[m]

    while True:
        try:
            
            OpenEphys_wrapper_IK.main(mouse_name)

            if True:  
                print("Code ran successfully!")
                break  # Exit the loop when successful
        except Exception as e:
            print(f"Error occurred: {e}. Retrying in 1 hour...")
            time.sleep(3600)  # Wait for an hour before retrying
   """ 
    
for m in range(len(mouse_names)):
    mouse_name = mouse_names[m]

    
    OpenEphys_wrapper_IK.main(mouse_name, directory=directory, source_folder_name = "kilosort", session_folder_pattern=mouse_name[0])

           

    