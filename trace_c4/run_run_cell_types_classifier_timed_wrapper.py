#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 17:01:58 2025

@author: Ilse Klinkhamer
"""

import time
import run_run_cell_types_classifier


mouse_names = ["ReserveMouse3", "Ana2", "Ana4", "Ana5"]

for m in range(len(mouse_names)):
    mouse_name = mouse_names[m]

    while True:
        try:
            # Your main code here
            run_run_cell_types_classifier.main(mouse_name)
            # Simulating successful execution (remove this in actual use)
            if True:  
                print("Code ran successfully!")
                break  # Exit the loop when successful
        except Exception as e:
            print(f"Error occurred: {e}. Retrying in 1 hour...")
            time.sleep(3600)  # Wait for an hour before retrying
    
    
