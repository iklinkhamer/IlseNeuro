#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 19:28:01 2025

@author: Ilse Klinkhamer
"""

import run_entire_C4_analysis_CSG

CSG_MICE_CBX = [
    "Helsinki_Cbx", "Kourou_Cbx", "Natal_Cbx", "Ocana_Cbx"
]

path_helsinki ="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/ContextMouseExperiments/Ilse/ephys/"

for mouse_name in CSG_MICE_CBX:
    try:
        if mouse_name == "Helsinki_Cbx":
            run_entire_C4_analysis_CSG.main(mouse_name, continuous_data_directory=path_helsinki)
        else:
            run_entire_C4_analysis_CSG.main(mouse_name)
    except Exception as e:
        print(f"Error processing {mouse_name}: {e}")
