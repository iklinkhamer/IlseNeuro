#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 19:28:01 2025

@author: Ilse Klinkhamer
"""

import run_entire_C4_analysis_CSG

CSG_MICE_CBX = [
    "Geneva_Cbx", "Georgetown_Cbx", "Helsinki_Cbx",
    "Ilo_Cbx", "Kourou_Cbx", "Limon_Cbx", "Natal_Cbx", "Ocana_Cbx"
]

for mouse_name in CSG_MICE_CBX:
    try:
        run_entire_C4_analysis_CSG.main(mouse_name)
    except Exception as e:
        print(f"Error processing {mouse_name}: {e}")
