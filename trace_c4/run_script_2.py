#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 14:22:35 2025

@author: Ilse Klinkhamer
"""

import OpenEphys_wrapper_IK
import run_run_cell_types_classifier
import stitch_CH_continuous_files_switch_sessions
import inspectPredictedCellTypes
import discharge_statistics
import copyC4ResultsToAnalysisOutput
from get_dropbox_path import get_dropbox_path
import time
import os

directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
folder_name = "SwitchSessionStitching/c4/c4_results_fpfnThreshold_0.1_confidenceRatio_1.5/"

switch_sessions = True
classify_again = True
mouse_names = ["Ana1", "Ana3"]

for mouse_name in mouse_names:
    try:
        OpenEphys_wrapper_IK.main(mouse_name, switch_sessions=switch_sessions)
        run_run_cell_types_classifier.main(mouse_name, classify_again=classify_again, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
        copyC4ResultsToAnalysisOutput.main(mouse_name=mouse_name)
    except:
        ik_var = 1
    try:
        inspectPredictedCellTypes.main(mouse_name, general_results=False, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
        inspectPredictedCellTypes.main(mouse_name, general_results=False, contamination_ratio=0.1, confidence_ratio_threshold=2)
    except:
        ik_var = 1
    try:
        discharge_statistics.main(mouse_name, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
    except:
        ik_var = 1

