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
import sys
import re
import numpy as np
import shutil


def main(mouse_name=None):
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
        else:
            print("Error: No mouse name provided")
            sys.exit(1)

    #directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
    #folder_name = "SwitchSessionStitching/c4/c4_results_fpfnThreshold_0.1_confidenceRatio_1.5/"

    switch_sessions = True
    classify_again = True

    try:
        print("Converting openephys output to continuous output")
        OpenEphys_wrapper_IK.main(mouse_name, switch_sessions=switch_sessions)
        print("Running c4 analysis")
        run_run_cell_types_classifier.main(mouse_name, classify_again=classify_again, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
        print("Copying c4 output files to Analysis Output")
        copyC4ResultsToAnalysisOutput.main(mouse_name=mouse_name)
    except:
        ik_var = 1
    try:
        print("Inspecting cell types 1.5")
        inspectPredictedCellTypes.main(mouse_name, general_results=False, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
        print("Inspecting cell types 2")
        inspectPredictedCellTypes.main(mouse_name, general_results=False, contamination_ratio=0.1, confidence_ratio_threshold=2)
    except:
        ik_var = 1
    try:
        print("Calculating discharge statistics")
        discharge_statistics.main(mouse_name, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
    except:
        ik_var = 1

if __name__ == "__main__":
    main()