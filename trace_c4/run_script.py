#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:22:35 2025

@author: Ilse Klinkhamer
"""

import OpenEphys_wrapper_IK
import run_run_cell_types_classifier
import stitch_CH_continuous_files_switch_sessions
import inspectPredictedCellTypes
import discharge_statistics
import time
import os

directory = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/"
folder_name = "SwitchSessionStitching/c4/c4_results_fpfnThreshold_0.1_confidenceRatio_1.5/"

runSwitchSessionsPipeline=False
switch_sessions = True
classify_again = False
mouse_names = ["Zurich"]

for mouse_name in mouse_names:
    if runSwitchSessionsPipeline:
        #Switch sessions pipeline
        stop = False
        while not stop:
            try:
                stitch_CH_continuous_files_switch_sessions.main(mouse_name)
            except:
                a = 1

            try:

                OpenEphys_wrapper_IK.main(mouse_name, switch_sessions=switch_sessions)
            except:
                a = 1
            try:
                run_run_cell_types_classifier.main(     mouse_name
                                               ,    classify_again=classify_again
                                               ,    switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
            except:
                a = 1

            if os.path.exists(os.path.join(directory, mouse_name, folder_name, "cluster_predicted_cell_type.tsv")):
                print(f"Switch sessions successfully analyzed for {mouse_name}")
                stop=True
            else:
                print(f"Failed analyzing switch sessions for {mouse_name}, waiting 60 seconds and trying again after that...")
                time.sleep(60)
    else:
        #OpenEphys_wrapper_IK.main(mouse_name, switch_sessions=switch_sessions)
        run_run_cell_types_classifier.main(mouse_name, classify_again=classify_again, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
        inspectPredictedCellTypes.main(mouse_name, general_results=False)
        discharge_statistics.main(mouse_name, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
