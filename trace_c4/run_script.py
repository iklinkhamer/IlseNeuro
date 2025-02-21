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

runSwitchSessionsPipeline=False
switch_sessions = False
classify_again = True
mouse_names = ["Seattle"]

for mouse_name in mouse_names:
    if runSwitchSessionsPipeline:
        #Switch sessions pipeline

        stitch_CH_continuous_files_switch_sessions.main(mouse_name)

        OpenEphys_wrapper_IK.main(  mouse_name
                                  , switch_sessions=switch_sessions)

        run_run_cell_types_classifier.main(     mouse_name
                                           ,    classify_again=classify_again
                                           ,    switch_sessions=switch_sessions)

    else:
        #OpenEphys_wrapper_IK.main(mouse_name, switchSessions=switchSessions)
        run_run_cell_types_classifier.main(mouse_name, classify_again=classify_again, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
        inspectPredictedCellTypes.main(mouse_name, general_results=False)
