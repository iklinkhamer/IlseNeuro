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
switchSessions = False
classify_again = True
mouse_names = ["Yosemite"]

for mouse_name in mouse_names:
    if runSwitchSessionsPipeline:
        #Switch sessions pipeline

        stitch_CH_continuous_files_switch_sessions.main(mouse_name)

        OpenEphys_wrapper_IK.main(  mouse_name
                                  , switchSessions=switchSessions)

        run_run_cell_types_classifier.main(     mouse_name
                                           ,    classify_again=classify_again
                                           ,    switchSessions=switchSessions)

    else:
        #OpenEphys_wrapper_IK.main(mouse_name, switchSessions=switchSessions)
        #run_run_cell_types_classifier.main(mouse_name, classify_again=classify_again, switchSessions=switchSessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5)
        inspectPredictedCellTypes.main(mouse_name, general_results=False)
