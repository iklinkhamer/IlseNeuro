#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:22:35 2025

@author: Ilse Klinkhamer
"""

import OpenEphys_wrapper_IK
import run_run_cell_types_classifier
import stitch_CH_continuous_files_switch_sessions

runSwitchSessionsPipeline=False
switchSessions = False
classify_again = False
mouse_name = "Amsterdam"


if runSwitchSessionsPipeline:
    #Switch sessions pipeline
    
    stitch_CH_continuous_files_switch_sessions.main(mouse_name)
    
    OpenEphys_wrapper_IK.main(  mouse_name
                              , switchSessions=switchSessions)
    
    run_run_cell_types_classifier.main(     mouse_name
                                       ,    classify_again=classify_again
                                       ,    switchSessions=switchSessions)
    
else:
    OpenEphys_wrapper_IK.main(mouse_name, switchSessions=switchSessions)
    run_run_cell_types_classifier.main(     mouse_name
                                       ,    classify_again=classify_again
                                       ,    switchSessions=switchSessions)

