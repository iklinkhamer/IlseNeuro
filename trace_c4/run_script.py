#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:22:35 2025

@author: Ilse Klinkhamer
"""

import OpenEphys_wrapper_IK
import run_run_cell_types_classifier
import stitch_CH_continuous_files_switch_sessions

runSwitchSessionsPipeline=True
mouse_name = "Ana2"


if runSwitchSessionsPipeline:
    #Switch sessions pipeline
    
    stitch_CH_continuous_files_switch_sessions(mouse_name)
    
    OpenEphys_wrapper_IK.main(  mouse_name
                              , switchSessions=True)
    
    run_run_cell_types_classifier.main(     mouse_name
                                       ,    classify_again=False
                                       ,    switchSessions=True)
    
else:
    #OpenEphys_wrapper_IK.main(mouse_name, switchSessions=False)
    run_run_cell_types_classifier.main(     mouse_name
                                       ,    classify_again=False
                                       ,    switchSessions=False)

