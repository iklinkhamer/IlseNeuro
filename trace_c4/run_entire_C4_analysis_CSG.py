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
import discharge_statistics2 as discharge_statistics
import copyC4ResultsToAnalysisOutput
from getBayesLabDropboxRoot import get_dropbox_path
import time
import os
import sys
import re
import numpy as np
import shutil
from pathlib import Path
import CSG_Env

CSG_MICE = ["Geneva_Cbx", "Georgetown_Cbx", "Georgetown_Ctx", "Helsinki_Cbx", "Ilo_Cbx", "Ilo_mPFC", "Kourou_Cbx", "Limon_Cbx", "Natal_Cbx", "Ocana_Cbx"]
CSG_MICE_CBX = ["Geneva_Cbx", "Georgetown_Cbx", "Helsinki_Cbx", "Ilo_Cbx", "Kourou_Cbx", "Limon_Cbx", "Natal_Cbx", "Ocana_Cbx"]
def main(mouse_name="Ocana_Cbx"
         , kilosort_output_directory = None #os.path.join(get_dropbox_path(), "ContextMouseExperiments/Ilse/ephys")
         , continuous_data_directory = None #os.path.join(get_dropbox_path(), "ContextMouseExperiments/Ilse/ephys")
         , source_folder_name = "kilosort"
         ):
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
        else:
            print("Error: No mouse name provided")
            sys.exit(1)
    if mouse_name in CSG_MICE: 
        if not kilosort_output_directory:
            kilosort_output_directory = CSG_Env.KILOSORT_OUTPUT
        if not continuous_data_directory:
            continuous_data_directory = CSG_Env.EPHYS

    #directory=os.path.join(get_dropbox_path(),"ExperimentOutput/Ephys4Trace1/MainFolder/")
    #folder_name = "SwitchSessionStitching/c4/c4_results_fpfnThreshold_0.1_confidenceRatio_1.5/"
    #continuous_data_mouse_directory = find_path(mouse_name, continuous_data_dir=continuous_data_directory)
    switch_sessions = True
    classify_again = True
    #OpenEphys_wrapper_IK.main(mouse_name, switch_sessions=switch_sessions, directory=kilosort_output_directory)
    run_run_cell_types_classifier.main(mouse_name, classify_again=classify_again, switch_sessions=switch_sessions,
                                       contamination_ratio=0.1, confidence_ratio_threshold=1.5, directory=kilosort_output_directory,
                                       dat_dir=continuous_data_directory, session_folder_pattern=mouse_name[0])

    try:
        print("Converting openephys output to continuous output")
        #OpenEphys_wrapper_IK.main(mouse_name, switch_sessions=switch_sessions, directory=kilosort_output_directory, source_folder_name=source_folder_name)
        print("Running c4 analysis")
        run_run_cell_types_classifier.main(mouse_name, classify_again=classify_again, switch_sessions=switch_sessions, contamination_ratio=0.1, confidence_ratio_threshold=1.5, directory = kilosort_output_directory, dat_dir = continuous_data_directory)
        print("Copying c4 output files to Analysis Output")
        #copyC4ResultsToAnalysisOutput.main(mouse_name=mouse_name)
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
        
def find_path(mouse_name, continuous_data_dir=None):
    #participants_to_exclude = ["Geneva", "Helsinki", "Georgetown", "Ilo", "Ilo_Cbx", "Ilo_mPFC", "Kourou", "Kourou_Cbx", "Limon", "Ocana", "testopto", "test_limon"]
    if not continuous_data_dir:
        continuous_data_dir = CSG_Env.EPHYS #Path(__file__).parent.parent.parent / 'data' / 'ephys'
    continuous_path = continuous_data_dir / 'continuous_data'
    kilosort_path = continuous_data_dir / 'kilosort_output'

    for participant_folder in continuous_path.iterdir():
        if str(participant_folder.name) is not mouse_name:
            continue
        for session_folder in participant_folder.iterdir():
            participant_name = participant_folder.name
            session = session_folder.name
            output_path = kilosort_path / participant_name / session            
            output_path = output_path / 'kilosort_input.bin'
            for item in session_folder.iterdir():
                if item.is_dir():
                    for item2 in item.iterdir():
                        if item2.is_dir():
                            for item3 in item2.iterdir():
                                if item3.isdir():
                                    for item4 in item3.iterdir():
                                        if item4.isdir():
                                            for item5 in item4.iterdir():
                                                if str(item2).endswith("continuous.dat"):
                                                    continuous_data_type = "binary"
                                                    input_path = item
                        elif str(item2).endswith("_ADC1.continuous"):
                            last_underscore_index = str(item2.name).rfind('_')
                            source = str(item2.name)[:last_underscore_index]
                            continuous_data_type = "openephys"
                            input_path = item
                else:
                    if str(item).endswith("_ADC1.continuous"):
                        last_underscore_index = str(item.name).rfind('_')
                        source = str(item.name)[:last_underscore_index]
                        continuous_data_type = "openephys"
                        input_path = item.parent
            if continuous_data_type == "openephys":
                a = 1
            if continuous_data_type == "binary":
                print(input_path)
                a = 1

if __name__ == "__main__":
    main()