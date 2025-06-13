#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 10:23:55 2024

@author: Ilse Klinkhamer
"""

"""
#from npyx.testing import test_npyx
#sys.path.insert(0, '/home/no1/anaconda3/envs/new_env/lib/python3.10/site-packages')  # Ensure this is the correct path to 'new_env'
#sys.path.insert(0, '/home/no1/anaconda3/envs/new_env/lib/python3.10/site-packages/npyx')  # Add npyx explicitly if needed
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf' # Path to the Arial font file (adjust the path based on your system)
font_manager.fontManager.addfont(font_path) # Add the font to Matplotlib's font manager
"""

import sys
import npyx

from get_dropbox_path_2 import get_dropbox_path
from npyx.c4.predict_cell_types import run_cell_types_classifier
from os import path
import os

def run_cell_types_classifier_wrapper(mouse_name
                                      ,classify_again=True
                                      ,switch_sessions=False
                                      ,contamination_ratio=0.1
                                      ,confidence_ratio_threshold=2
                                      ,directory=os.path.join(get_dropbox_path(),"TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/")
                                      ,cache_dir=os.path.join(get_dropbox_path(), "C4 Cache")
                                      , dat_dir=os.path.join(get_dropbox_path(), "TraceExperiments/C4_conversion")
                                      , session_folder_pattern = ""
                                      ,skip_without_continuous=True):

    if not session_folder_pattern:
        session_folder_pattern = mouse_name

    # any spike sorted recording compatible with phy
    # (e.g. kilosort output)
    dp_base = os.path.join(directory, mouse_name)
    if "ReserveMouse" in mouse_name:
        dp_base = dp_base.replace("MainFolder", "ReserveFolder")

    # Get all folders with mouse_name in their name and without "copy"
    mouse_folders = [
        folder for folder in os.listdir(dp_base)
        if os.path.isdir(os.path.join(dp_base, folder))
           and session_folder_pattern in folder
           and "copy" not in folder
           and "Copy" not in folder
    ]
    mouse_folders.sort()

    if switch_sessions:
        switch_folder = os.path.join(dp_base, "SwitchSessionStitching")
        if os.path.exists(switch_folder):
            switch_folder_name = "SwitchSessionStitching"
            mouse_folders.append(switch_folder_name)

    phy_folder = "c4"

    for sess_oebin in mouse_folders:

        print(f"Processing folder: {sess_oebin}")

        #if sess_oebin != "Ana3_20190531193955":
        #   continue

        dp = path.join(dp_base, sess_oebin, phy_folder)
        cla_res_path = path.join(dp_base, sess_oebin, phy_folder, "cell_type_classification")
        save_path = os.path.join(dp, f"c4_results_fpfnThreshold_{contamination_ratio}_confidenceRatio_{confidence_ratio_threshold}")
        results_file_path = os.path.join(save_path, "cluster_predicted_cell_type.tsv")

        if not os.path.exists(dp):
            print(f"Folder {dp} does not exist, skipping c4 analysis.")
            continue  # Skip to the next iteration if the folder doesn't exist
        if not os.path.exists(f"{dp}/params.py"):
            print(f"params.py folder not found in folder {sess_oebin}, skipping c4 analysis...")
            continue
        if not classify_again and os.path.isfile(results_file_path):
            print(f"Session {dp} has already been classified and classify again is false, skipping c4 analysis.")
            continue
        os.makedirs(save_path, exist_ok=True)

        cache_path = os.path.join(cache_dir, mouse_name, sess_oebin, ".NeuroPyxels")
        os.makedirs(cache_path, exist_ok=True)

        dat_path = os.path.join(dat_dir, mouse_name, sess_oebin)
        if not os.path.exists(os.path.join(dat_path, "continuous")):
            dat_path = os.path.join(dat_path, phy_folder)
            if not os.path.exists(os.path.join(dat_path, "continuous")):
                print("continuous folder not found in data path folder, so please check.")


        run_cell_types_classifier(dp, quality = 'all', parallel = False, fp_threshold = contamination_ratio, fn_threshold = contamination_ratio, threshold = confidence_ratio_threshold, save_path = save_path, cache_path=cache_path, dat_path=dat_path)

        # if any test fails, re-run them with the following to print the error log, and try to fix it or post an issue on github:
        #run_cell_types_classifier(dp, raise_error=True)


def main(mouse_name="ReserveMouse3", classify_again=True, switch_sessions=True, contamination_ratio=0.1, confidence_ratio_threshold=1.5
         , directory=os.path.join(get_dropbox_path(),"TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/")
         , cache_dir=os.path.join(os.path.dirname(get_dropbox_path().rstrip("/")), "C4 Cache")
         , dat_dir=os.path.join(get_dropbox_path(), "C4_conversion")
         , session_folder_pattern = ""
         ):
    if mouse_name is None:
        if len(sys.argv) > 1:
            mouse_name = sys.argv[1]
        else:
            print("Error: No mouse name provided")
            sys.exit(1)
    run_cell_types_classifier_wrapper(mouse_name, classify_again, switch_sessions, contamination_ratio, confidence_ratio_threshold, directory=directory, cache_dir=cache_dir, dat_dir=dat_dir, session_folder_pattern = session_folder_pattern)


if __name__ == "__main__":
    main()


