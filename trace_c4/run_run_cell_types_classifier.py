#from npyx.testing import test_npyx
import sys
sys.path.insert(0, '/home/no1/anaconda3/envs/new_env/lib/python3.10/site-packages')  # Ensure this is the correct path to 'new_env'
sys.path.insert(0, '/home/no1/anaconda3/envs/new_env/lib/python3.10/site-packages/npyx')  # Add npyx explicitly if needed

from matplotlib import font_manager
font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf' # Path to the Arial font file (adjust the path based on your system)
font_manager.fontManager.addfont(font_path) # Add the font to Matplotlib's font manager

from predict_cell_types_IK import run_cell_types_classifier
from os import path
import os

mouse_name = "Yosemite"
classify_again = True

# any spike sorted recording compatible with phy
# (e.g. kilosort output)
dp_base = f"/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/{mouse_name}/"

# Get all folders with mouse_name in their name and without "copy"
mouse_folders = [
    folder for folder in os.listdir(dp_base)
    if os.path.isdir(os.path.join(dp_base, folder))
    and mouse_name in folder
    and "copy" not in folder
    and "Copy"  not in folder
]
mouse_folders.sort()

phy_folder = "c4"

for sess_oebin in mouse_folders:
    
    print(f"Processing folder: {sess_oebin}")
    
    dp = path.join(dp_base, sess_oebin, phy_folder)
    cla_res_path = path.join(dp_base, sess_oebin, phy_folder, "cell_type_classification")
    if not os.path.exists(dp):
        print(f"Folder {dp} does not exist, skipping.")
        continue  # Skip to the next iteration if the folder doesn't exist
    if not classify_again and os.path.exists(cla_res_path):
        print(f"Session {dp} has already been classified and classify again is false, skipping.")
        continue

    run_cell_types_classifier(dp, quality = 'all', parallel = False, fp_threshold = 0.1, fn_threshold = 0.1)

    # if any test fails, re-run them with the following to print the error log, and try to fix it or post an issue on github:
    #run_cell_types_classifier(dp, raise_error=True)