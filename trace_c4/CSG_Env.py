#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 16:27:48 2025

@author: Ilse Klinkhamer
"""

from pathlib import Path
import inspect
from getBayesLabDropboxRoot import get_dropbox_path

DATA_ROOT = Path(get_dropbox_path(), "ContextMouseExperiments", "pythonpipeline", "data")

ANALYSIS_OUTPUT = Path(DATA_ROOT, "AnalysisOutput")
EXPERIMENT_OUTPUT = Path(DATA_ROOT, "ExperimentOutput")

BEHAVIOUR = Path(EXPERIMENT_OUTPUT, "behaviour")
BEHAVIOUR_TRAINING = Path(BEHAVIOUR, "training")
BEHAVIOUR_RECORDING = Path(BEHAVIOUR, "recording")

EPHYS = Path(EXPERIMENT_OUTPUT, "ephys", "continuous_data")

EYEBLINK = Path(ANALYSIS_OUTPUT, "eyeblink")
EYEBLINK_TRAINING = Path(EYEBLINK, "training")
EYEBLINK_RECORDING = Path(EYEBLINK, "recording")

KILOSORT_OUTPUT = Path(ANALYSIS_OUTPUT, "kilosort_output")

C4 = Path(ANALYSIS_OUTPUT, "c4")
CELL_COUNTS = Path(ANALYSIS_OUTPUT, "cell_type_counts")
FIRING_RATES = Path(ANALYSIS_OUTPUT, "smoothed_firing_rates")
PCA = Path(ANALYSIS_OUTPUT, "PCA")
SCALING_INDEX = Path(ANALYSIS_OUTPUT, "scaling_index")

RASTERS = Path(ANALYSIS_OUTPUT, "PSTH_Rasters")

RASTERS_FULL = Path(RASTERS, "figures_full")
RASTERS_FULL_NO_SCALING = Path(RASTERS_FULL, "no_scaling")
RASTERS_FULL_AMP = Path(RASTERS_FULL, "amplitude_scaling")
RASTERS_FULL_TEMP = Path(RASTERS_FULL, "temporal_scaling")

RASTERS_SIG_R2 = Path(RASTERS, "figures_sig_R2")
RASTERS_SIG_R2_NO_SCALING = Path(RASTERS_SIG_R2, "no_scaling")
RASTERS_SIG_R2_AMP = Path(RASTERS_SIG_R2, "amplitude_scaling")
RASTERS_SIG_R2_TEMP = Path(RASTERS_SIG_R2, "temporal_scaling")

RASTERS_OPTO = Path(RASTERS, "figures_opto")

# Individual subfolder names
TRIALPARAMETERS = "TrialParameters"
CAMERA = "CameraData"
WHEEL = "WheelData"
OPTO = "OptoData"

EXPERIMENT1 = Path("experiment1", "recording1", "continuous", "Acquisition_Board-108.acquisition_board")


def ensure_dirs():
    """Prompt for each missing Path object defined in this module whether to create it.
    Only directories inside DATA_ROOT can be created.
    Options: y/yes (create), n/no (skip), a/all (create all remaining without asking).
    """
    current_module = inspect.getmodule(ensure_dirs)
    create_all = False

    for name, value in vars(current_module).items():
        if isinstance(value, Path) and not value.exists():
            # Skip silently if outside DATA_ROOT
            try:
                value.relative_to(DATA_ROOT)
            except ValueError:
                continue  # silently skip

            if not create_all:
                answer = input(f"Directory {value} does not exist. Create it? [y/N/a]: ").strip().lower()
                if answer in ["a", "all"]:
                    create_all = True
                    answer = "y"

            if create_all or answer in ["y", "yes"]:
                value.mkdir(parents=True, exist_ok=True)
                print(f"Created {value}")
            else:
                print(f"⚠️  Skipped creating {value}")

