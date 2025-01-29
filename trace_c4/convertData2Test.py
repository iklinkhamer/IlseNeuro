#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:13:46 2025

@author: Ilse Klinkhamer
"""

from open_ephys.analysis import Session
folder = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Quimper/Copied folders. Trash?/Quimper_20230727115928 (Copy)/Record Node 105/"
session = Session(folder)
data = session.recordnodes[0].recordings[0].continuous[0].samples
data = data.astype('int16') # convert to 16-bit integers before export
data.tofile('continuous_test.dat') # write samples to a .dat file

#import spikeinterface.extractors as se

#folder = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Quimper/Copied folders. Trash?/Quimper_20230727115928 (Copy)/Data/"
#rec = se.OpenEphysLegacyRecordingExtractor(folder, stream_id="CH")

#data = rec.astype('int16') # convert to 16-bit integers before export
#data.tofile('continuous_test.dat') # write samples to a .dat file