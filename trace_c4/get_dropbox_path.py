#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 10:23:55 2025

@author: Ilse Klinkhamer
"""


import os

import socket


def get_dropbox_path():
    HOSTNAME = socket.gethostname()
    
    DROPBOX_ENV = {
        "sphinx": "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/",
        "hydra": "/home/devika/BayesLab Dropbox/Julius Koppen/TraceExperiments/"}
    
    DROPBOX_PATH = DROPBOX_ENV[HOSTNAME]
    
    return DROPBOX_PATH
"""
def extract_name_from_path(path):
    #Extract the Dropbox user name from the given path.
    parts = path.split(os.sep)  # Split by folder separator (/)
    
    if "BayesLab Dropbox" in parts and "TraceExperiments" in parts:
        bayeslab_index = parts.index("BayesLab Dropbox")
        traceexp_index = parts.index("TraceExperiments")
        
        if bayeslab_index + 1 == traceexp_index - 1:  # Ensure there's exactly one name in between
            return parts[bayeslab_index + 1]

    return None  # If the pattern isn't found

def get_dropbox_path():
    #Constructs and returns the corrected Dropbox path.
    relative_path = "Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/"
    
    # Construct the full path dynamically
    directory = os.path.join(os.path.expanduser("~"), relative_path)
    
    dropbox_user = extract_name_from_path(directory)
    
    if dropbox_user:
        # Replace first occurrence of "Lucas Bayones" with the extracted user
        dropbox_path = directory.replace("Lucas Bayones", dropbox_user, 1)
        
        if os.path.exists(dropbox_path):
            return dropbox_path
        else:
            # Remove the first instance of "Lucas Bayones" before "BayesLab Dropbox"
            parts = directory.split(os.sep)
            if dropbox_user in parts:
                parts.remove(dropbox_user)  # Remove only the first occurrence
                new_path = os.sep.join(parts)

                if os.path.exists(new_path):
                    return new_path

    return None  # Return None if the Dropbox user or path cannot be determined
"""
# When used as a script, print the path
if __name__ == "__main__":
    path = get_dropbox_path()
    if path:
        print("Updated Dropbox Path:", path)
    else:
        print("Could not determine Dropbox path.")
