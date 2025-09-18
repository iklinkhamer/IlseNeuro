#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 11:38:10 2025

@author: Ilse Klinkhamer
"""
import os
import socket

def get_dropbox_path():
    HOSTNAME = socket.gethostname()
    
    DROPBOX_ENV = {
        "sphinx": "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/",
        "hydra": "/home/devika/BayesLab Dropbox/Julius Koppen/",
        "DESKTOP-BSHMJ1M": "D:\BayesLab Dropbox\Ilse Klinkhamer",
        "LAPTOP-I93JEBDJ": r"C:\Users\31681\BayesLab Dropbox\Ilse Klinkhamer"}
    
    DROPBOX_PATH = DROPBOX_ENV[HOSTNAME]
    
    return DROPBOX_PATH

# When used as a script, print the path
if __name__ == "__main__":
    path = get_dropbox_path()
    if path:
        print("Updated Dropbox Path:", path)
    else:
        print("Could not determine Dropbox path.")
