#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 12:14:13 2025

@author: Ilse Klinkhamer
"""
import numpy as np

def c4_colors_rgb():
    C4_COLORS = {
        "PkC_ss": [28, 120, 181],
        "PkC_cs": [0, 0, 0],
        "MLI": [224, 85, 159],
        "MFB": [214, 37, 41],
        "GrC": [143, 103, 169],
        "GoC": [56, 174, 62],
        "laser": [96, 201, 223],
        "drug": [239, 126, 34],
        "background": [244, 242, 241],
        "MLI_A": [224, 85, 150],
        "MLI_B": [220, 80, 160],
    }
    return C4_COLORS

def c4_colors_rgb_normalized():
    C4_COLORS_RGB = c4_colors_rgb()  # Assuming this returns a dict of RGB values
    C4_COLORS_RGB_N = {key: [v / 255 for v in value] for key, value in C4_COLORS_RGB.items()}
    return C4_COLORS_RGB_N


def sRgb2linRgb(sRgb):
    """
    Converts sRGB to linear RGB.
    
    Parameters:
    - sRgb: A list, tuple, or NumPy array of 3 values (0-1 or 0-255).
    
    Returns:
    - Linear RGB values as a NumPy array.
    """
    sRgb = np.array(sRgb, dtype=float) / 255  # Normalize to 0-1 if needed
    
    def decode(s):
        return s / 12.92 if s <= 0.04045 else ((s + 0.055) / 1.055) ** 2.4

    return np.array([decode(s) for s in sRgb])

def linRgb2sRgb(cieRgb):
    """
    Converts linear RGB to sRGB.
    
    Parameters:
    - cieRgb: A NumPy array of linear RGB values.
    
    Returns:
    - sRGB values (0-255 range).
    """
    def encode(cie):
        return 12.92 * cie if cie <= 0.0031308 else 1.055 * (cie ** (1 / 2.4)) - 0.055

    return np.clip([encode(cie) for cie in cieRgb], 0, 1) * 255

def mkGradientFn(sRgbFrom, sRgbTo):
    """
    Creates a function that interpolates between two colors in linear RGB space.
    
    Parameters:
    - sRgbFrom: Starting sRGB color.
    - sRgbTo: Target sRGB color.
    
    Returns:
    - A function that interpolates between the two colors based on input x (0 to 1).
    """
    linFrom = sRgb2linRgb(sRgbFrom)
    linTo = sRgb2linRgb(sRgbTo)

    def gradientFn(x):
        return linRgb2sRgb(linFrom + x * (linTo - linFrom))
    
    return gradientFn

def lighten(sRgb, percentage=10):
    """
    Lightens an sRGB color by a given percentage towards white.
    
    Parameters:
    - sRgb: A list or tuple of 3 numerical values (RGB) in range [0, 255].
    - percentage: Percentage to lighten (default is 10).
    
    Returns:
    - Lightened RGB color (0-255 range).
    """
    gradientFn = mkGradientFn(sRgb, [255, 255, 255])
    return gradientFn(percentage / 100).astype(int)

# Example usage
lightened_color = lighten([128, 64, 32], 20)
#print(lightened_color)  # Output: [179 115  83]

def normalize_RGB_dict(color_dict):
    normalized_color_dict = {key: [v / 255 for v in value] for key, value in color_dict.items()}
    return normalized_color_dict

    
