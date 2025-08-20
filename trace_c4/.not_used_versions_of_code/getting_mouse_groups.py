#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:28:36 2025

@author: Ilse Klinkhamer
"""

def get_mouse_groups():
    """Returns a dictionary containing categorized mouse groups."""
    return {
        "Switch": [
            "ReserveMouse3", "Dallas", "Flint", "Greene", "Houston", "Iowa", "Jackson",
            "Lincoln", "Newark", "Missouri?", "Pittsburg", "Queens?", "Orleans"
        ],
        "WideExperts": ["Reno", "Seattle", "Yosemite", "Zachary", "Kyiv", "Istanbul", "Copenhagen"],
        "Narrow": ["Rotterdam", "Tallinn", "Quimper", "Porto", "Lisbon", "Madrid"],
        "Bimodal": ["Uppsala", "Venice", "Willemstad", "Zurich", "York", "Xanthi"],
        "Naive": ["Ana1", "Ana2", "Ana3", "Ana4", "Ana5"]
    }

def main():
    mouse_groups = get_mouse_groups()
    
    return mouse_groups
    
    
if __name__ == "__main__":
    main()
    
    
    """
# Example usage:
mouse_groups = get_mouse_groups()
print(mouse_groups["Switch"])  # Prints the list of "Switch" mice
"""