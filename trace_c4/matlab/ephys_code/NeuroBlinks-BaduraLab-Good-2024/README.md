Neuroblinks 
===========

A program for doing eyeblink classical conditioning using Matlab and inexpensive microcontrollers. 

Author: Shane Heiney, University of Iowa

This version is a fork of the code originally written by Shane Heiney, Shogo Ohmae, and Olivia Kim in the Medina lab at Baylor College of Medicine

Copyright &copy; 2018


Requirements
------------

* Modern Windows computer (tested on Win 10) with lots of RAM and multi-core processor (e.g. 16-32 GB RAM, Intel i7).
* Matlab v 2015 or later with the Image Acquisition and Image Processing toolboxes. Note that the code used to control the camera varies based on the particular release of Matlab because of changes in the Image Acquisition Toolbox. 
* A high frame rate video camera (at least 200 FPS) supported by the Image Acquisition toolbox, for instance the Basler [Ace ac1300](https://www.baslerweb.com/en/products/cameras/area-scan-cameras/ace/aca1300-200um/). If you use a different camera you may need to modify the camera control code in Neuroblinks. We have only tested the Ace with this version. The Basler Pulse is also used in this version for pupil measurements.
* [Teensy microcontroller](https://www.pjrc.com/store/teensy35.html). Should also work with 32-bit Arduino microcontrollers such as Due and Zero. 


Instructions
------------

1. Configure your camera(s) to work with Matlab using the GenTL adapter. 

2. Configure your hardware (Teensy or other microcontroller). The necessary C++ source code files are included in the "private" directory under "neuroblinks". You will also need to configure the inputs and outputs. The code is written to work with the Teensy 3.5 but should be relatively easy to extend to other devices. Contact us if you have any questions about how to do this. 

3. Download or clone the "neuroblinks" project to your computer and add the main directory of the project (i.e., "neuroblinks") to your Matlab path. Do not include subfolders as these will be added automatically for you. Manually adding subfolders below "neuroblinks" could result in unexpected behavior. 

4. Modify the three configuration files (`neuroblinks_config.m`, `rig_config.m`, and `user_config.m`) to match your particular configuration (see comments within the m-files).

5. Create a directory for your mouse, e.g. `<data root>\<animalID>` and optionally add a copy of `condparams.csv` to this directory (an example is included in the Neuroblinks distribution). For example, we use `data\MXXX`, where XXX is the unique mouse ID. Note that if you follow our naming conventions, the metadata.mouse field will be automatically populated for you and the session directory will be automatically created as "YYMMDD", which is the current date (`datestr(now,'yymmdd')`. If you deviate from this naming convention and directory structure you will likely need to modify the code and some of our analysis code might not work for you. 

6. Run `neuroblinks()` to start the main program (with optional arguments).

This [document](https://docs.google.com/document/d/1InIuTQ_H1JthY9_0v9BHR0_naf4_ief3BlzZlXbdtBc/edit?usp=sharing)  provides some brief instructions on using the Neuroblinks GUI to help get you started. 

You'll want to set up the camera and IR light source so that the mouse's face is in frame and the entire eye is in focus. The IR light should be positioned to give good contrast between the iris/pupil and surrounding fur. See [Heiney et al, 2014b](http://www.ncbi.nlm.nih.gov/pubmed/25378152) for more information. 

The data for each trial, including a raw video of the trial and the metadata, will be stored in separate .MAT files within the session directory using the base name that you specify plus an auto-incremented suffix corresponding to the trial number. Contact us if you would like guidance on how to process these files to get calibrated eyelid traces and metadata in a format that's easier to analyze. 