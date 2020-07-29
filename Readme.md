# Analysis Routines for Optical Reconstruction of Chromatin Architecture 

## Welcome
This repository contains code written primarily in Matlab for the analysis of microscopy data for optical reconstruction of chromatin architecture as described in our manuscript on Optical Reconstruction of Chromatin Architecture (Mateo et al 2019, Nature) and the updated version of this code provided with our protocol paper (in review).  

Additional code and tutorials on use will be added in the future to facilitate adoption, so please check back!  Meanwhile all analysis scripts used for our current work may be found here. 

## Licensed CC BY, Boettiger Lab 2018, 2020
You are free to use, modify and distribute all elements of this repository provided the original authors are attributed, as described by the Creative Commons license. 

## Installation
1. Install Matlab   
    - This software is written in Matlab and requires and existing installation of Matlab (TM) R2016 or later.
2. Pull the MERFISH_analysis repo: https://github.com/ZhuangLab/MERFISH_analysis. 
    - Some of the probe construction files use matlab classes from the MERFISH probe construction software written by Jeff Moffitt, which must be downloaded from the following repository due to licensing restrictions described in the readme file included with that repository. 
3. Add the MERFISH_analysis repo to your matlab filepath
    - e.g. `addpath(genpath('C:\code\MERFISH_analysis\'))`
    - Tip, if you add this command to your Matlab startup script, it will execute automatically when you start matlab.
4. Add this folder to your matlab filepath 
    - e.g. `addpath(genpath('C:\code\ORCA-public\'))`

### Note on data format
Currently the software assumes the data are saved in ".dax" format - a flat binary form with header information specified in an accompanying ".inf" file with the same name. This is the output of our microscope software [storm-control](https://github.com/Boettiger-lab/storm-control).  Please note, this is a fork of the popular "storm-control" branch developed by Hazen Babcock and colleagues in the Zhuang lab at Harvard University: [storm-control](https://github.com/ZhuangLab/storm-control).  

### Dependencies
* Matlab (TM) R2016 or later
* The code has been primarily used on Windows. Minor adjustments of some filesep symbols are sufficient to migrate to Linux or Mac.
* It is recommended to have >16 Gb of RAM for running analysis.

## Demo Data
* Demo data can be found here: https://bit.ly/2S6eCjk.  Warning: while dramatically reduced from the full dataset, this file is still large (~350 Gb). We aim to provide sufficient data to explore the rich features of the highly multiplexed ORCA data and illustrate its analysis.

### Walk-throughs
* Detailed explanations of the different Graphical User Interface Tools in this repository are provided in the protocols paper, please stay tuned for links to the text and figures. 

## Questions
 Please use the github site for this project to log bug reports or questions.  You may also contact Prof. Boettiger directly by sending email to boettiger *at* stanford dot edu. 
