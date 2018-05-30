# Analysis Routines for Optical Reconstruction of Chromatin Architecture

## Welcome
This repository contains code written primarily in Matlab for the analysis of microscopy data for optical reconstruction of chromatin architecture as described in our current manuscript on Optical Reconstruction of Chromatin Architecture (submission pending).  

Additional code and tutorials on use will be added in the future to facilitate adoption, so please check back!  Meanwhile all analysis scripts used for our current work may be found here. 

## Licensed CC BY, Boettiger Lab 2018
You are free to use, modify and distribute all elements of this repository provided the original authors are attributed, as described by the Creative Commons license. 

## Installation
Simply add this folder and its children to you matlab filepath, e.g.: `addpath(genpath('C:\code\ORCA-public'))`.  

### Note on data format
Currently the software assumes the data are saved in ".dax" format - a flat binary form with header information specified in an accompanying ".inf" file with the same name. This is the output of our microscope software [storm-control](https://github.com/Boettiger-lab/storm-control).  Please note, this is a fork of the popular "storm-control" branch developed by Hazen Babcock and colleagues in the Zhuang lab at Harvard University: [storm-control](https://github.com/ZhuangLab/storm-control).  

## Questions
This is our first release. Please use the github site for this project to log bug reports or questions.  You may also contact Prof. Boettiger directly by sending email to boettiger *at* stanford dot edu. 
