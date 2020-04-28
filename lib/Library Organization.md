


## ChrTracer2
This is a graphical user interface to guide the reader through data analysis as described in our current manuscript on Optical Reconstruction of Chromatin Architecture (submission pending). 

#### Input
It requires a data folder which contains the raw images and an excel table which records how the data files are named and organized in subfolders. 

#### Output
* ChrTracer2 converts raw image data into vector data of barcode positions, and saves these as files as comma separated value (csv) text files that record the positions and associated statistics for all molecules localized.  
* Additionally, images of the field of view, and 4D renderings of the raw data for each chromosome are saved in the target folder.  
* For speed and data backup, the current version also saves new copies of all raw data with x-y drift removed.

#### Procedure
The first step loads the data and computes an x-y alignment to pixel accuracy using a single frame from each z-stack.  The user may select a variety of contrast parameters to facilitate accurate correlation based alignment of the data, and change the size of field of view used for the alignment to increase the speed or decrease the probability of failed registration (for sparse data).  

## GeneralDependencies 
The directory includes commonly used functions, including those integral to managing variable input to our custom functions.

## ImageProcessing
Some custom routines for general image processing required by ChrTracer2 and broadly useful in other routines as well are included here. 