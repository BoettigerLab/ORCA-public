function designStruct = NameOligoArrayFiles(designStruct, filePath, vargin)
%--------------------------------------------------------------------------
% designStruct = NameOligoArrayFiles(designStruct, filePath)
% This function populates the various file path fields in a design
% structure with a name derived from the various elements
%
%--------------------------------------------------------------------------
% Outputs:
%
% designStruct/struct: A structure containing fields with the required
% parameters for OligoArray2.0
%
%--------------------------------------------------------------------------
% Inputs:
%
%--------------------------------------------------------------------------
% Variable Inputs:
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% November 21, 2013
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Return design struct with basic defaults
%--------------------------------------------------------------------------
designStruct.inputFile = '';
designStruct.databaseFile = '';
designStruct.outputFile = '';
designStruct.failedFile = '';
designStruct.logFile = '';
designStruct.maxOligos = 30;
designStruct.minLength = 32;
designStruct.maxLength = 32;
designStruct.maxDistance = 7000;
designStruct.minTm = 75;
designStruct.maxTm = 90;
designStruct.secStructTm = 70;
designStruct.crossHybTm = 70;
designStruct.minGC = 35;
designStruct.maxGC = 80;
designStruct.mask = '"GGGG;CCCC;TTTT;AAAA"';
designStruct.minSep = 40;
designStruct.numParallel = 4;
