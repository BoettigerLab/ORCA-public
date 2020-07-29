% RunPipeline5Demo
% Alistair Boettiger, 05/25/15
% 
% Goal: Re-analyze all previous data sets with new pixel-based pipeline.


clc; startup; 

%% Input parameters
% File Info
MERFISHdrive = '\\Morgan\MorganData2\';
MERFISHdata = [MERFISHdrive,'MERFISHdataArchive\'];
analysisXLS = [MERFISHdata,'MERFISHexperiments.xls'];
dataToAnalyze = LoadDatasets(analysisXLS, 'verbose', true);

% Parameters from other functions
% These fields can all be commented out / deleted to use the defaults built
% into the individual functions. 

% LoadMERFISHdax pars
defaults = cell(0,3);
defaults(end+1,:) = {'imageTag', 'string', 'STORM'}; % Base tag for all images
defaults(end+1,:) = {'numFrames', 'integer', 1}; %  Number of frames in raw dax movie (could be determined automatically).
defaults(end+1,:) = {'bitFrames', 'integer', 1}; %   vector of the frames in the movie containing images of RNA bits  
% shared pars
defaults(end+1,:) = {'overwriteDax', 'boolean', false}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; % 

% SaveCellImagesByBit pars
psf = fspecial('gaussian',10,2);
upsample = 3;
defaults(end+1,:) = {'upsampleWarp', 'positive', 3}; % for correlation based drift alignment 
defaults(end+1,:) = {'minBeadContrast', 'fraction', .5}; % for correlation based drift alignment 
defaults(end+1,:) = {'upsample', 'positive', upsample}; % 
defaults(end+1,:) = {'psf', 'array', psf}; % 
defaults(end+1,:) = {'numIters', 'positive', 5}; % 
defaults(end+1,:) = {'corrAlignImages', 'boolean', true}; % 
defaults(end+1,:) = {'showImages', 'boolean', false}; % 
% shared pars
defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
% defaults(end+1,:) = {'overwriteDax', 'boolean', false}; % 

% DecodePixels pars
defaults(end+1,:) = {'minDistFromZeros', 'nonnegative', .001}; % .004 Distance from all zero bit (speeds things up to make this cut before brightness cut0.
defaults(end+1,:) = {'minSeparationFromNext', 'nonnegative', .0001}; % Distance from next nearest bit
defaults(end+1,:) = {'quantileBlankBitRatioCut', 'fraction', .9}; % 
defaults(end+1,:) = {'quantileGeneBrightnessCut', 'fractions', .3}; % 
defaults(end+1,:) = {'borderSize', 'nonnegative', 15*upsample}; %  % down scale this if we drop upsampling.
defaults(end+1,:) = {'minPixels', 'nonnegative', ceil(5*(upsample/2))}; %  % down scale this if we drop upsampling.
defaults(end+1,:) = {'calcSpotMat', 'boolean', false}; %  
defaults(end+1,:) = {'rerun', 'boolean', false}; %  

% defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
% defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
% defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
% defaults(end+1,:) = {'verbose', 'boolean', true}; % 
defaults(end+1,:) = {'moreVerbose', 'boolean', true}; % 
defaults(end+1,:) = {'overwrite', 'boolean', true}; % 
pars = ParseVariableArguments([], defaults);

%% Loop over all 
tic
for d = 1:length(dataToAnalyze); % d = 3
    [basePath,folderName] = fileparts( dataToAnalyze(d).dataPath);
    dataPath = [dataToAnalyze(d).dataPath,filesep];
    
    newPath = [MERFISHdrive,'AdditionalAnalysis\',folderName,'\'];
    dataToAnalyze(d).processedDataFolder = SetFigureSavePath( [newPath,'ProcessedData\'],'makeDir',true); 
    dataToAnalyze(d).analyzedDataFolder = SetFigureSavePath([newPath,'AnalyzedData\'],'makeDir',true); 

    disp(['Analyzing data: ',folderName,'. Dataset ',num2str(d),' of ',num2str(length(dataToAnalyze))]);
    disp(['total runtime = ',num2str(toc/60/60,2),' hours']);
    disp(['saving processed data to: ',dataToAnalyze(d).processedDataFolder]);
    disp(['saving analysis to: ',dataToAnalyze(d).analyzedDataFolder]);
    
    
    %% Find all raw daxfiles
    % Parameters for parsing file names
    dataTag = dataToAnalyze(d).imageTag;
    defaults1 = cell(0,3);
    defaults1(end+1,:) = {'imageTag', 'string', dataTag}; %#ok<*SAGROW> % Base tag for all images
    defaults1(end+1,:) = {'fileExt', 'string', 'dax'}; % 
    defaults1(end+1,:) = {'fieldNames', 'cell', {'movieType', 'hybNum', 'cellNum', 'isFiducial', 'binType'}}; 
    defaults1(end+1,:) = {'fieldConv', 'cell', {@char, @str2num, @str2num, @(x)strcmp(x, 'c2'), @char}};
    parsForFindFiles = ParseVariableArguments([], defaults1);
    foundFiles = BuildFileStructure(dataPath,'parameters', parsForFindFiles,'requireFlag',dataTag); % 
    if isempty(foundFiles);
        error('No files found! Check to see data paths are entered correctly.');
    end

    %% Load dax files
    [movieN,infoPerCell,movieFid,pars] = LoadMERFISHdax(foundFiles,'saveDataPath',dataToAnalyze(d).processedDataFolder,'imageTag',dataTag,'parameters',pars);

    %% Compute fiducial warp and background leveling, save new dax file
    [daxNormName,fidName,daxNoNormName,pars] = SaveCellImagesByBit(movieN,infoPerCell,movieFid,dataToAnalyze(d).processedDataFolder,'parameters',pars);
    
    %% Identify RNA and save cell data and analysis images in current figure 'saveFolder'
    DecodePixels(dataToAnalyze(d),'parameters',pars);    
    
    %% Save a copy of this Script with data
    saveMfileName = [dataToAnalyze(d).analyzedDataFolder,mfilename,'.m'];
    saveMfileName = IncrementSaveName(saveMfileName);
    try
        copyfile( [mfilename('fullpath'),'.m'],saveMfileName);
        display('------------------------------------------------------------------');
        display(['Copied analysis script to ' saveMfileName]);
    catch
        warning('Could not find mfile path');
    end
end