function [moleculeLists, binFileNames, parameters] = NormalizeDataForThresholdScan(dataPath, varargin)
% ------------------------------------------------------------------------
% moleculeLists, parameters = NormalizeDataForThresholdScan(dataPath, varargin)
% This function prepares the data in dataPath for threshold scanning. 
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Alistair Boettiger, Jeffrey Moffitt
% November 26, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', false}; 
defaults(end+1,:) = {'imageTag', 'string', []}; 
defaults(end+1,:) = {'printedUpdates', 'boolean', true};

% Parameters for DaoSTORM
defaults(end+1,:) = {'maxCPU', 'freeType', 95}; 
defaults(end+1,:) = {'batchsize', 'integer', 12}; 
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
defaults(end+1,:) = {'parametersDefaultFile', 'string', ''};
% Thresholds to scan with DaoSTORM
defaults(end+1,:) = {'thresholds','array',round(logspace(log10(500),log10(3000),20))};

% Parameters for parsing file names
defaults(end+1,:) = {'fileExt', 'string', 'dax'}; % Delimiters for bin files
defaults(end+1,:) = {'fieldNames', 'cell', {'movieType', 'hybNum', 'cellNum', 'isFiducial'}}; 
defaults(end+1,:) = {'fieldConv', 'cell', {@char, @str2num, @str2num, @(x)strcmp(x, 'c2')}};
defaults(end+1,:) = {'appendExtraFields', 'bool', true};
defaults(end+1,:) = {'fiducialMListType', 'string', 'mlist'};
defaults(end+1,:) = {'justCombineMoleculeLists', 'boolean', false};

% Parameters for saving
defaults(end+1,:) = {'savePath', 'string', []};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~(exist(dataPath) == 7) % 7=folder
    error('matlabFunctions:invalidArguments', 'A valid data path is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Rename paths
% -------------------------------------------------------------------------
thresholds = parameters.thresholds; % Threshold values to try
if isempty(parameters.savePath)
    parameters.savePath = dataPath;
end
newPath = parameters.savePath; 

% -------------------------------------------------------------------------
% Find files in data path
% -------------------------------------------------------------------------
foundFiles = BuildFileStructure(dataPath, 'parameters', parameters,...
    'requireFlag',parameters.imageTag);

% -------------------------------------------------------------------------
% Determime properites
% -------------------------------------------------------------------------
isImageData = strcmp({foundFiles.movieType},parameters.imageTag);
numHybs = length(unique([foundFiles.hybNum]));
numCells = length(unique([foundFiles.cellNum]));

% Catch problems with fiducials 
for i=1:length(foundFiles)
    if isempty(foundFiles(i).isFiducial)
        foundFiles(i).isFiducial = 0;
    end
end

% -------------------------------------------------------------------------
% Create composite dax files
% -------------------------------------------------------------------------
if parameters.printedUpdates
    display('-------------------------------------------------------------');
    display('Creating composite dax');
end
compositeDaxTimer = tic;

% Prepare cell arrays for the names of these files 
daxNames = cell(numHybs,1); 
hybNames = cell(numHybs,1); 

% Loop over the number of hybridization/imaging rounds
for h=1:numHybs  
    % Create composite dax name
    daxNames{h} = ['AllCells_Hybe_',num2str(h),'.dax'];
    
    % Check to see if file exists and skip if overwrite is false
    if exist([newPath daxNames{h}])
        if parameters.verbose
            display(['Found existing file: ' daxNames{h}]);
            if parameters.overwrite
                display(['... overwriting']);
            else
                display(['... skipping file']);
                continue; % Break the execution of this iteration of the loop
            end
        end
    end
    
    % Index files that will be part of this image stack
    hybeN = [foundFiles.hybNum]==h-1;
    notFid = [foundFiles.isFiducial] ~= true;
    hybNames{h} = {foundFiles.name}';
    hybNames{h} = hybNames{h}(hybeN & isImageData & notFid );

    numCells = length(hybNames{h});

    % Load the first movie to determine the needed memory
    if h==1 
        tempImage = ReadDax([dataPath, hybNames{h}{1}], 'verbose', false);
        movieN = zeros(size(tempImage,1), size(tempImage,2), numCells);
    end
    
    % Loop over all FOV/cells, load, and insert into stack
    for c=1:numCells
        [movieN(:,:,c),infoN] = ReadDax([dataPath,hybNames{h}{c}],'verbose',false); 
    end
    
    % Define properties of this stack for the new information file
    infoN.number_of_frames = numCells;
    infoN.localPath = newPath; 
    infoN.localName = regexprep(daxNames{h},'.dax','.inf'); 
    
    % Write the stack in dax format
    WriteDAXFiles(movieN,infoN, ...
        'verbose', parameters.verbose); 
end

% Display progress
elapsedTime = toc(compositeDaxTimer);
if parameters.printedUpdates
    display(['... finished in ' num2str(elapsedTime) ' s']);
end

% -------------------------------------------------------------------------
% Create Parameter Files for Dax
% -------------------------------------------------------------------------
numThresholds = length(thresholds);
paramsNames = cell(numThresholds,1);
for t=1:numThresholds
    fitParams = ReadParsFile(parameters.parametersDefaultFile);
    fitParams.threshold = num2str(thresholds(t));
    paramsNames{t} = [newPath,'convpars_t',num2str(thresholds(t)),'.xml'];
    if ~exist(paramsNames{t})
        WriteParsFile(paramsNames{t} ,fitParams);
    end
end

% -------------------------------------------------------------------------
% Run Spot Finder
% -------------------------------------------------------------------------
% Display progress
if parameters.printedUpdates
    display('-------------------------------------------------------------');
    display('Fitting molecules');
end
dotFinderTimer = tic;
for t = 1:numThresholds
%     RunDotFinder('daxnames', strcat(newPath, daxNames),...
%                  'parsfile', paramsNames{t},...
%                  'binname',['DAX_',num2str(thresholds(t))],...
%                  'maxCPU',parameters.maxCPU,...
%                  'batchsize',parameters.batchsize,...
%                  'hideterminal',false,...
%                  'overwrite',parameters.overwrite);
             
   % AnalyzeSTORM
   AnalyzeSTORM('path', newPath, ...
       'config', paramsNames{t}, ...
       'method', 'daoSTORM', ...
       'verbose', parameters.verbose, ...
       'overwrite', parameters.overwrite, ...
       'numParallel', parameters.batchsize, ...
       'outputInMatlab', false, ...
       'prefix', num2str(thresholds(t))); 
end
% Display progress
elapsedTime = toc(dotFinderTimer);
if parameters.printedUpdates
    display(['... finished in ' num2str(elapsedTime) ' s']);
end

% -------------------------------------------------------------------------
% Load molecule lists
% -------------------------------------------------------------------------
if parameters.printedUpdates
    display('-------------------------------------------------------------');
    display('Loading and combining molecule lists');
end
moleculeListTimer = tic;

display('loading molecule lists...'); 
binFileNames = cell(numHybs, numThresholds);
moleculeLists = cell(numCells,numHybs,numThresholds); 
for t =1:numThresholds;
    if parameters.verbose
        display(['Loading threshold: ' num2str(thresholds(t)) ]);
    end
    for h = 1:numHybs
        binName = regexprep(  daxNames{h}, '\.dax',['_',num2str(thresholds(t)),'_mlist\.bin']);
        try
            imLists = ReadMasterMoleculeList([newPath,binName],'verbose',false, ...
                'fieldsToLoad', {'x', 'y', 'a', 'xc', 'yc', 'frame'}, ...
                'loadAsStructArray', true);
            binFileNames{h,t} = binName;
        catch er
            warning(er.getReport);
        end
        for c=1:numCells
            moleculeLists{c,h,t} = imLists(imLists.frame == c);
        end
    end
end

moleculeTime = toc(moleculeListTimer);
if parameters.printedUpdates
    display(['... finished in ' num2str(moleculeTime) ' s']);
end

% -------------------------------------------------------------------------
% Compute and Apply Warps
% -------------------------------------------------------------------------
if parameters.printedUpdates
    display('-------------------------------------------------------------');
    display('Calculating and applying warps');
end
warpTimer = tic;
for c=1:numCells
    % Aligning fiducial data
    fiducialData = struct;
    fiducialData(numHybs).mList = []; 
    for h=1:numHybs
        fedNames = {foundFiles.name}'; 
        hybeN = [foundFiles.hybNum]==h-1;  % belongs to this hybe h
        isFid = [foundFiles.isFiducial] == true; % is a fiducial
        currentCell = [foundFiles.cellNum] == c-1; % belongs to current cell
        hybeName = fedNames{hybeN & isImageData & isFid & currentCell};
        beadBinName = regexprep(hybeName,'\.dax',...
            ['_' parameters.fiducialMListType '\.bin']);
        try
            fiducialData(h).mList = ReadMasterMoleculeList([dataPath,beadBinName],'verbose',false, ...
                'fieldsToLoad', {'xc', 'yc', 'frame', 'w'});
            fiducialData(h).cellNum = c-1; % Mark cell number for report generation in AlignFiducials
        catch er
            warning(er.getReport);
        end
    end
    
    [fiducialData, parameters] = AlignFiducials(fiducialData, ...
        'parameters', parameters);
    
    % Warp mLists
    for t=1:numThresholds
        for h=1:numHybs
            localList = moleculeLists{c,h,t};
            if ~isempty(localList)
                [localList.xc, localList.yc] = tforminv(fiducialData(h).tform, double(localList.x), double(localList.y));
                moleculeLists{c,h,t} = localList;
            end
        end
    end
end

warpTime = toc(warpTimer);
if parameters.printedUpdates
    display(['... finished in ' num2str(warpTime) ' s']);
end

