function [moleculeLists, binFileNames, parameters] = FitMultiThresholdsForScan(dataPath, varargin)
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
defaults(end+1,:) = {'overwrite', 'integer', 0}; 
defaults(end+1,:) = {'parsfile', 'string', '\\MORGAN\TSTORMdata2\AdditionalAnalysis\141021_ParametersFiles\default_convpars.xml'}; 

% Thresholds to scan with DaoSTORM
defaults(end+1,:) = {'thresholds','array',round(logspace(log10(150),log10(2500),30))};

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
STORMmovies = strcmp({foundFiles.movieType},parameters.imageTag);
numHybes = length(unique([foundFiles.hybNum]));
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

% Create and write composite dax 
daxNames = cell(numHybes,1); 
hybeNames = cell(numHybes,1); 
for h=1:numHybes  
    daxNames{h} = ['AllCells_Hybe_',num2str(h),'.dax'];
    fileExists = exist([newPath,daxNames{h}],'file') == 2;  
    hybeN = [foundFiles.hybNum]==h-1;
    notFid = [foundFiles.isFiducial] ~= true;
    hybeNames{h} = {foundFiles.name}';
    hybeNames{h} = hybeNames{h}(hybeN & STORMmovies & notFid );

    info1 = ReadInfoFile([dataPath,hybeNames{h}{1}],'verbose',false);
    frameH = info1.frame_dimensions(1); 
    frameW = info1.frame_dimensions(2); 
    
    if (~fileExists || parameters.overwrite == 1)
        numCells = length(hybeNames{h});
        movieN = zeros(frameH,frameW,numCells); 
        for c=1:numCells
            [movieN(:,:,c),infoN] = ReadDax([dataPath,hybeNames{h}{c}],'verbose',false); 
        end
        infoN.number_of_frames = numCells;
        infoN.localPath = newPath; 
        infoN.localName = regexprep(daxNames{h},'.dax','.inf'); 
        WriteDAXFiles(movieN,infoN, 'verbose', parameters.verbose); 
   end
end

compositeTime = toc(compositeDaxTimer);
if parameters.printedUpdates
    display(['... finished in ' num2str(compositeTime) ' s']);
end

% -------------------------------------------------------------------------
% Create Parameter Files for Dax
% -------------------------------------------------------------------------
numThresholds = length(thresholds);
if ~parameters.justCombineMoleculeLists
    parsNames = cell(numThresholds,1);
    for t=1:numThresholds
        fitPars = ReadParsFile(parameters.parsfile);
        fitPars.threshold = num2str(thresholds(t));
        parsNames{t} = [newPath,'convpars_t',num2str(thresholds(t)),'.xml'];
        if ~exist(parsNames{t})
            WriteParsFile(parsNames{t} ,fitPars);
        end
    end
end

% -------------------------------------------------------------------------
% Run Spot Finder
% -------------------------------------------------------------------------
if ~parameters.justCombineMoleculeLists
    for t = 1:numThresholds
        RunDotFinder('daxnames', strcat(newPath, daxNames),...
                     'parsfile', parsNames{t},...
                     'binname',['DAX_',num2str(thresholds(t))],...
                     'maxCPU',parameters.maxCPU,...
                     'batchsize',parameters.batchsize,...
                     'hideterminal',true,...
                     'overwrite',parameters.overwrite);
    end
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
binFileNames = cell(numHybes, numThresholds);
moleculeLists = cell(numCells,numHybes,numThresholds); 
for t =1:numThresholds;
    if parameters.verbose
        display(['Loading threshold: ' num2str(thresholds(t)) ]);
    end
    for h = 1:numHybes
        binName = regexprep(  daxNames{h}, '\.dax',['_',num2str(thresholds(t)),'_mlist\.bin']);
        try
            imLists = ReadMasterMoleculeList([newPath,binName],'verbose',false, ...
                'fieldsToLoad', {'x', 'y', 'a', 'xc', 'yc', 'frame'});
            binFileNames{h,t} = binName;
        catch er
            warning(er.getReport);
        end
        for c=1:numCells
            moleculeLists{c,h,t} = IndexStructure(imLists,imLists.frame == c);
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
    fiducialData(numHybes).mList = []; 
    for h=1:numHybes
        fedNames = {foundFiles.name}'; 
        hybeN = [foundFiles.hybNum]==h-1;  % belongs to this hybe h
        isFid = [foundFiles.isFiducial] == true; % is a fiducial
        currentCell = [foundFiles.cellNum] == c-1; % belongs to current cell
        hybeName = fedNames{hybeN & STORMmovies & isFid & currentCell};
        beadBinName = regexprep(hybeName,'\.dax','_list\.bin');
        try
            fiducialData(h).mList = ReadMasterMoleculeList([dataPath,beadBinName],'verbose',false, ...
                'fieldsToLoad', {'xc', 'yc', 'frame', 'w'});
        catch er
            warning(er.getReport);
        end
    end
    
    [fiducialData, parameters] = AlignFiducials(fiducialData, ...
        'parameters', parameters);
    
    % Warp mLists
    for t=1:numThresholds
        for h=1:numHybes
            localList = moleculeLists{c,h,t};
            [localList.xc, localList.yc] = tforminv(fiducialData(h).tform, double(localList.x), double(localList.y));
            moleculeLists{c,h,t} = localList;
        end
    end
end

warpTime = toc(warpTimer);
if parameters.printedUpdates
    display(['... finished in ' num2str(warpTime) ' s']);
end

