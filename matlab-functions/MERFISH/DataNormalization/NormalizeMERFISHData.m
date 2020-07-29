function [foundFiles, parameters] = NormalizeMERFISHData(dataPath, outputPath, varargin)
% ------------------------------------------------------------------------
% [foundFiles, parameters] = NormalizeMERFISHData(dataPath, varargin)
% This function handles the normalization of MERFISH data contained within
% the specified dataPath. 
%--------------------------------------------------------------------------
% Necessary Inputs: dataPath
%--------------------------------------------------------------------------
% Outputs: 
%   foundFiles -- an array of file structures as constructed by
%       BuildFileStructure
%--------------------------------------------------------------------------
% External files written:
% This function will write files to the specified directory. These will be
% files of the form: bitName_data_allFOV.dax and the corresponding .inf
%
% In addition, for each combined dax, a .pos file will be written. This is
% a csv file in which the first two columns correspond to the stage x and y
% positions for each frame (FOV) in the corresponding dax file. The third
% and fourth columns correspond to the actual z position in volts and the
% autofocus target, respectively. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%   'schemaPath' -- a path to a valid data scheme xls
%   'dataLayouts' -- a structure array. The acceptable fields are described
%       in LoadDataSchema.m
%--------------------------------------------------------------------------
% Guiping Wang and Jeffrey Moffitt
% guiping.w.d@gmail.com
% jeffmoffitt@gmail.com
% July 21, 2015
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for setting up data schema
defaults(end+1,:) = {'schemaPath', 'filePath', []};

% Parameters for the default parsing of data names: Parameters for
% BuildFileStructure
defaults(end+1,:) = {'fileExt', 'string', 'dax'}; % File extension of identified files
defaults(end+1,:) = {'delimiters', 'cell', {'_'}}; % Delimiters to use to split the dax names
defaults(end+1,:) = {'fieldNames', 'cell', {'movieType', 'fovID', 'hybID'}}; % FieldNames
defaults(end+1,:) = {'fieldConv', 'cell', {@char, @(x)int32(str2num(x)), @(x)int32(str2num(x))}}; % Conversion functions
defaults(end+1,:) = {'appendExtraFields', 'boolean', true}; 
defaults(end+1,:) = {'containsDelimiters','positive', 1}; % An integer specifying a field that might have internal delimiters, such as the movieType

% Parameters for the dataLayout
defaults(end+1,:) = {'dataLayouts', 'array', []}; % Input for predefined data layout structures. See LoadDataSchema.

% Parameters to handle existing data
defaults(end+1,:) = {'overwrite', 'boolean', false}; % Do not overwrite by default

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2 || ~(exist(dataPath) == 7) % 7=folder
    error('matlabFunctions:invalidArguments', 'A valid data path is required.');
end

if ~exist(outputPath) == 7
    warning('matlabFunctions:path', 'Specified output path does not exist. It will be created.');
    if ~mkdir(outputPath)
        error('matlabFunctions:invalidArguments', 'The specified output path cannot be created');
    end
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Define default data schema path if not provided and if no data layouts
% are provided
% -------------------------------------------------------------------------
if isempty(parameters.schemaPath) && isempty(parameters.dataLayouts)
    parameters.schemaPath = [dataPath 'dataSchema.xls'];
    % Check for the presence of the default schema
    if exist(parameters.schemaPath) ~= 2 % 2=file
        parameters.schemaPath = [parameters.schemaPath 'x']; % Test for xlsx version
        if exist(parameters.schemaPath) ~= 2
            error('matlabFunctions:dataNormalization', 'No data schema was provided or found');
        end
    end    
end

% -------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
if parameters.verbose
    PageBreak();
    display('Normalizing MERFISH Data');
    PageBreak();
    mainTimer = tic; % Start main progress timer
end

% -------------------------------------------------------------------------
% Load data schema if no data layouts are provided
% -------------------------------------------------------------------------
if isempty(parameters.dataLayouts)
    % Display progress
    PageBreak();
    display(['Loading data schema from ' parameters.schemaPath]);
    
    % Load data layouts and schema parameters
    [parameters.dataLayouts, schemaParameters] = LoadDataSchema(parameters.schemaPath);
    
    % Transfer NAME_SCHEMA parameters
    parameters.fieldNames = schemaParameters.fieldNames;
    parameters.fieldConv = schemaParameters.fieldConv;
end
% Copy parameters.dataLayouts to local variable for convenience
dataLayouts = parameters.dataLayouts;

% -------------------------------------------------------------------------
% Build file structure for files in data path
% -------------------------------------------------------------------------
foundFiles = BuildFileStructure(dataPath, 'parameters', parameters);

if parameters.verbose
    display(['Found ' num2str(length(foundFiles)) ' files in ' ...
        dataPath]);
end

% -------------------------------------------------------------------------
% Loop over data structures
% -------------------------------------------------------------------------
for i=1:length(dataLayouts)
    % Display progress
    timer = tic;
    display(['Writing ' dataLayouts(i).dataType ' bit ' dataLayouts(i).bit]);
    
    % Identify files that fit this format    
    dataToLoad = foundFiles(strcmp({foundFiles.movieType}, dataLayouts(i).movieType) & ... % Find identical movieType
        [foundFiles.hybID] == dataLayouts(i).hybID); % Find identical hybID
    
    % Handle case if no data are found
    if isempty(dataToLoad)
        warning('matlabFunctions:fileNotFound', 'No data files were found for this bit and data type');
        continue; % Skip the remaining part of this execution of this loop
    end
    
    % Sort based on fovID
    [~, sind] = sort([dataToLoad.fovID], 'ascend');
    dataToLoad = dataToLoad(sind);
        
    % Construct new inf file
    infFile = ReadInfoFile(dataToLoad(1).filePath);
    newInfFile = infFile;
    newInfFile.localPath = outputPath;
    newInfFile.localName = [dataLayouts(i).bit '_' dataLayouts(i).dataType '_allFOV.inf'];
    
    % Check to see if file already exists
    if ~parameters.overwrite
        if exist([outputPath newInfFile.localName]) == 2
            warning('matlabFunctions:existingFile', 'Found existing file. Skipping.');
            continue;
        end
    end
    
    % Determine subregion parameters
    subregion = zeros(4,1);
    if isfield(dataLayouts(i), 'ROI')
        % RESERVED FOR FUTURE USE
    end
    
    % Load first frame to determine memory requirements
    dax = ReadDax(dataToLoad(1).filePath, ...
        'startFrame', dataLayouts(i).frame, ...
        'endFrame', dataLayouts(i).frame, ...
        'subregion', subregion, ...
        'verbose', false);
    
    % Allocate memory
    newDax = repmat(dax, [1 1 length(dataToLoad)]);
    
    % Preallocate position array
    stagePos = zeros(length(dataToLoad),4);
    
    % Loop over all dax
    for j=1:length(dataToLoad) % Reload the first frame to simplify code
        % Load dax
        [newDax(:,:,j), localInf] = ReadDax(dataToLoad(j).filePath, ...
            'startFrame', dataLayouts(i).frame, ...
            'endFrame', dataLayouts(i).frame, ...
            'subregion', subregion, ...
            'verbose', false);
        
        % Load offset file
        offData = tdfread([dataToLoad(j).filePath(1:(end-3)) 'off'], ' ');

        % Save stage position
        stagePos(j,:) = [...
            localInf.Stage_X, ... 
            localInf.Stage_Y, ...
            offData.stage0x2Dz(dataLayouts(i).frame), ...
            localInf.Lock_Target];
    end
    
    % Update dimensions of the new dax
    newInfFile.number_of_frames = size(newDax,3);
    
    % UNDER CONSTRUCTION %
    % Need to handle potention cropping of dax %
    
    % Write new dax
    WriteDAXFiles(newDax, newInfFile, 'verbose', true);
    
    % Write stage pos
    positionFile = [dataLayouts(i).bit '_' dataLayouts(i).dataType '_allFOV.pos'];
    csvwrite([outputPath positionFile], stagePos);
    
    % Display progress
    if parameters.verbose
        display(['... completed in ' num2str(toc(timer)) ' s']);
    end
    
end

% -------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
if parameters.verbose
    PageBreak();
    display(['Completed data normalization in ' num2str(toc(mainTimer)) ' s']);
end


