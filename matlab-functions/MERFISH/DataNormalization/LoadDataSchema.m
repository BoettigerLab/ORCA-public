function [dataLayouts, parameters] = LoadDataSchema(schemaPath, varargin)
% ------------------------------------------------------------------------
% [dataLayouts, parameters] = LoadDataSchema(schemaPath, varargin)
% This function loads a data schema file, extracting parameters associated
% with parsing file names as well as parsing various types of data. 
%--------------------------------------------------------------------------
% Necessary Inputs: 
%   schemaPath - a path to a valid data scheme xls(x) file.
%--------------------------------------------------------------------------
% Outputs: 
%   dataLayouts - an array of dataLayout structures. Each structure
%   contains the following fields
%       bit: The name of the bit
%       movieType: The name portion of the movie that contains this bit
%       hybID: The ID number of the hybridization round in which the bit
%           was collected
%       frame: The frame(s) of that movie in which this bit was imaged
%       color: An identifier of the color channel
%       dataType: An identifier specifying the type of data, i.e.
%           fiducials, cell contrast marker such as GFP, or the MERFISH data
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 
%--------------------------------------------------------------------------
% Example Calls
%
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

% Parameters for the default parsing of data names: Parameters for
% BuildFileStructure
defaults(end+1,:) = {'fileExt', 'string', 'dax'}; % File extension of identified files
defaults(end+1,:) = {'delimiters', 'cell', {'_'}}; % Delimiters to use to split the dax names
defaults(end+1,:) = {'fieldNames', 'cell', {'movieType', 'fovID', 'hybID'}}; % FieldNames
defaults(end+1,:) = {'fieldConv', 'cell', {@char, @(x)int32(str2num(x)), @(x)int32(str2num(x))}}; % Conversion functions
defaults(end+1,:) = {'appendExtraFields', 'boolean', true}; 
defaults(end+1,:) = {'containsDelimiters','positive', 1}; % An integer specifying a field that might have internal delimiters, such as the movieType

% Parameters for user updates
defaults(end+1,:) = {'verbose', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~exist(schemaPath) == 2 % 2==file
    error('matlabFunctions:invalidArguments', 'A valid schema path is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Load data schema file
% -------------------------------------------------------------------------
% Display progress
if parameters.verbose
    PageBreak();
    display(['Loading data schema from ' schemaPath]);
end

% Load xls
[~, ~, raw] = xlsread(schemaPath);

% Pad with nan if the final entry is not a nan: needed to parse final
% DATA_LAYOUT
if ~isnan(raw{end,1})
    raw(end+1,:) = num2cell(nan(1, size(raw,2)));
end

% -------------------------------------------------------------------------
% Parse movie name schema if requested
% -------------------------------------------------------------------------
% Find the NAME_SCHEMA entry
nameSchemaInd = find(strcmp(raw(:,1), 'NAME_SCHEMA'));
if length(nameSchemaInd)>1
    nameSchemaInd = nameSchemaInd(1);
    warning('matlabFunctions:dataSchema', 'Found more than one NAME_SCHEMA flag. Using the first');
end

% Load NAME_SCHEMA if found
if ~isempty(nameSchemaInd)

    % Check for all needed fields: fieldNames and fieldTypes
    if ~all([strcmp(raw{nameSchemaInd+1,1}, 'fieldNames') strcmp(raw{nameSchemaInd+2,1}, 'fieldTypes')])
        error('matlabFunctions:dataSchema', 'Corrupted NAME_SCHEMA entry. Could not find either fieldNames or fieldTypes');
    end

    % Transfer field names
    potentialFieldNames = raw(nameSchemaInd+1,2:end);  
    parameters.fieldNames = potentialFieldNames(cellfun(@(x)~any(isnan(x)), potentialFieldNames)); % Remove NaN

    % Parse field types
    fieldTypes = raw(nameSchemaInd+2,2:end);
    convFuncs = {};
    for i=1:length(parameters.fieldNames)
        switch lower(fieldTypes{i})
            case {'char', 'string'}
                convFuncs{i} = @char;
            case {'int'}
                convFuncs{i} = @(x)int32(str2num(x));
            otherwise
                convFuncs{i} = @char;
        end
    end
    parameters.fieldConv = convFuncs;

    % Display progress
    if parameters.verbose
        display(['Loaded name schema']);
    end
else
    display(['Using default name schema']);
end

% -------------------------------------------------------------------------
% Parse data layouts
% -------------------------------------------------------------------------
% Find DATA_LAYOUT flags
dataLayoutInds = find(strcmp(raw(:,1), 'DATA_LAYOUT'));

% Raise error if none are found
if isempty(dataLayoutInds)
    error('matlabFunctions:dataSchema', 'No DATA_LAYOUT flags found in provided data schema');
end

% Display progress
if parameters.verbose
    display(['Found ' num2str(length(dataLayoutInds)) ' data layouts']);
end

% Add Inf to end of the inds--useful for finding the finish of each block
dataLayoutInds(end+1) = size(raw,1);

% Find empty entries
emptyInds = find(cellfun(@(x)any(isnan(x)), raw(:,1)))';

% Prepare empty dataLayouts
dataLayouts = [];

% Loop over found data layouts
for i=1:(length(dataLayoutInds)-1) % -1 because of the added final entry
    % Find the indices corresponding to this DATA_LAYOUT
    startInd = dataLayoutInds(i)+1;
    finishInd = min([dataLayoutInds(i+1) emptyInds(emptyInds>dataLayoutInds(i))])-1;
    
    % Coerce this DATA_LAYOUT into a structure array
    T = cell2table(raw((startInd+1):finishInd,:), 'VariableNames', raw(startInd,:));
    localLayouts = table2struct(T);
    [localLayouts.dataType] = deal(raw{dataLayoutInds(i),2});
    dataLayouts = cat(1,dataLayouts, localLayouts);
end

