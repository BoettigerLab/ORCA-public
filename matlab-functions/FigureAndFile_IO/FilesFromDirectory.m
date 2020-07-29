function fileName = FilesFromDirectory(directories, varargin)
% ------------------------------------------------------------------------
% fileName = FilesFromDirectory(directories, varargin)
% This function returns a list of file names of the file type provided 
% contained within the provided directory or directories 
%--------------------------------------------------------------------------
% Necessary Inputs
% directories/string or cell array: Path to a valid directory or
%    directories
%
%--------------------------------------------------------------------------
% Outputs
% fileNames/cell array: Cell array of found dax file names
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 'longName'/boolean(true): Determines if the returned dax file names
%    contain the full path (long) or not (short)
% 'exclude'/cell array: A set of regexp strings to exclude
% 'required'/cell array: A set of regexp strings that must be present in
%   all files that are returned. 
% 'required'/cell array: A set of regexp strings for which at least one 
%   entry must be present in all files that are returned. 
% 'fileExt'/string('dax'): The default extension for the files to find
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% March 25, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Define default values
% ------------------------------------------------------------------------
longName = true;
exclude = {};
required = {};
fileExt = 'dax';
requiredOr = {};

% ------------------------------------------------------------------------
% Parse variable arguments
% ------------------------------------------------------------------------
if nargin < 1
    error([mfilename ' requires at least one directory']);
end
if nargin > 2
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'longName'
                longName = CheckParameter(parameterValue, 'boolean', parameterName);
            case 'exclude'
                exclude = CheckParameter(parameterValue, 'cell', parameterName);
            case 'required'
                required = CheckParameter(parameterValue, 'cell', parameterName);
            case 'requiredOr'
                requiredOr = CheckParameter(parameterValue, 'cell', parameterName);
            case 'fileExt'
                fileExt = CheckParameter(parameterValue, 'string', parameterName);
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename])
        end
    end
end

% ------------------------------------------------------------------------
% Converting provided directory to cell
% ------------------------------------------------------------------------
if ~iscell(directories)
    temp = directories;
    directories = {};
    directories{1} = temp;
end

% ------------------------------------------------------------------------
% Check directories
% ------------------------------------------------------------------------
for i=1:length(directories)
    if ~(exist(directories{i}) == 7)
        warning([directories{i} ' is not a valid directory']);
    end
    temp = directories{i};
    if temp(end) ~= filesep;
        temp(end+1) = filesep;
        directories{i} = temp;
    end
end

% ------------------------------------------------------------------------
% Find all dax (or other) files
% ------------------------------------------------------------------------
fileName = {};
for i=1:length(directories)
    dirData = dir(directories{i});
    extInds = cellfun(@(x) ~isempty(regexp(x, ['.' fileExt])), {dirData.name});
    
    % Excluded features
    excludeInds = false(length(exclude), length(extInds));
    for j=1:length(exclude)
        excludeInds(j, :) = cellfun(@(x) ~isempty(regexp(x, exclude{j})), {dirData.name});
    end
    excludeInds = any(excludeInds, 1);
    
    % Required features
    requiredInds = false(length(required), length(extInds));
    for j=1:length(required)
        requiredInds(j,:) = cellfun(@(x) ~isempty(regexp(x, required{j})), {dirData.name});
    end
    requiredInds = all(requiredInds, 1);
    
    % Required features but with or logic
    requiredOrInds = false(length(requiredOr), length(extInds));
    for j=1:length(requiredOr)
        requiredOrInds(j,:) = cellfun(@(x) ~isempty(regexp(x, requiredOr{j})), {dirData.name});
    end
    if ~isempty(requiredOrInds)
        requiredOrInds = any(requiredOrInds, 1);
    else
        requiredOrInds = true(1, length(extInds));
    end
    
    indsToKeep = (extInds & requiredInds & ~excludeInds & requiredOrInds);
    
    foundFileNames = {dirData(indsToKeep).name};
    for j=1:length(foundFileNames)
        if longName
            fileName{end+1} = [directories{i} foundFileNames{j}];
        else
            fileName{end+1} = foundFileNames{j};
        end
    end
end

    