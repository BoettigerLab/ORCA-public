function foundDirectories = FindAllDirectories(startingDirectory, varargin)
% ------------------------------------------------------------------------
% function directoryNames = FindAllDirectories(startingDirectory)
% This function returns all directoryNames in the startingDirectory
%--------------------------------------------------------------------------
% Necessary Inputs
% startingDirectory/string: The path of the initial directory. 
%
%--------------------------------------------------------------------------
% Outputs
% directoryNames/cell array of paths: All paths to directories within the
%   starting directory
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 8, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Define default values
% ------------------------------------------------------------------------
recursionDepth = Inf;

% ------------------------------------------------------------------------
% Parse variable arguments
% ------------------------------------------------------------------------
if nargin < 1
    error([mfilename ' requires a starting directory']);
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
            case 'recursionDepth'
                recursionDepth = parameterValue;
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename])
        end
    end
end

% ------------------------------------------------------------------------
% Check validity of starting directory
% ------------------------------------------------------------------------
if ~(exist(startingDirectory) == 7)
    error('Provided directory is not a valid directory');
end
if startingDirectory(end) ~= filesep
    startingDirectory(end+1) = filesep;
end

% ------------------------------------------------------------------------
% Check recursion depth
% ------------------------------------------------------------------------
if recursionDepth < 1
    foundDirectories = {};
    return
end

% ------------------------------------------------------------------------
% Find all directories within the starting directory
% ------------------------------------------------------------------------
dirContents = dir(startingDirectory);
directoryNames = {dirContents([dirContents.isdir]).name};
directoryNames = setdiff(directoryNames, {'.', '..'});

% ------------------------------------------------------------------------
% Build full paths
% ------------------------------------------------------------------------
for i=1:length(directoryNames)
    directoryNames{i} = [startingDirectory directoryNames{i} filesep];
end

% ------------------------------------------------------------------------
% Recursively find all directories within these directories
% ------------------------------------------------------------------------
foundDirectories = {};
for i=1:length(directoryNames)
    enclosedDirectories = FindAllDirectories(directoryNames{i}, ...
        'recursionDepth', recursionDepth - 1);
    foundDirectories = [foundDirectories directoryNames{i} enclosedDirectories ];
end

