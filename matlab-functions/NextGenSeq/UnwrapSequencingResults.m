function UnwrapSequencingResults(varargin)
%--------------------------------------------------------------------------
% UnwrapSequencingResults(varargin)
% Unwrap the folder structure from Beckman Coulter Genomics
%--------------------------------------------------------------------------
% Necessary Inputs
%
%--------------------------------------------------------------------------
% Outputs
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% May 24, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = false;

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'verbose'
                verbose = CheckParameter(parameterValue,'boolean','verbose');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Find folder
%--------------------------------------------------------------------------
unwrapPath = uigetdir();

dataPaths = {};
if includeSubdir
    dirContents = dir(analysisPath);
    dataPaths = {dirContents([dirContents.isdir]).name};
    dataPaths = dataPaths(~ismember(dataPaths, {'.', '..'}));
    % Build full path
    for i=1:length(dataPaths)
        dataPaths{i} = [analysisPath dataPaths{i} '\'];
    end
end
dataPaths{end+1} = analysisPath; % Append analysisPath
