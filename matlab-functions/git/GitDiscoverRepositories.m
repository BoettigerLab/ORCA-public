function gitRepoPaths = GitDiscoverRepositories(varargin)
%--------------------------------------------------------------------------
% gitRepoPaths = GitDiscoverRepositories(varargin)
% This function searches the current matlab path for git repositories and
% returns a cell array of paths for the found repositories
%--------------------------------------------------------------------------
% Outputs:
% gitRepoPaths/cell array: The found paths to git repositories
%
%--------------------------------------------------------------------------
% Necessary Inputs:
% 
%--------------------------------------------------------------------------
% Variable Inputs:
% 
% 'verbose'/boolean(False): Display search information
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% March 28, 2014
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = false;
%--------------------------------------------------------------------------
% Define default parameters
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
                verbose = CheckParameter(parameterValue, 'boolean', parameterName);
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Compile current matlab search paths
%--------------------------------------------------------------------------
pathString = path; 
currentPaths = strread(pathString,'%s','delimiter',';'); 

%--------------------------------------------------------------------------
% Search for .git repositories
%--------------------------------------------------------------------------
folderID = 7;
containsGit = find(cellfun(@(x) exist([x '\.git']) == folderID, currentPaths));

gitRepoPaths = {};
for i=1:length(containsGit)
    gitRepoPaths{end+1} = [currentPaths{containsGit(i)} '\'];
end

%--------------------------------------------------------------------------
% Report results
%--------------------------------------------------------------------------
if verbose
    display('--------------------------------------------------------------');
    display(['Found ' num2str(length(containsGit)) ' repositories in ' num2str(length(currentPaths)) ' paths']);
    for j=1:length(gitRepoPaths)
        display(['    ' gitRepoPaths{j}]);
    end
end
