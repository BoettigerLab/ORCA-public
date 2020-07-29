function GitRegisterRepositories(gitRepoPaths, varargin)
%--------------------------------------------------------------------------
% gitRepoPaths = GitRegisterRepositories(gitRepoPaths, varargin)
% This function registers the provided repository paths with git
%--------------------------------------------------------------------------
% Necessary Inputs:
% gitRepoPaths/cell array: Paths to git repositories
%
%--------------------------------------------------------------------------
% Variable Inputs:
% 
% 'verbose'/boolean(False): Display paths
% 'defaultRepoPath'/string(''): Define the default git repository
%
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
defaultRepoPath = [];
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
            case 'defaultRepoPath'
                defaultRepoPath = CheckParameter(parameterValue, 'string', parameterName);
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Check validity of provided paths
%--------------------------------------------------------------------------
folderID = 7;
for i=1:length(gitRepoPaths)
    if ~(exist([gitRepoPaths{i} '.git']) == folderID)
        error([gitRepoPaths{i} ' does not contain a .git folder']);
    end
end
gitRepoPaths = unique(gitRepoPaths); % Remove duplicates

%--------------------------------------------------------------------------
% Move default path to the front of the list
%--------------------------------------------------------------------------
if ~isempty(defaultRepoPath)
    if ~ismember(defaultRepoPath, gitRepoPaths)
        error('Provided default path is not a member of the provided repository paths');
    end
    gitRepoPaths = setdiff(gitRepoPaths, defaultRepoPath); % Remove path
    gitRepoPaths = fliplr(gitRepoPaths);
    gitRepoPaths{end+1} = defaultRepoPath;
    gitRepoPaths = fliplr(gitRepoPaths);
end

%--------------------------------------------------------------------------
% Register Git Repositories
%--------------------------------------------------------------------------
global gitRepositoryPaths
gitRepositoryPaths = gitRepoPaths;

%--------------------------------------------------------------------------
% Display Progress
%--------------------------------------------------------------------------
if verbose
    display(['Registered ' num2str(length(gitRepoPaths)) ' git repositories']);
    for j=2:length(gitRepoPaths)
        display(['    ' gitRepoPaths{j}]);
    end
    display(['    (default) ' gitRepoPaths{1}]);
end
