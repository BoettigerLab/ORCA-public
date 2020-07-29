function response = Git(gitString, varargin)
%--------------------------------------------------------------------------
% response = Git(gitString, varargin)
% This function is an interface between git and matlab. All standard git
% commands can be passed to this function
%--------------------------------------------------------------------------
% Outputs:
%
% response/string: The output returned by git
%
%--------------------------------------------------------------------------
% Inputs:
% 
% gitString/string: The command to send to git
% 
%--------------------------------------------------------------------------
% Variable Inputs:
% 
% 'gitDirectory': Path to the git directory. The default is the first 
%    registered directory
% 'repoTag'/string([]): A string tag used by regexp to identify a registered
%    git repo.
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% January 14, 2014
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global gitPath;
global gitRepositoryPaths

%--------------------------------------------------------------------------
% Define default git directory: (The first registered directory)
%--------------------------------------------------------------------------
if isempty(gitRepositoryPaths)
    error('No registered git repositories');
end
gitDirectory = gitRepositoryPaths{1};

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
repoTag = [];

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
            case 'gitDirectory'
                gitDirectory = CheckParameter(parameterValue, 'string', parameterName);
            case 'repoTag'
                repoTag = CheckParameter(parameterValue, 'string', parameterName);
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end
%--------------------------------------------------------------------------
% Switch directory if needed
%--------------------------------------------------------------------------
if ~isempty(repoTag)
    repoInd = find(cellfun(@(x) ~isempty(regexp(lower(x), lower(repoTag))), ...
        gitRepositoryPaths));
    if isempty(repoInd)
        error(['Could not find repository containing tag: ' repoTag]);
    end
    gitDirectory = gitRepositoryPaths{repoInd(1)}; % Use first match
end

%--------------------------------------------------------------------------
% Compose nonlocal git call
%--------------------------------------------------------------------------
gitCommand = ['"' gitPath 'git"' ' --git-dir=' gitDirectory '.git' ' --work-tree=' gitDirectory];

%--------------------------------------------------------------------------
% Add git command
%--------------------------------------------------------------------------
switch gitString
    case 'current commit'
        gitString = 'rev-parse HEAD';
    case 'current branch'
        gitString = 'rev-parse --abbrev-ref HEAD';
    case 'commit message'
        gitString = 'log -1 --pretty=%B';
end
%--------------------------------------------------------------------------
% Determine if dos prompt is required
%--------------------------------------------------------------------------
ampersand = '';
if ~isempty(regexp(gitString, 'push'))
    ampersand = ' &';
end

%--------------------------------------------------------------------------
% Send git command
%--------------------------------------------------------------------------
[~, response] = dos([gitCommand ' ' gitString ampersand]);
