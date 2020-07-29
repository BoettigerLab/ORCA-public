function gitStates = GitState(varargin)
%--------------------------------------------------------------------------
% function gitStates = GitState(varargin)
% This function returns an array of status objects for each of the
% registered git repositories. These objects contain the repository name,
% its local location, the current branch name, the current commit hash, and
% the status of the directory, i.e. clean or not. 
%--------------------------------------------------------------------------
% Outputs:
% gitStates/structure array: A structure array of git state structures. Each
%   structure can contain the following fields:
%   --repoName: The folder name containing the repository
%   --repoPath: The local path of the repository
%   --branch: The name of the current branch
%   --commit: The hash of the current commit
%   --commitMessage: The current commit message
%   --clean: A boolean representing whether there were no
%     uncommited changes
%   By default repoName, branch, commit, and clean are included
%--------------------------------------------------------------------------
% Variable Inputs:
% 'fields'/cell array: List of the desired fields to include instead of the
%   default fields
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% March 29, 2014
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global gitRepositoryPaths;
if isempty(gitRepositoryPaths)
    warning('matlabFunctions:noGitRepositories', 'No registered git repositories');
    gitStates = [];
end

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
fields = {'repoName', 'branch', 'commit', 'clean'};
possibleFields = {'repoName', 'repoPath', 'branch', 'commit', 'commitMessage', 'clean'};
verbose = false;
allFields = false;

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
            case 'fields'
                fields = CheckParameter(parameterValue, 'cell', parameterName);
            case 'allFields'
                allFields = CheckParameter(parameterValue, 'boolean', parameterName);
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end
%--------------------------------------------------------------------------
% Handle all fields request
%--------------------------------------------------------------------------
if allFields
    fields = possibleFields;
end

%--------------------------------------------------------------------------
% Loop over git repositories
%--------------------------------------------------------------------------
for i=1:length(gitRepositoryPaths)
    gitState = [];
    for j=1:length(fields)
        switch fields{j}
            case 'repoName'
                parts = regexp(gitRepositoryPaths{i}, filesep, 'split');
                parts = parts(cellfun(@(x) ~isempty(x), parts));
                gitState.repoName = parts{end};
            case 'repoPath'
                gitState.repoPath = gitRepositoryPaths{i};
            case 'branch'
                gitState.branch = strtrim(Git('current branch', 'gitDirectory', gitRepositoryPaths{i}));
            case 'commit'
                gitState.commit = strtrim(Git('current commit', 'gitDirectory', gitRepositoryPaths{i}));
            case 'commitMessage'
                gitState.commitMessage = strtrim(Git('commit message', 'gitDirectory', gitRepositoryPaths{i}));
            case 'clean'
                gitState.clean = isempty(Git('status --porcelain', 'gitDirectory', gitRepositoryPaths{i}));
            otherwise
                error([fields{j} ' is not a valid field for a git status structure']);
        end
    end
    gitStates(i) = gitState;
end

%--------------------------------------------------------------------------
% Display results
%--------------------------------------------------------------------------
if verbose
    PageBreak();
    for i=1:length(gitStates)
        if isfield(gitStates(i), 'repoName')
            display(['Repository Name: ' gitStates(i).repoName]);
        end
        if isfield(gitStates(i), 'repoPath')
            display(['   Path: ' gitStates(i).repoPath]);
        end
        if isfield(gitStates(i), 'branch')
            display(['   Current Branch: ' gitStates(i).branch]);
        end 
        if isfield(gitStates(i), 'commit')
            display(['   Commit: ' gitStates(i).commit]);
        end 
        if isfield(gitStates(i), 'commitMessage')
            display(['   Message: ' gitStates(i).commitMessage]);
        end 
        if isfield(gitStates(i), 'clean')
            if gitStates(i).clean
                display(['   Status: Up to date'])
            else
                display(['   Status: Contains uncommitted changes'])
            end
        end
    end
end

