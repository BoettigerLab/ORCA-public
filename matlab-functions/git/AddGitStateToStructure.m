function structArray = AddGitStateToStructure(structArray, varargin)
%--------------------------------------------------------------------------
% structArray = AddGitStateToStructure(structArray, varargin)
% This function adds the current state of all registered git repositories 
% to the .git field in each element of the provided structure array.
%
% A warning will be issued if there are uncommitted changes in any of the
% registered repositories. Turn off this warning with 
% warning('off', 'matlabfunctions:GitUncommitedChanges'). 
%--------------------------------------------------------------------------
% Outputs:
% structArray/structure array: The input structure array but with the git
% state added to the .git field of each element
%--------------------------------------------------------------------------
% Variable Inputs:
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 1, 2014
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'gitState', 'struct', []};
defaults(end+1,:) = {'verbose', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabFunctions:invalidArguments', 'A structure array is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
% Generate current git state
%--------------------------------------------------------------------------
if isempty(parameters.gitState)
    gitState = GitState();
else
    gitState = parameters.gitState;
end

%--------------------------------------------------------------------------
% Check for uncommitted changes
%--------------------------------------------------------------------------
for i=1:length(gitState)
    if ~gitState(i).clean
        warning('matlabfunctions:GitUncommitedChanges', ...
            ['The ' gitState(i).repoName ' repository has uncommitted changes!']);
    end
end

%--------------------------------------------------------------------------
% Add git state to each element in structure array
%--------------------------------------------------------------------------
for i=1:length(structArray)
    structArray(i).git = gitState;
end

