function clean = GitIsClean(varargin)
%--------------------------------------------------------------------------
% clean = GitIsClean(varargin)
% This function returns a boolean which indicates whether all registered 
% git repositories are clean or not.
%--------------------------------------------------------------------------
% Outputs:
% clean/boolean: A boolean indicating whether all of the git repositories
% are clean, i.e. do not have uncommited changes. 
%--------------------------------------------------------------------------
% Variable Inputs:
% 'verbose'/boolean (True): Controls the verbosity of this function
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% March 30, 2014
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Default variables
%--------------------------------------------------------------------------
verbose = true;

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global gitRepositoryPaths;
if isempty(gitRepositoryPaths)
    error('No registered git repositories');
end

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
% Get Git State
%--------------------------------------------------------------------------
gitState = GitState('fields', {'repoName', 'clean'});
clean = all([gitState.clean]);

%--------------------------------------------------------------------------
% Display clean status
%--------------------------------------------------------------------------
if verbose
    if clean
        display('All repositories are clean');
    else
        display('The following repositories contain uncommitted changes');
        uncleanInds = find(~[gitState.clean]);
        for j=1:length(uncleanInds)
            display(['   ' gitState(uncleanInds(j)).repoName]);
        end
    end
end
