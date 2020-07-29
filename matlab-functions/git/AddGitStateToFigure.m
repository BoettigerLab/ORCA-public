function [figHandle, gitState] = AddGitStateToFigure(figHandle, gitState, varargin)
% ------------------------------------------------------------------------
% function [figHandle, gitState] = AddGitStateToFigure(figHandle, gitState, varargin)
% This function adds a git state structure to the UserData field of a figure.  
%--------------------------------------------------------------------------
% Necessary Inputs
% figHandle/handle: Handle to a valid figure. 
% gitState/A git state structure. See GitState. If not provided, then the
%   state is automatically determined by GitState. 
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
% April 8, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Default variables
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Parse necessary input
% ------------------------------------------------------------------------
if nargin < 1
    error('No figure handle provided');
elseif nargin < 2
    gitState = GitState();
end
if nargin >= 2
    if ~isstruct(gitState)
        temp = varargin;
        varargin{1} = gitState;
        varargin{2:length(temp)} = temp(:);
    end
end
% ------------------------------------------------------------------------
% Parse optional input
% ------------------------------------------------------------------------
if ~isempty(varargin)
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename])
        end
    end
end

% ------------------------------------------------------------------------
% Check validity of figure handle
% ------------------------------------------------------------------------
if ~(ishghandle(figHandle) && strcmp(get(figHandle, 'type'), 'figure'))
    error('matlabFunctions:invalidHandle', 'Provided figure handle is not valid');
end

% ------------------------------------------------------------------------
% Access existing userdata and add git state
% ------------------------------------------------------------------------
userData = get(figHandle, 'UserData');
if ~isempty(userData)
    if iscell(userData)
        isGitState = false(1, length(userData));
        for i=1:length(userData)
            if isstruct(userData{i}) && isfield(userData{i}, 'repoName')
                userData{i} = gitState;
                isGitState(i) = true;
            end
        end
        if ~any(isGitState)
            userData{end+1} = gitState;
        end
    else
        temp = userData;
        if isstruct(temp) && isfield(temp, 'repoName')
            userData = gitState;
        else
            userData = {};
            userData{1} = temp;
            userData{2} = gitState;
        end
    end
else
    userData = gitState;
end

% ------------------------------------------------------------------------
% Add back to figure
% ------------------------------------------------------------------------
set(figHandle, 'UserData', userData);
