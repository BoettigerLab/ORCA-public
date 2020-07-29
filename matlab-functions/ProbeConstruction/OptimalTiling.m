function [probeInds, parameters] = OptimalTiling(startPos, separations, varargin)
% ------------------------------------------------------------------------
% [probeInds, parameters] = OptimalTiling(startPos, separations, varargin)
% This function returns the indices of the set of probes that contains the
% most members but with no members that overlap. The starting position of
% the probes is provided with startPos and the distance between the
% starting position of each probe the next probe is provided in separations. 
%--------------------------------------------------------------------------
% Necessary Inputs
% startPos -- An array of the starting positions of all possible probes. 
% separations -- An array of the distance between the start of each probe and
%   the possible start of the downstream probe. If the length of separations
%   is not equal to the length of startPos, all separations are assumed to
%   be set by the first element of separations. 
%--------------------------------------------------------------------------
% Outputs
% probeInds -- The index of all probes that should be kept. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 16, 2015
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for parsing file names
defaults(end+1,:) = {'verbose', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2 || any(startPos < 0) ...
        || any(separations < 1)
    error('matlabFunctions:invalidArguments', ...
        ['A valid list of starting positions and separations must be provided.']);
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Coerce dimensions of input
% -------------------------------------------------------------------------
startPos = startPos(:)'; % Force row vector
separations = separations(:)'; 
numPos = length(startPos);

% -------------------------------------------------------------------------
% Check dimensionality of separations
% -------------------------------------------------------------------------
if length(separations) ~= length(startPos)
    separations = separations(1)*ones(1, length(startPos));
end

% -------------------------------------------------------------------------
% Build distance matrix and graph
% -------------------------------------------------------------------------
% Define uniform separation
separations = repmat(separations, [numPos 1]);

% Compute distance matrix
D = repmat(startPos', [1 numPos]) - repmat(startPos, [numPos 1]);

% Define graph  
graph = D;
graph(D > separations) = -1; % Only connect nodes that have a positive distance and that are separated by more than the appropriate value of separations
graph(D <= separations) = 0; % All other nodes are not connected

% -------------------------------------------------------------------------
% Identify potential starting nodes
% -------------------------------------------------------------------------
startingInds = find(graph(:,1)==0);

% -------------------------------------------------------------------------
% Loop over starting nodes and find minimum (maximum) paths
% -------------------------------------------------------------------------
% Initialize registers
maxPathLength = 0;
probeInds = [];
% Loop over starting inds
for i=1:length(startingInds)
    % Find all distances and paths
    [dist, path] = graphshortestpath(sparse(graph'), startingInds(i), 'Directed', true, 'Method', 'Acyclic');
    
    % Find minimum path
    [localMaxPath, pathInd] = min(dist);
    localMaxPath = -localMaxPath + 1; % Change to positive path length
    
    % Compare to previous starting inds and update registers
    if localMaxPath > maxPathLength
        maxPathLength = localMaxPath;
        probeInds = path{pathInd};
        %probeInds = path;
    end
end
