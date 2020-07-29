function [localizationsPerCluster, parameters] = SimulateLocalizationsInCluster(varargin)
% ------------------------------------------------------------------------
% [localizationsPerCluster, parameters] = SimulateLocalizationsInCluster(varargin)
% This function simulates the number of localizations observed for a series
% of clustered fluorophores. For example, it can be used to simulate the binding of 
% multiple labeled primary or secondary probes to an mRNA. 
%--------------------------------------------------------------------------
% Necessary Inputs
%
%--------------------------------------------------------------------------
% Outputs
% localizationsPerCluster/An Nx1 array of the number of localizations per
% cluster
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% July 14, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'numClusters', 'positive', 1000};
defaults(end+1,:) = {'numLocalizationsPerDye', 'positive', 10};
defaults(end+1,:) = {'dyeLocalizationDistributionFunction', 'function', @geornd};
defaults(end+1,:) = {'bindingProbabilityPerNestedProbe', 'nonnegative', .9};
defaults(end+1,:) = {'numberOfNestedProbesPerProbe', 'positive', [1]};
defaults(end+1,:) = {'detectionOffset', 'nonnegative', 1};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse length of binding probability and number of dyes per probe
% parameters
% -------------------------------------------------------------------------
if length(parameters.bindingProbabilityPerNestedProbe) ~= length(parameters.numberOfNestedProbesPerProbe)
    error('matlabFunctions::invalidParameters', ...
        'Length of binding probabilities per probe and number of dyes per probe must be the same');
end

% -------------------------------------------------------------------------
% Determine total number of binding sites
% -------------------------------------------------------------------------
numNestedMolecules = cumprod(parameters.numberOfNestedProbesPerProbe);

% -------------------------------------------------------------------------
% Determine number of localizations per binding site
% -------------------------------------------------------------------------
localizationsPerProbe = parameters.dyeLocalizationDistributionFunction(...
    1/parameters.numLocalizationsPerDye, [numNestedMolecules(end) parameters.numClusters]);

% -------------------------------------------------------------------------
% Determine Binding Probability Masks
% -------------------------------------------------------------------------
% Allocate memory
isBound = zeros([numNestedMolecules(end) parameters.numClusters length(parameters.numberOfNestedProbesPerProbe)]);

% Build a binding probability mask for each set of nested probes
for i=1:(length(parameters.numberOfNestedProbesPerProbe)-1)
    isBound(:,:,i) = repmat( ...
        rand([numNestedMolecules(i) parameters.numClusters]) <= parameters.bindingProbabilityPerNestedProbe(i), ...
        [prod(parameters.numberOfNestedProbesPerProbe((i+1):end)) 1]);
end
isBound(:,:,end) = rand([numNestedMolecules(end) parameters.numClusters]) <= parameters.bindingProbabilityPerNestedProbe(end);

% Build final mask: Require that for a dye to contribute, all of its nested
% probes must be bound
isBound = all(isBound, 3);

% -------------------------------------------------------------------------
% Determine cluster localizations
% -------------------------------------------------------------------------
localizationsPerProbe = localizationsPerProbe.*isBound;
localizationsPerCluster = sum(localizationsPerProbe, 1);

% -------------------------------------------------------------------------
% Return only clusters that pass the detection threshold
% -------------------------------------------------------------------------
localizationsPerCluster = localizationsPerCluster(localizationsPerCluster>=parameters.detectionOffset);
