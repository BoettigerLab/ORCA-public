function [polyCistrons, connectedMat, parameters] = IdentifyPolycistrons(coverageVector, genomicFeatures, varargin)
% ------------------------------------------------------------------------
% [polyCistrons, connectedMat, parameters] = IdentifyPolycistrons(coverageVector, features, varargin)
% This function identifies polycistronic messages from the relationship 
% between normalized counts between adjacent genomicFeatures as well as
% the geometry of those features, contained in features. 
%--------------------------------------------------------------------------
% Necessary Inputs
% coverageVector/structure. A structure with the following fields:
%   --top: The counts per bp on the top strand.
%   --bottom: The counts per bp on the bottom strand. 
% genomicFeatures/structure. See LoadGenbank for an example. This structure
%   must have the following fields:
%   --isComplement. A boolean that determines the strand of the feature.
%   --indices. A list of the initial and final bp of the feature
%   --name. The name of the feature. 
%
%--------------------------------------------------------------------------
% Outputs
% polyCistrons/structure array. Each element contains the following fields
%  --name: The name of the polycistron
%  --length: The number of genomic features in the polycistron
%  --isComplement: The strand of the polycistron
%  --names: The names of all features within the cistron. 
%  --inds: The starting and stopping inds of the polycistron on the genome
%  --featureInds: The indices of the features in genomicFeatures contained
%  in the polycistron.
%
% connectedMat: NxN boolean matrix. The ijth element is true if the ith element
%   of genomicFeatures is in a polycistron with the jth element.  
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 29, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'adjacentLogRatio', 'positive', 2};
defaults(end+1,:) = {'intergenicLogRatio', 'positive', 2};
defaults(end+1,:) = {'intergenicPad', 'positive', 100};
defaults(end+1,:) = {'verbose', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabFunctions:invalidArguments', 'Both a list of counts and a genomic features structure are required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Determine feature counts
% -------------------------------------------------------------------------
countStructure = GenerateFeatureCounts(coverageVector, genomicFeatures, 'verbose', parameters.verbose);
normCountsPerFeature = [countStructure.senseMean]; % The mean counts on the sense strand

% -------------------------------------------------------------------------
% Determine boundaries based on counts ratios
% -------------------------------------------------------------------------
if parameters.adjacentLogRatio ~= 0    
    genomicRatio = log2(normCountsPerFeature./circshift(normCountsPerFeature, [0 -1])) > parameters.adjacentLogRatio | ...
        log2(normCountsPerFeature./circshift(normCountsPerFeature, [0 -1])) < -parameters.adjacentLogRatio;
else
    genomicRatio = false(1, length(genomicFeatures));
end

% -------------------------------------------------------------------------
% Determine boundaries based on change in strand
% -------------------------------------------------------------------------
directionSwitch = xor([genomicFeatures.isComplement], ...
    circshift([genomicFeatures.isComplement], [0 -1]));

% -------------------------------------------------------------------------
% Compute intergenic counts
% -------------------------------------------------------------------------
intergenicCounts = Inf(1, length(genomicFeatures));
indices = [genomicFeatures.indices];
for i=1:(length(indices)/2)
    if i~=(length(indices)/2)
        intergenicInds = indices(2*i):indices(2*i+1);
    else
        intergenicInds = [indices(end):length(coverageVector.top) ...
            1:indices(1)];
    end
    if length(intergenicInds) > parameters.intergenicPad
        if genomicFeatures(i).isComplement
            intergenicCounts(i) = mean(coverageVector.bottom(intergenicInds));
        else
            intergenicCounts(i) = mean(coverageVector.top(intergenicInds));
        end
    end
end

% -------------------------------------------------------------------------
% Determine intergenic ratios
% -------------------------------------------------------------------------
if parameters.intergenicLogRatio ~=0
    intergenicRatio = log2(intergenicCounts./normCountsPerFeature) < parameters.intergenicLogRatio;
else
    intergenicRatio = false(1, length(genomicFeatures));
end

% -------------------------------------------------------------------------
% Combine boundaries
% -------------------------------------------------------------------------
boundaries = genomicRatio | directionSwitch | intergenicRatio;
switchPositions = find(boundaries);

% -------------------------------------------------------------------------
% Build polycistron structure
% -------------------------------------------------------------------------
count = 1;
for i=0:(length(switchPositions)-1)
    % ---------------------------------------------------------------------
    % Determine putative inds for combined features
    % ---------------------------------------------------------------------
    if i==0
        featureInds = 1:switchPositions(1);
    else
        featureInds = (switchPositions(i)+1):switchPositions(i+1);
    end   
    
    % ---------------------------------------------------------------------
    % Build polycistron element
    % ---------------------------------------------------------------------
    polyCistrons(count).featureInds = featureInds;
    polyCistrons(count).length = length(featureInds);
    polyCistrons(count).isComplement = genomicFeatures(featureInds(1)).isComplement;
    polyCistrons(count).names = {genomicFeatures(featureInds).name};
    polyCistrons(count).indices = [genomicFeatures(featureInds(1)).indices(1), ...
        genomicFeatures(featureInds(end)).indices(2)];
    name = genomicFeatures(featureInds(1)).name;
    for j=2:length(featureInds)
        name = [name '-' genomicFeatures(featureInds(j)).name];
    end
    polyCistrons(count).name = name;

    % ---------------------------------------------------------------------
    % Display info
    % ---------------------------------------------------------------------
    if parameters.verbose
        display(['Polycistron ' num2str(i) ': ' num2str(polyCistrons(count).inds(1)) '-' ...
            num2str(polyCistrons(count).inds(2)) ]);
        for j=1:length(polyCistrons(count).featureInds)
            display(genomicFeatures(polyCistrons(count).featureInds(j)).name);
        end
    end
    count = count + 1;
end

% -------------------------------------------------------------------------
% Compute connected matrix
% -------------------------------------------------------------------------
connectedMat = false(length(genomicFeatures));

for i=1:length(polyCistrons)
    featureInds = polyCistrons(i).featureInds;
    connectedMat(featureInds, circshift(featureInds, [0 -(i-1)])) = true(length(featureInds));
end
    
    
    