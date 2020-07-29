function [wordsDetected,brightnessPerSpot] = AssignConvSpotsToCentroids(flists,mRNAcents,varargin)
% [wordsDetected,datPos] = AssignLocsToMRNAcentroids(flists,mRNAcents,varargin)

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'minPhotsPerStain', 'nonnegative', 1};
defaults(end+1,:) = {'maxDtoCentroid', 'positive', 1};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires flists, mRNAcents');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function

minPhotsPerStain = parameters.minPhotsPerStain;
numHybes = length(flists); 

% minPhotsPerStain = [5 ,3];
if length(minPhotsPerStain) < numHybes
    minPhotsPerStain = [minPhotsPerStain, repmat(minPhotsPerStain(end),1, numHybes - length(minPhotsPerStain)) ];
end

% Assign localizations to mRNA centroids
wordsDetected = false(length(mRNAcents),numHybes);
brightnessPerSpot = zeros(length(mRNAcents),numHybes);
for h=1:numHybes

    [idx,di] = knnsearch([flists{h}.xc,flists{h}.yc],mRNAcents,'K',1);
    try
        brightnessPerSpot(:,h) = flists{h}.a(idx);
        brightnessPerSpot(di >= parameters.maxDtoCentroid,h) = 0; 
    wordsDetected(:,h) = brightnessPerSpot(:,h) > minPhotsPerStain(h); 
    catch er
        disp(er.message)
    end
    
end





