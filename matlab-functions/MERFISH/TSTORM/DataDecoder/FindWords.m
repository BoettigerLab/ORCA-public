function [wordsDetected,wordLocations,brightnessPerSpot] = FindWords(flists,spotPositions,varargin)
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
% parameters = ParseVariableArguments([], defaults, mfilename);

%% Main Function


minPhotsPerStain = parameters.minPhotsPerStain;
numHybes = length(flists); 

% minPhotsPerStain = [5 ,3];
if length(minPhotsPerStain) < numHybes
    minPhotsPerStain = [minPhotsPerStain, repmat(minPhotsPerStain(end),1, numHybes - length(minPhotsPerStain)) ];
end

% Assign localizations to mRNA centroids
wordsDetected = false(length(spotPositions),numHybes);
brightnessPerSpot = NaN(length(spotPositions),numHybes);
xPerSpot = NaN(length(spotPositions),numHybes);
yPerSpot = NaN(length(spotPositions),numHybes);
for h=1:numHybes
    try
        if ~isempty(flists{h})
            [idx,di] = knnsearch([flists{h}.xc,flists{h}.yc],spotPositions,'K',1);
            brightnessPerSpot(:,h) = flists{h}.a(idx);
            brightnessPerSpot(di >= parameters.maxDtoCentroid,h) = 0; 
            validIdx = brightnessPerSpot(:,h) > minPhotsPerStain(h);
            wordsDetected(:,h) = validIdx; 
            xPerSpot(:,h) = flists{h}.xc(idx); 
            xPerSpot(~validIdx,h) = NaN; 
            yPerSpot(:,h) = flists{h}.yc(idx);
            yPerSpot(~validIdx,h) = NaN; 
        end
    catch er
        warning(er.message)
    end
    
end

meanX = nanmean(xPerSpot,2);
meanY = nanmean(yPerSpot,2);
wordLocations = [meanX,meanY];
[wordLocations,uniqueIdx] = unique(wordLocations,'rows','stable');
% [wordLocations,uniqueIdx] = unique(wordLocations,'rows');

wordsDetected = wordsDetected(uniqueIdx,:);
brightnessPerSpot = brightnessPerSpot(uniqueIdx,:); 




