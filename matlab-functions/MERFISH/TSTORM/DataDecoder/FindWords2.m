function [wordsDetected,wordLocations,brightnessPerSpot] = FindWords2(flists,spotPositions,varargin)
% [wordsDetected,datPos] = AssignLocsToMRNAcentroids(flists,mRNAcents,varargin)

error('function does not work yet')

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
        
        
%         tempLists = [flists{[1:h-1,h+1:end]}];
%         molsPerHybe = cellfun(@(x) length(x.xc),flists([1:h-1,h+1:end]));
%         molIdx = cumsum(molsPerHybe)
%         [idx,di] = knnsearch([cat(1,tempLists.xc),cat(1,tempLists.yc)],[flists{h}.xc,flists{h}.yc],'K',numHybes-1);
%         hybeMat = idx;
%         hybeMat(di>= parameters.maxDtoCentroid) = NaN;
%         figure(3); clf; imagesc(hybeMat)
        
        notMapped = length(flists{h}.xc)+1;
        [idx,di] = knnsearch([flists{h}.xc,flists{h}.yc],spotPositions,'K',1);
        tooFar = di >= parameters.maxDtoCentroid;
        self = di==0;
        idx(tooFar | self) = notMapped;
        di(tooFar | self) = inf;
        inds = unique(idx); 
        
        hybeInds = notMapped*ones(length(spotPositions),1);
        for i=1:length(inds)
            allMatches = idx == inds(i);
            matchIndex = find(allMatches);
            [~,closestIdx] = min( di(allMatches) );
            hybeInds(matchIndex(closestIdx)) = inds(i);
        end
             flists{h}.xc = [flists{h}.xc; NaN];
             flists{h}.yc = [flists{h}.yc; NaN];
             
        xPerSpot(:,h) = flists{h}.xc(hybeInds);       
        yPerSpot(:,h) = flists{h}.yc(hybeInds);
        xPerSpot(1:length(flists{h}.xc(1:end-1)),h) = flists{h}.xc(1:end-1);
        yPerSpot(1:length(flists{h}.yc(1:end-1)),h) = flists{h}.yc(1:end-1);
    catch er
        disp(er.message)
    end
    
    figure(3); clf; imagesc(xPerSpot); pause(1); 
    
end
figure(3); clf; imagesc(xPerSpot);
meanX = nanmean(xPerSpot,2);
meanY = nanmean(yPerSpot,2);
wordLocations = [meanX,meanY];
[wordLocations,uniqueIdx] = unique(wordLocations,'rows','stable');
% [wordLocations,uniqueIdx] = unique(wordLocations,'rows');

wordsDetected = wordsDetected(uniqueIdx,:);
brightnessPerSpot = brightnessPerSpot(uniqueIdx,:); 




