function [missedRate,gainedRate,detectionRate] = GetBitFlipRates(libCodes,spotMatrices,varargin)
%  [missedRate,gainedRate,detectionRate] = GetBitFlipRates(libCodes,spotMatrices)

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'flipBits', 'boolean', false};

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


if parameters.flipBits
    libCodes = fliplr(libCodes);
end

%% Main Function
[numGenes,numHybes] = size(libCodes); 
minPhotsPerStain = 1;
if length(minPhotsPerStain) == 1;
minPhotsPerStain = repmat(minPhotsPerStain,1,numHybes);
end



correctLocsPerHybe = NaN*zeros(numHybes,numGenes); 
allMissedBits = NaN*zeros(numHybes,numGenes); 
allGainedBits = NaN*zeros(numHybes,numGenes); 
expectedBits = zeros(numHybes,1); 
cntsInHybe = zeros(numHybes,1);
detectRate  = NaN*zeros(numHybes,numGenes); 
detectRateN  = NaN*zeros(numHybes,numGenes); 
totSpots = 0;
totLocs = cell(numGenes,1); 
i = 0;
for n=1:numGenes % i = 3
    i=i+1;
    totLocs{i} = libCodes(n,:)*spotMatrices{n};
    numSpots = size(spotMatrices{n},2);
    correctLocs = repmat(libCodes(n,:)',1,numSpots) .* spotMatrices{n};  % numHybes by numSpots 
    locsTrueBits = correctLocs(libCodes(n,:),:);
    expectedBits = expectedBits + libCodes(n,:)'*numSpots; 
    
   % Find missed bits from corrected words
    missedBits = sum(correctLocs < repmat(minPhotsPerStain',1,numSpots),2);
    missedBits( ~libCodes(n,:),:) = 0;
    allMissedBits(:,i) = missedBits; 
    
    % find gained bits from corrected words
    gainedBits = repmat(~libCodes(n,:)',1,numSpots) .* spotMatrices{n};
    allGainedBits(:,i) = sum(gainedBits >= repmat(minPhotsPerStain',1,numSpots),2);
    correctLocsPerHybe(:,i) = mean(correctLocs, 2);
    
    onBits = find(libCodes(i,:));
    for j=onBits % i = 3; j = 1
        detectRate(j,i) = sum(spotMatrices{i}(j,:)>= minPhotsPerStain(j))/length(spotMatrices{i}(j,:));
        detectRateN(j,i) = numSpots*detectRate(j,i); 
        cntsInHybe(j) = cntsInHybe(j) + numSpots;
    end
    totSpots = totSpots + numSpots;
end
missedRate = nansum(allMissedBits,2)./expectedBits;
gainedRate = nansum(allGainedBits,2)./expectedBits;
detectionRate = nansum(detectRateN,2)./cntsInHybe;

