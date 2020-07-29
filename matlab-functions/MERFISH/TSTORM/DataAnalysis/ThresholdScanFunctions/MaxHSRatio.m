function objectiveValue = MaxHSRatio(words,codeWords, geneInds, blankInds,varargin)
% ------------------------------------------------------------------------
% objectiveValue = NumberAboveBlank(words, wordIDsNonBlank, wordIDBlank)
% This function determines the number of species greater than the largest
% blank. 
%--------------------------------------------------------------------------
% Necessary Inputs
% words -- integer valued codewords
% codeWords -- logical version of codewords assigned to genes
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 'intCodewords' -- integer form of codewords assigned to genes
% (precomputed for speed). 
% 
% 
%--------------------------------------------------------------------------
% Alistair Boettiger, Jeffrey Moffitt
% November 26, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Allow logical output (faster) 
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'intCodewords', 'nonnegative', []};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function

intCodewords = parameters.intCodewords;

% precompute for speed
if isempty(intCodewords)
    intCodewords = bi2de(codeWords);
end
[numGenes,numHybes] = size(codeWords);

% -------------------------------------------------------------------------
% Compute integer codeword histogram
% -------------------------------------------------------------------------
[wordCounts,x] = hist(words, 1:2^numHybes);

% -------------------------------------------------------------------------
% Compute hamming sphere counts
% -------------------------------------------------------------------------
hammingSphereCounts = zeros(2,numGenes);
hammingSphereCounts(1,:) = wordCounts(intCodewords);
for j=1:numGenes
    surroundingWords = GenerateSurroundingCodewords(codeWords(j,:), 1,'logical',true);
    surroundingWordIntegers = bi2de(surroundingWords);
    hammingSphereCounts(2,j) = sum(wordCounts(surroundingWordIntegers));
end

geneHSRatio = hammingSphereCounts(1,geneInds)./hammingSphereCounts(2,geneInds);
blankHSRatio= hammingSphereCounts(1,blankInds)./hammingSphereCounts(2,blankInds);
objectiveValue = sum(geneHSRatio > max(blankHSRatio));
