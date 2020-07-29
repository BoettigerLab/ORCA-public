function objectiveValue = Ratio4To345andNumAboveBlank(words, numHybes, wordIDsNonBlank, wordIDsBlank)
% ------------------------------------------------------------------------
% objectiveValue = NumberAboveBlank(words, wordIDsNonBlank, wordIDBlank)
% This function determines the number of species greater than the largest
% blank. 
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Alistair Boettiger, Jeffrey Moffitt
% November 26, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------
binaryWords = de2bi(words, numHybes);

num3s = sum( sum(binaryWords,2) == 3);
num4s = sum( sum(binaryWords,2) == 4);
num5s = sum( sum(binaryWords,2) == 5);

numNonBlank = arrayfun(@(x) sum(words==x), wordIDsNonBlank);
numBlank = arrayfun(@(x) sum(words==x), wordIDsBlank);

fracAboveBlank = sum(numNonBlank > max(numBlank))/1000 ;
fracFours = num4s/(num5s + num4s + num3s);

objectiveValue = fracAboveBlank + fracFours;

