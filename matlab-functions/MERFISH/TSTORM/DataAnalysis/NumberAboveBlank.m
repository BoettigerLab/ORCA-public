function objectiveValue = NumberAboveBlank(words, wordIDsNonBlank, wordIDsBlank)
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

numNonBlank = arrayfun(@(x) sum(words==x), wordIDsNonBlank);
numBlank = arrayfun(@(x) sum(words==x), wordIDsBlank);

objectiveValue = sum(numNonBlank > max(numBlank));