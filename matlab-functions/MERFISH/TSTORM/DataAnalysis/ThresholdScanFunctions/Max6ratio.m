function objectiveValue = Max6ratio(words,numHybes)
% ------------------------------------------------------------------------
% objectiveValue = Max6ratio(words)
% this function mesures the number 6s relative to 5s 6s and 7s. 
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% February 12, 2016
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

biWords = de2bi(words,numHybes);
num5s = sum(sum(biWords,2) == 5);
num6s = sum(sum(biWords,2) == 6);
num7s = sum(sum(biWords,2) == 7);

objectiveValue = num6s/ (num5s+num6s+num7s);