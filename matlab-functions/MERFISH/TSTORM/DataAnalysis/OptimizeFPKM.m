function objectiveValue = OptimizeFPKM(words, geneID2fpkm)
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


[uniqueWords,wordCnts] = occurrences(words); 
geneFPKMs = arrayfun(@(x) geneID2fpkm(x),uniqueWords,'ErrorHandler',@ErrorFunZero);
realGeneIdx = geneFPKMs~=0;
% geneFPKMs = arrayfun(@(x) geneID2fpkm(x),uniqueWords','ErrorHandler',@ArrayErrors);
objectiveValue = corr(geneFPKMs(realGeneIdx)',wordCnts(realGeneIdx)');


function output = ArrayErrors(input)
output = 0;