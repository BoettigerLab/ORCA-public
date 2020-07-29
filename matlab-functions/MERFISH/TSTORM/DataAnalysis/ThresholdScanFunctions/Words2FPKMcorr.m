function objectiveValue = Words2FPKMcorr(words, decCodes, libExpect)
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


% % This is WAY slow!
% [uniqueWords,wordCnts] = occurrences(words); 
% geneFPKMs = arrayfun(@(x) geneID2fpkm(x),uniqueWords,'ErrorHandler',@ErrorFunZero);
% realGeneIdx = geneFPKMs~=0;
% objectiveValue = corr(geneFPKMs(realGeneIdx)',wordCnts(realGeneIdx)');


% disp('decoding...');

perCnts = DecodeConvImage(words,decCodes,'correctErrors',false,'codeType','decimal');
objectiveValue = LogCorr(libExpect,perCnts);



