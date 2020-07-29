function objectiveValue = Max4ratio(barcodes,numBits)
% ------------------------------------------------------------------------
% objectiveValue = Max4ratio(barcodes, numBits)
% This function measures the number of barcodes with a Hamming Weight of 4
% relative to the total number of barcodes with a HW of 3, 4, or 5. 
%--------------------------------------------------------------------------
% Necessary Inputs
% --words/An array of integer words (logical array) 
% --numBits/The number of bits in the binary barcodes. 
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% March 13, 2016
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

biWords = de2bi(barcodes,numBits);
num3s = sum(sum(biWords,2) == 3);
num4s = sum(sum(biWords,2) == 4);
num5s = sum(sum(biWords,2) == 5);

objectiveValue = num4s/ (num3s+num4s+num5s);