function objectiveValue = HybeAccuracy(words,brightnessPerSpot, decCodes, allCorCodes,varargin)
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
defaults(end+1,:) = {'allCorCodes', 'cell', {}};
defaults(end+1,:) = {'biCodes', 'array', []};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function
[~,~,spotMatrix] = DecodeConvImage(words,decCodes,...
                                'codeType','decimal',...
                                'correctErrors',true,...
                                'allCorCodes',allCorCodes,...
                                'brightnessPerSpot',brightnessPerSpot);

if isempty(parameters.biCodes)
    [~,numHybes] = size(brightnessPerSpot);
    parameters.biCodes = de2bi(decCodes,numHybes); 
end

[perHybeMissRate,perHybeGainRate] = GetBitFlipRates(parameters.biCodes,spotMatrix);
objectiveValue = 1 - mean(perHybeMissRate) - mean(perHybeGainRate);


