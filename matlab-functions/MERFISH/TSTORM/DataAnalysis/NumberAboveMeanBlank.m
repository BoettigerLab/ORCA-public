function objectiveValue = NumberAboveMeanBlank(words, wordIDsNonBlank, wordIDsBlank,varargin)
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

% general parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'numStd', 'float', 1};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


%%

numNonBlank = arrayfun(@(x) sum(words==x), wordIDsNonBlank);
numBlank = arrayfun(@(x) sum(words==x), wordIDsBlank);

maxBlank = mean(numBlank) + parameters.numStd*std(numBlank);
objectiveValue = sum(numNonBlank > maxBlank);