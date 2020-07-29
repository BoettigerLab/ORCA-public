function [corrM,bitPairs,parameters] = BitCorrelation(Mcents,varargin)
% compute pair-wise correlation between bits

%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger.alistair@gmail.com
% April 16, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'maxCorr', 'nonnegative', .3}; % threshold to ID super-correlated bits


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'A MList is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


%% Main Function

[~,numHybes] = size(Mcents); 

corrM = zeros(numHybes); 
for i=1:numHybes 
    for j=1:numHybes  % i = 1; j = 2
        corrM(i,j) = sum((Mcents(:,i)+Mcents(:,j))==2)/sum((Mcents(:,i)+Mcents(:,j))>0);
    end
end

% Plotting command
imagesc(corrM);
colorbar; caxis([0,.5]);
xlabel('bit number'); ylabel('bit number');
title('frequency of pairwise bit co-localization');
PresentationPlot;

% return unique bitpairs thar are highly correlated
corrM2 = corrM.*tril(ones(numHybes),-1);
[x,y] = ind2sub([numHybes,numHybes], find(corrM2> parameters.maxCorr & corrM2<1));
bitPairs = [x,y];
