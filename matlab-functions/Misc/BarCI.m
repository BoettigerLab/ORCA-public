function BarCI(data,varargin)
% BarCI(data)
% Inputs
% data -  to take median of (vector, matrix, or cell array)
% 
% returns the median value of data as a first output 
% and the 95% confidence intervals by bootstrapping as the second output. 
%
% MedianWithCI(data,'cI',99) 
% returns the 99% confifidence intervals.  
%
% MedianWithCI(data,'iters',1E5)
% changes the number of iterations.  
% 
% works with:
% vector input - 1 x n observations
% matrix input - m datasets x n observations
% cell input - m datasets with variable number of observations. 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'cI', 'positive', 95};
defaults(end+1,:) = {'iters', 'positive', 1000};
defaults(end+1,:) = {'x', 'array', []};
defaults(end+1,:) = {'faceColor', 'colormap', [.5 .5 .5]};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);

% -------------------------------------------------------------------------

[m,cI] = MedianWithCI(data,'cI',parameters.cI,'iters',parameters.iters);
if isempty(parameters.x)
    bar(m,1,'FaceColor',parameters.faceColor); hold on;
    ploterr(1:length(m),m,[],{cI(:,1),cI(:,2)},'.','color','k');
else
    bar(parameters.x,m,1,'FaceColor',parameters.faceColor); hold on;
    ploterr(parameters.x,m,[],{cI(:,1),cI(:,2)},'.','color','k');
end

