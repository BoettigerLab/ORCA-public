function [dataFun,dataCI] = FunWithCI(data,fun,varargin)
% MedianWithCI(data)
% Inputs
% data -  to take median of (vector, matrix, or cell array)
% fun - inline function to analyze.   fun = @(x) median(x,2);
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
% A little more parameter parsing

%%
lowLim = (100 -parameters.cI)/2;
highLim = 100 - lowLim;

% fun = @(x) median(x,2);

% vector
if ~iscell(data) && min(size(data)) == 1
    data = data(~isnan(data));
    numDataPts = sum(~isnan(data));
    randSample = randi(numDataPts,parameters.iters,numDataPts);
    medianBoots =  sort( fun( data(randSample) ) );
    dataCI = medianBoots(round([lowLim,highLim]*parameters.iters/100))';
    dataFun = fun(data);
end

% matrix
if ~iscell(data) && min(size(data)) > 1
    numSamples = size(data,1);
    dataFun = zeros(numSamples,1);
    dataCI = zeros(numSamples,2);
    for i=1:numSamples
        dat = data(i,:); % could be done with matrix algebra instead of loops
        dat = dat(~isnan(dat));
        numDataPts = sum(~isnan(dat));
        randSample = randi(numDataPts,parameters.iters,numDataPts);
        medianBoots =  sort( fun( data(randSample) ) );
        dataCI(i,:) = medianBoots(round([lowLim,highLim]*parameters.iters/100))';
        dataFun(i) = fun(dat);
    end
end

% cell
if iscell(data)
    dataFun = NaN(length(data),1);
    dataCI = NaN(length(data),2);
    for i=1:length(data)
        dat = data{i};
        dat = dat(~isnan(dat));
        numDataPts = sum(~isnan(dat));
        if numDataPts > 0
            randSample = randi(numDataPts,parameters.iters,numDataPts);
            medianBoots =  sort( fun( data(randSample) ) );
            dataCI(i,:) = medianBoots(round([lowLim,highLim]*parameters.iters/100))';
            dataFun(i) = fun(dat);
        end
    end
end



