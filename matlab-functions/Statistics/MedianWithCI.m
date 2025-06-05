function [dataMed,dataCI] = MedianWithCI(data,varargin)
% MedianWithCI(data)
% Inputs
% data -  to take median of (vector, matrix, or cell array)
% 
% returns the median value of data as a first output 
% and the 95% confidence intervals by bootstrapping as the second output. 
%
% MedianWithCI(data,'cI',99) 
% returns the 99% condfidence intervals.  
% 1 - (.05)^(1/2) ~ 78% condfidence interval. Non-overlapping are <.05 different.   
%
% MedianWithCI(data,'iters',1E3)
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
defaults(end+1,:) = {'cI', 'positive', 100*(1 - (.05)^(1/2))};
defaults(end+1,:) = {'iters', 'positive', 1000};
defaults(end+1,:) = {'verbose', 'boolean', true};
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
dataCI = NaN;
dataMed = NaN;

% vector
if ~iscell(data) && min(size(data)) == 1
    data = data(~isnan(data));
    numDataPts = sum(~isnan(data));
    if numDataPts > 0
        randSample = randi(numDataPts,parameters.iters,numDataPts); % resample from data WITH replacement a total of "iters" times  
        medianBoots = sort(nanmedian( data(randSample),2));
        dataCI = medianBoots(round([lowLim,highLim]*parameters.iters/100))';
        dataMed = nanmedian(data);
    else
        dataCI = NaN;
        dataMed = NaN;
        if parameters.verbose
           warning('no non-NaN data to measure median'); 
        end
    end
end

% matrix
if ~iscell(data) && min(size(data)) > 1
    numSamples = size(data,1);
    dataMed = zeros(numSamples,1);
    dataCI = zeros(numSamples,2);
    for i=1:numSamples
        dat = data(i,:); % could be done with matrix algebra instead of loops
        dat = dat(~isnan(dat));
        if ~isempty(dat)
        numDataPts = sum(~isnan(dat));
        try
            randSample = randi(numDataPts,parameters.iters,numDataPts);
        catch er
            error(er.getReport)
        end
        medianBoots = sort(nanmedian( dat(randSample),2));
        dataCI(i,:) = medianBoots(round([lowLim,highLim]*parameters.iters/100))';
        dataMed(i) = nanmedian(dat);
        end
    end
end

% cell
if iscell(data)
    dataMed = NaN(length(data),1);
    dataCI = NaN(length(data),2);
    for i=1:length(data)
        dat = data{i};
        dat = dat(~isnan(dat));
        numDataPts = sum(~isnan(dat));
        if numDataPts > 0
            randSample = randi(numDataPts,parameters.iters,numDataPts);
            medianBoots = sort(nanmedian( dat(randSample),2));
            dataCI(i,:) = medianBoots(round([lowLim,highLim]*parameters.iters/100))';
            dataMed(i) = nanmedian(dat);
        end
    end
end



