function [axesHandle, objHandles, boxQuantileValues, whiskersQuantileValues] = PlotBoxAndWhiskers(data, xPos, axesHandle, varargin)
% ------------------------------------------------------------------------
% [axesHandle, objHandles] = PlotBoxAndWhiskers(data, xPos, axesHandle, varargin)
% This creates a box with whiskers for the provided data and plots them on
% the provided axesHandle. 
%--------------------------------------------------------------------------
% Necessary Inputs
% data/An array of data points.
% xPos/The position of box and whiskers to be plotted. 
% axesHandle/A handle to an axes object. The box and whiskers will be
% plotted there. If this is not provided, then the axesHandle is the
% current axis. 
%--------------------------------------------------------------------------
% Outputs
% axesHandle/The handle to object in which the box and whiskers were
% plotted
% objHandles/An array of handles to all of the objects that represent the
% box and whiskers. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% October 28, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'boxQuantileProbabilities','array', [0.25 0.75]};
defaults(end+1,:) = {'whiskerQuantileProbabilities','array', [0.09 0.91]};
defaults(end+1,:) = {'centralLine', {'mean', 'median', 'none'}, 'mean'};
defaults(end+1,:) = {'boxColor', 'string', 'k'};
defaults(end+1,:) = {'whiskersColor', 'string', 'k'};
defaults(end+1,:) = {'centralLineColor', 'string', 'k'};
defaults(end+1,:) = {'boxWidth', 'nonnegative', 0.5};
defaults(end+1,:) = {'whiskerEndWidth', 'nonnegative', 0.25};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2 || length(xPos) > 1
    error('matlabFunctions:invalidArguments', 'Invalid values for data and xPos.');
end

% -------------------------------------------------------------------------
% Identify axes handle 
% -------------------------------------------------------------------------
if nargin < 3 
    if isempty(varargin)
        axesHandle = gca;
    else
        axesHandle = varargin{1};
        varargin = varargin(2:end);
        if ~ishghandle(axesHandle)
            error('matlabFunctions:invalidArguments', 'Invalid axes handle.');
        end
    end
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Remove nan
% -------------------------------------------------------------------------
data = data(~isnan(data));

% -------------------------------------------------------------------------
% Compute quantiles, SEM, median, and mean
% -------------------------------------------------------------------------
boxQuantileValues = quantile(data, parameters.boxQuantileProbabilities);
whiskersQuantileValues = quantile(data, parameters.whiskerQuantileProbabilities);

meanValue = mean(data);
medianValue = median(data);
SEM = std(data)/sqrt(length(data));

% -------------------------------------------------------------------------
% Create box and whiskers
% -------------------------------------------------------------------------
hold on;
% Quantile Box
objHandles(1) = rectangle('Position', ...
    [-parameters.boxWidth/2 + xPos, ....
    boxQuantileValues(1), ...
    parameters.boxWidth, ...
    boxQuantileValues(2) - boxQuantileValues(1)], ...
    'EdgeColor', parameters.boxColor);

% Whiskers
objHandles(2) = plot(ones(1,2)*xPos, whiskersQuantileValues, ...
    parameters.whiskersColor);

% Whiskers Caps
if parameters.whiskerEndWidth > 0
    objHandles(3) = plot([-parameters.whiskerEndWidth/2 parameters.whiskerEndWidth/2]+xPos, ...
        whiskersQuantileValues(1)*ones(1,2), parameters.whiskersColor);
    objHandles(4) = plot([-parameters.whiskerEndWidth/2 parameters.whiskerEndWidth/2]+xPos, ...
        whiskersQuantileValues(2)*ones(1,2), parameters.whiskersColor);
end

% Middle line
switch parameters.centralLine
    case 'mean'
        objHandles(5) = plot([-parameters.boxWidth/2 parameters.boxWidth/2]+ xPos, ...
            ones(1,2)*meanValue, parameters.centralLineColor);
    case 'median'
        objHandles(5) = plot([-parameters.boxWidth/2 parameters.boxWidth/2]+ xPos, ...
            ones(1,2)*medianValue, parameters.centralLineColor);
    otherwise
end

