function [axesHandles] = PlotCoverageVector(varargin)
%--------------------------------------------------------------------------
% [axesHandles] = PlotCoverageVector(coverVectorStructure, varargin)
% This function generates a figure from a coverage vector structure
%--------------------------------------------------------------------------
% Necessary Inputs
% coverageVectorStructure: A structure containing counts for top and bottom
% strand of a set of sequences. Multiple reference sequences can be
% includes as an array of coverage vector structures. 
%
%--------------------------------------------------------------------------
% Outputs
% axesHandles/array. A 2x1 array of the handles to the axes. The first is the
% coverage vector and the second is to the labels. 
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% February 9, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
requiredFields = {'refName', 'fileName', 'top', 'bottom'};
%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = true;
YMode = 'log';
colorStyle = {'g', 'r'};
logBase = 2;
continuous = true;
readLength = 50;
countAllBases = false;
titleString = '';
indices = [];
weight = 1;
figHandle = [];
filterLength = [];
toPlot = {'top', 'bottom'};
plotSign = {1, -1};
genomicFeatures = [];
axesHandles = -ones(1,2);

%--------------------------------------------------------------------------
% Confirm cvStruct input
%--------------------------------------------------------------------------
if ~strcmp(class(varargin{1}), 'struct')
    error(['The first argument to ''' mfilename ''' must be a cvStruct']);
end
foundFields = fields(varargin{1});
for i=1:length(requiredFields)
    if ~strcmp(requiredFields{i}, foundFields)
        error(['A valid cvStruct must contain ''' requiredFields{i} ''' ' ]);
    end
end
cvStruct = varargin{1};
varargin = varargin(2:end);

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
presentationPlotArgIn = {};
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'verbose'
                verbose = CheckParameter(parameterValue,'boolean','verbose');
            case 'YMode'
                YMode = CheckList(parameterValue, {'log', 'linear'}, 'YMode');
            case 'countAllBases'
                countAllBases = CheckParameter(parameterValue,'boolean','countAllBases');
            case 'titleString'
                titleString = CheckParameter(parameterValue,'string','titleString');
            case 'indices'
                indices = CheckParameter(parameterValue,'array','indices');
            case 'logBase'
                logBase = CheckParameter(parameterValue,'positive','logBase');
            case 'weight'
                weight = CheckParameter(parmeterValue, 'positive', 'weight');
            case 'figHandle'
                figHandle = CheckParameter(parameterValue, 'handle', 'figHandle');
            case 'filterLength'
                filterLength = CheckParameter(parameterValue, 'positive', 'filterLength');
            case 'toPlot'
                toPlot = CheckParameter(parameterValue, 'cell', 'toPlot');
            case 'plotSign'
                plotSign = CheckParameter(parameterValue, 'cell', 'plotSign');
            case 'colorStyle'
                colorStyle = CheckParameter(parameterValue, 'cell', 'colorStyle');
            case 'genomicFeatures'
                genomicFeatures = CheckParameter(parameterValue, 'struct', 'genomicFeatures');
            case 'axesHandles'
                axesHandles = CheckParameter(parameterValue, 'array', 'axesHandles');
            otherwise
                presentationPlotArgIn(end+1:end+2) = varargin( (parameterIndex*2 - 1):(parameterIndex*2) );
                %error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Create Figure(s)
%--------------------------------------------------------------------------
if isempty(figHandle)
    figHandle = figure('Name', ['Coverage Vector for ' cvStruct.fileName ' reference ' cvStruct.refName]);
end


%--------------------------------------------------------------------------
% Plot top and bottom strands
%--------------------------------------------------------------------------
figure(figHandle);
if ~isempty(genomicFeatures)
    if axesHandles(1) < 0
        axesHandles(1) = subplot(5,1, [1 2 3 4]);
    else
        axes(axesHandles(1));
    end
end

topStrand = double(cvStruct.top)/weight;
bottomStrand = double(cvStruct.bottom)/weight;

chromPos = 1:length(cvStruct.top);

if countAllBases
    kernel = ones(1,readLength);
    topStrand = conv(topStrand, kernel);
    topStrand = topStrand(1:(end-readLength+1));
    bottomStrand = conv(bottomStrand, kernel);
    bottomStrand = bottomStrand(1:(end-readLength+1));
end

if ~isempty(indices)
    topStrand = topStrand(indices(1):indices(2));
    bottomStrand = bottomStrand(indices(1):indices(2));
    chromPos = chromPos(indices(1):indices(2));
end

if ~isempty(filterLength)
    topStrand = filter(ones(1, filterLength), filterLength, topStrand);
    bottomStrand = filter(ones(1, filterLength), filterLength, bottomStrand);
end

switch YMode
    case 'log'
        topStrand = 1/log(logBase)*log(topStrand);
        bottomStrand = 1/log(logBase)*log(bottomStrand);

        if continuous
            topStrand(isinf(topStrand)) = 0;
            bottomStrand(isinf(bottomStrand)) = 0;
        end
        
        yLabel = ['Counts (log_{' num2str(logBase) '})'];

    case 'linear'
        yLabel = 'Counts';
        xlabel('Position (bp)');

end

for i=1:length(toPlot)
    switch toPlot{i}
        case 'top'
            plot(chromPos, plotSign{i}*topStrand, colorStyle{i}); hold on;
        case 'bottom'
            plot(chromPos, plotSign{i}*bottomStrand, colorStyle{i}); hold on;
    end
end
ylabel(yLabel);

%--------------------------------------------------------------------------
% Clean up plots
%--------------------------------------------------------------------------
if ~isempty(indices)
    xlim(indices);
else
    xlim([1 length(cvStruct.top)]);
end
title(titleString);

%--------------------------------------------------------------------------
% Plot Genomic context
%--------------------------------------------------------------------------
if ~isempty(genomicFeatures)
    if axesHandles(1) < 0
        axesHandles(1) = subplot(5,1, [1 2 3 4]);
    else
        axes(axesHandles(1));
    end
    set(gca, 'XTick', []);
    
    if axesHandles(2) < 0
        axesHandles(2) =  subplot(5,1,5);
    else
        axes(axesHandles(2));
    end
    
    featureInds = [genomicFeatures.indices];
    startInds = featureInds(1:2:end);
    finishInds = featureInds(2:2:end);
    case1Inds = startInds >= indices(1) & startInds <= indices(2);
    case2Inds = finishInds >= indices(1) & finishInds <= indices(2);
    case3Inds = startInds < indices(1) & finishInds > indices(2);
    
    featuresToDisplay = find(case1Inds | case2Inds | case3Inds);
    for i=1:length(featuresToDisplay)
        localInds = genomicFeatures(featuresToDisplay(i)).indices;
        if genomicFeatures(featuresToDisplay(i)).isComplement
            plot(localInds, zeros(1,2), 'b-'); hold on;
            plot(localInds(1), 0,  'b<'); 
        else
            plot(localInds, zeros(1,2), 'r-'); hold on;
            plot(localInds(2), 0, 'r>');
        end
        textLocation = mean(localInds);
        if textLocation < indices(1)
            textLocation = indices(1) - 0.05*diff(indices);
        end
        if textLocation > indices(2)
            textLocation = indices(2) + 0.05*diff(indices);
        end
        text(textLocation, .5, genomicFeatures(featuresToDisplay(i)).name);
    end
    set(gca, 'YTick', []);
    xlim(indices);
end
xlabel('Genomic Position (bp)');
linkaxes(axesHandles, 'x');
PresentationPlot(presentationPlotArgIn);
