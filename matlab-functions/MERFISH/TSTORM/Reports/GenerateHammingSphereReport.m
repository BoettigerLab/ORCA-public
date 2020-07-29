function [reportStruct, parameters] = GenerateHammingSphereReport(words, exactMap, varargin)
% ------------------------------------------------------------------------
% [reportStruct, parameters] = GenerateHammingSphereReport(words, exactMap, varargin)
% This function generates the counts in different hamming spheres for each
% gene. The default is to calculate hamming sphere 0 and 1.  
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A structure array with an element for each word. See
%   CreateWordsStructure for information about fields. 
%--------------------------------------------------------------------------
% Outputs
% figHandles/Handles for the generated figures. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% November 30, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'hammingSphereReport', 'on'}; %
defaultReports(end+1,:) = {'confidenceRatioReport', 'on'};
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate','cell', defaultReports};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'subFolder', 'string', ''};
defaults(end+1,:) = {'maxHammingSphere', 'positive', 1};
defaults(end+1,:) = {'blankWordIdentifiers', 'cell', {'blank', 'notarget'}};
defaults(end+1,:) = {'colorMap', 'function', @jet};
defaults(end+1,:) = {'numHistogramBins', 'positive', 25};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 2 || ...
        ~isstruct(words) || ...
        isempty(intersect('intCodeword', fields(words))) || ...
        ~strcmp(class(exactMap), 'containers.Map') )
    error('matlabFunctions:invalidArguments', 'Invalid words structure or exactMap.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Define useful function for converting logical to string 
% -------------------------------------------------------------------------
removeWs = @(x) x(~isspace(x));

% -------------------------------------------------------------------------
% Extract properties of map
% -------------------------------------------------------------------------
codeWords = keys(exactMap);
geneNames = values(exactMap);

intCodewords = cellfun(@(x) bin2dec(x), codeWords);
numHybs = length(codeWords{1});

% -------------------------------------------------------------------------
% Compute integer codeword histogram
% -------------------------------------------------------------------------
[wordCounts,x] = hist([words.intCodeword], 1:2^numHybs);

% -------------------------------------------------------------------------
% Compute hamming sphere counts
% -------------------------------------------------------------------------
reportStruct.geneNames = geneNames;
reportStruct.hammingSphereCounts = zeros(parameters.maxHammingSphere + 1, ...
    length(geneNames));

reportStruct.hammingSphereCounts(1,:) = wordCounts(intCodewords);

for i=1:parameters.maxHammingSphere
    for j=1:length(geneNames)
        surroundingWords = GenerateSurroundingCodewords(codeWords{j}=='1', i);
        surroundingWordIntegers = cellfun(@(x) bin2dec(removeWs(num2str(x))), ...
            surroundingWords);
        
        reportStruct.hammingSphereCounts(i+1,j) = sum(wordCounts(surroundingWordIntegers));
    end
end

% -------------------------------------------------------------------------
% Identify blank words
% -------------------------------------------------------------------------
% Identify blank word geneIDs
isBlank = zeros(1, length(geneNames));
for i=1:length(parameters.blankWordIdentifiers)
    isBlank = cellfun(@(x) ~isempty(regexp(x, parameters.blankWordIdentifiers{i})), ...
        geneNames) | isBlank;
end
blankNames = geneNames(isBlank);
blankIDs = find(isBlank);
nonBlankIDs = setdiff(1:length(geneNames), blankIDs);

reportStruct.blankNames = blankNames;
reportStruct.blankIDs = blankIDs;
reportStruct.nonBlankIDs = nonBlankIDs;

% -------------------------------------------------------------------------
% Calculate histograms and CDF for 0/(1+0) ratio
% -------------------------------------------------------------------------
hammingSphereRatio = reportStruct.hammingSphereCounts(1,:)./ ...
    sum(reportStruct.hammingSphereCounts([1 2],:),1);

[nNonBlank, x] = hist(hammingSphereRatio(~isBlank), ...
    linspace(0, max(hammingSphereRatio), parameters.numHistogramBins));
[nBlank, x] = hist(hammingSphereRatio(isBlank), x);

reportStruct.hammingSphere01Ratio = hammingSphereRatio;
reportStruct.nNonBlank = nNonBlank;
reportStruct.nBlank = nBlank;
reportStruct.histBins = x;

reportStruct.NonBlankCDF = cumsum(nNonBlank)/sum(nNonBlank);
reportStruct.blankCDF = cumsum(nBlank)/sum(nBlank);

[~, sind] = sort(reportStruct.hammingSphere01Ratio, 'descend');
reportStruct.sortedInd = sind;

% -------------------------------------------------------------------------
% Display Report
% -------------------------------------------------------------------------
reportStruct.figHandles = [];
reportID = find(strcmp('hammingSphereReport', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    reportStruct.figHandles(end+1) = figure('Name','hammingSphereReport', ...
        'visible', parameters.reportsToGenerate{reportID, 2}, ...
        'Position', [275         434        1419         548]);
    
    subplot(2,2, [1 2]);
    plot(reportStruct.hammingSphere01Ratio(sind), 'b'); hold on;
    revSind(sind) = 1:length(sind);
    plot(revSind(blankIDs), reportStruct.hammingSphere01Ratio(blankIDs), 'rx');
    xlabel('Sorted Gene ID');
    ylabel('Fraction');
    xlim([0 length(geneNames)+1]);
    ylim([0 1.2*max(reportStruct.hammingSphere01Ratio)]); 
    
    subplot(2,2,3);
    bar(reportStruct.histBins, reportStruct.nNonBlank, 'FaceColor', 'b', 'EdgeColor', 'none'); hold on;
    bh = bar(reportStruct.histBins, reportStruct.nBlank, 'FaceColor', 'r', 'EdgeColor', 'none'); hold on;
    alpha(get(bh, 'Children'), .5);
    xlabel('Hamming Sphere Fraction');
    ylabel('Counts');
    xlim([reportStruct.histBins(1) - mean(diff(reportStruct.histBins)) ...
        reportStruct.histBins(end) + mean(diff(reportStruct.histBins))]);
    ylim([0 1.2*max([reportStruct.nNonBlank reportStruct.nBlank])]);
    
    subplot(2,2,4);
    plot(reportStruct.histBins, reportStruct.NonBlankCDF, 'b'); hold on;
    plot(reportStruct.histBins, reportStruct.blankCDF, 'r'); 
    ind = find(reportStruct.blankCDF == 1, 1, 'first');
    plot(reportStruct.histBins(ind)*ones(1,2), [0 1], 'k--'); 

    xlabel('Hamming Sphere Fraction');
    ylabel('CDF');
    
    numAbove = sum(reportStruct.hammingSphere01Ratio(reportStruct.nonBlankIDs) > ...
        max(reportStruct.hammingSphere01Ratio(reportStruct.blankIDs)));
    
    xlim([min(reportStruct.histBins) max(reportStruct.histBins)]);
    ylim([0 1]);
    
    title(['Number Above: ' num2str(numAbove)]);
    
    PresentationPlot();
    
    if parameters.saveAndClose
        SaveFigure(reportStruct.figHandles(end),'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats);

        close(reportStruct.figHandles(end));
        reportStruct.figHandles(end) = -1;
    end

end
% -------------------------------------------------------------------------
% Display Confidence Ratio Report
% -------------------------------------------------------------------------
reportID = find(strcmp('confidenceRatioReport', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    % Generate figure handle
    reportStruct.figHandles(end+1) = figure('Name','confidenceRatioReport', ...
    'visible', parameters.reportsToGenerate{reportID, 2}, ...
    'Position', [275         434        1419         548]);

    % Compute ratios for real and blank words
    allRatio = reportStruct.hammingSphereCounts(1,:)./ ...
        reportStruct.hammingSphereCounts(2,:);
    exact2CorrectedGenes = allRatio(reportStruct.nonBlankIDs);
    blank2CorrectedGenes = allRatio(reportStruct.blankIDs);

    % Sort
    geneRatio = sort(exact2CorrectedGenes,'descend');
    blankRatio = sort(blank2CorrectedGenes,'descend');

    % Plot
    bar(1:length(geneRatio), geneRatio/max(blankRatio), 1, 'EdgeColor', 'none', 'FaceColor', 'b'); hold on;
    bar([1:length(blankRatio)]+length(geneRatio), blankRatio/max(blankRatio), 1, 'EdgeColor', 'none', 'FaceColor', 'r'); hold on;

    maxBlank = max(blankRatio);
    plot([1 length(allRatio)], ones(1,2), 'k--');
    xlim([0 length(allRatio)+1]);
    xlabel('Gene ID');
    ylabel('Normalized Confidence Ratio');
    set(gca, 'XTick', [1 35 70 105 140]);
    PresentationPlot();
    
    if parameters.saveAndClose
        SaveFigure(reportStruct.figHandles(end),'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats);

        close(reportStruct.figHandles(end));
        reportStruct.figHandles(end) = -1;
    end

end


