function [reportStruct, parameters] = GenerateBlankWordReport(words, geneNames, varargin)
% ------------------------------------------------------------------------
% [reportStruct, parameters] = GenerateBlankWordReport(words, varargin)
% This function creates a report about the relative abundance of blank
% words
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
% October 3, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'blankWordsReport', 'on'}; %

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate','cell', defaultReports};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'blankWordIdentifiers', 'cell', {'blank', 'notarget'}};
defaults(end+1,:) = {'numBins', 'positive', 10};
defaults(end+1,:) = {'subFolder', 'string', []};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~isstruct(words) || ...
        ~isempty(setdiff(fields(StripWords(CreateWordsStructure(0,0))), fields(words)) ))
    error('matlabFunctions:invalidArguments', 'Invalid words structure.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Blank Word Report
% -------------------------------------------------------------------------
words = words([words.isExactMatch] | [words.isCorrectedMatch]);

% -------------------------------------------------------------------------
% Create gene name map and determine gene IDs
% -------------------------------------------------------------------------
geneName2ID = containers.Map(geneNames, 1:length(geneNames));
geneIDs = cellfun(@(x) geneName2ID(x), {words.geneName}, 'ErrorHandler', @ErrorFunc);

% Identify blank word geneIDs
isBlank = zeros(1, length(geneNames));
for i=1:length(parameters.blankWordIdentifiers)
    isBlank = cellfun(@(x) ~isempty(regexp(x, parameters.blankWordIdentifiers{i})), ...
        geneNames) | isBlank;
end
blankNames = geneNames(isBlank);
blankIDs = cellfun(@(x) geneName2ID(x), blankNames);
nonBlankIDs = setdiff(1:length(geneNames), blankIDs);

[numCountsPerGene, x] = hist(geneIDs, 0:length(geneNames));
numCountsPerGene = numCountsPerGene(2:end); % Remove unlisted words

% Archive results
reportStruct.blankNames = blankNames;
reportStruct.blankIDs = blankIDs;
reportStruct.nonBlankNames = geneNames(~isBlank);
reportStruct.nonBlankIDs = nonBlankIDs;
reportStruct.numCountsPerGene = numCountsPerGene;
reportStruct.nonBlankCounts = numCountsPerGene(nonBlankIDs);
reportStruct.blankCounts = numCountsPerGene(blankIDs);

% Find number non blank above most abundant blank
[maxBlank] = max(numCountsPerGene(blankIDs));
numAbove = sum(numCountsPerGene(nonBlankIDs) > maxBlank);
reportStruct.maxBlank = maxBlank;
reportStruct.numAbove = numAbove;
reportStruct.fractAbove = numAbove/length(nonBlankIDs);

reportID = find(strcmp('blankWordsReport', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    reportStruct.figHandles = figure('Name','blankWordsReport', ...
        'visible', parameters.reportsToGenerate{reportID, 2}, ...
        'Position', [275         562        1419         420]);
    
    subplot(1,2,1);
    localCounts = numCountsPerGene;
    localCounts(localCounts==0) = 0.5;
    numberOfGenes = length(localCounts);
    numberOfBlanks = length(blankIDs);
    bar(1:(numberOfGenes-numberOfBlanks), localCounts(nonBlankIDs), 1, 'FaceColor', 'b', 'EdgeColor', 'none', 'BaseValue', 0.5); hold on;
    bar((numberOfGenes-numberOfBlanks + 1):numberOfGenes, localCounts(blankIDs), 1, 'FaceColor', 'r', 'EdgeColor', 'none','BaseValue', 0.5);
    xlabel('Gene ID');
    ylabel('Counts');
    set(gca, 'YScale', 'log');
    xlim([0 length(localCounts)+1]);
    ylim([0.5 2*max(localCounts)]);
    
    subplot(1,2,2);
    data1 = log10(numCountsPerGene(nonBlankIDs));
    data2 = log10(numCountsPerGene(blankIDs));
    
    goodInds = ~isinf(data1);
    data1 = data1(goodInds);
    
    [n1, x] = hist(data1, linspace(min(data1), max(data1), parameters.numBins));
    bh = stairs(x, n1, 'b'); hold on;
    
    goodInds = ~isinf(data2);
    data2 = data2(goodInds);

    [n2, x] = hist(data2, x);
    bh = stairs(x, n2, 'r');

    xlabel('Counts (log_{10})');
    ylabel('Counts');
    
    title(['Number above: ' num2str(reportStruct.numAbove)]);
    % 
    
    PresentationPlot();
    
    if parameters.saveAndClose
        SaveFigure(reportStruct.figHandles,'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats, ...
            'subFolder', parameters.subFolder);

        close(reportStruct.figHandles);
        reportStruct.figHandles = -1;
    end

end

end


function output = ErrorFunc(errStruct, geneName)
output = 0;
end
