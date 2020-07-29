function [reportStruct, parameters] = GenerateDetectionEfficiencyReport(words, exactMap, probabilities, varargin)
% ------------------------------------------------------------------------
% [reportStruct, parameters] = GenerateDetectionEfficiencyReport(words, exactMap, probabilities, varargin)
% This function uses bit flip probabilities and words and exactMap to
% estimate the detection efficiency for every word in the codebook as a
% function of the measured abundance. 
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A structure array with an element for each word. See
%   CreateWordsStructure for information about fields. 
% exactMap/A container map: see CodebookToMap for details.
% probabilities: A 2xN vector of the per hyb (N) probability of making a
% 1->0 or 0->1 transition. 1->0 corresponds to probabilities(1,:). 
%--------------------------------------------------------------------------
% Outputs
% reportStruct/A structure containing information on the report. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% December 11, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'detectionEfficiency', 'on'}; %

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate','cell', defaultReports};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'probToUse', {'exact', 'firstOrder'}, 'exact'};
defaults(end+1,:) = {'errCorrFunc', 'function', []};
defaults(end+1,:) = {'numHistBins', 'positive', 10};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Defined required words fields
% -------------------------------------------------------------------------
requiredWordsFields = {'intCodeword'};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~isstruct(words) || ...
        ~isempty(setdiff(requiredWordsFields, fields(words))))
    error('matlabFunctions:invalidArguments', 'Words is missing some required fields.');
end

% -------------------------------------------------------------------------
% Normalize codeword type
% -------------------------------------------------------------------------
geneNames = values(exactMap);
codewords = keys(exactMap);

if ischar(codewords{1})
    codewords = cellfun(@(x)(x=='1'), codewords, 'UniformOutput', false);
else
    codewords = cellfun(@(x)de2bi(x,parameters.numHybs), codewords, 'UniformOutput', false);
end

% -------------------------------------------------------------------------
% Compute abundance for all genes in codebook
% -------------------------------------------------------------------------
% Create geneID map
geneIDMap = containers.Map(geneNames, 1:length(geneNames));

% Cut words on those that have been mapped to a gene
words = words([words.isExactMatch] | [words.isCorrectedMatch]);

% Find GeneID
geneIDs = cellfun(@(x) geneIDMap(x), {words.geneName});

% Calculate abundance
geneCounts = hist(geneIDs, 1:length(geneNames));

reportStruct.geneNames = geneNames;
reportStruct.geneCounts = geneCounts;

% -------------------------------------------------------------------------
% Compute detection efficiency
% -------------------------------------------------------------------------
numHybes = size(probabilities,1);
noErrorProb = zeros(1,length(codewords));
oneErrorProb = zeros(1, length(codewords));
for i=1:length(codewords)    
    inds = sub2ind(size(probabilities), 1:numHybes, 2-codewords{i}); % 1: 1->0 error; 2: 0->1 error
    localProbability = probabilities(inds);
    
    noErrorProb(i) = prod(1-localProbability);
    for j=1:numHybes
       oneErrorProb(i) = oneErrorProb(i) + ...
           localProbability(j)*prod(1-localProbability(setdiff(1:numHybes, j)));
    end
end
reportStruct.noErrorProb = noErrorProb;
reportStruct.oneErrorProb = oneErrorProb;

% -------------------------------------------------------------------------
% Compute detection efficiency
% -------------------------------------------------------------------------
reportStruct.exactDetectionEff = noErrorProb;
reportStruct.correctedDetectionEff = noErrorProb + oneErrorProb;

% -------------------------------------------------------------------------
% Compute histograms
% -------------------------------------------------------------------------
[n, x] = hist(reportStruct.exactDetectionEff, parameters.numHistBins);
reportStruct.exactHist = [x; n];

[n, x] = hist(reportStruct.correctedDetectionEff, parameters.numHistBins);
reportStruct.correctedHist = [x; n];

% -------------------------------------------------------------------------
% Display Probabilites
% -------------------------------------------------------------------------
figCount = 1;
reportID = find(strcmp('detectionEfficiency', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    reportStruct.figHandles(figCount) = figure('Name',['detectionEfficiency'], ...
        'visible', parameters.reportsToGenerate{reportID, 2}, ...
        'Position', [680         558        1188         420]);
    
    subplot(1,2,1);
    semilogx(reportStruct.geneCounts, reportStruct.exactDetectionEff, 'b.'); hold on;
    semilogx(reportStruct.geneCounts, reportStruct.correctedDetectionEff, 'r.');
    xlabel('Counts');
    ylabel('Detection Efficiency');
    xlim([min(reportStruct.geneCounts)*0.8 max(reportStruct.geneCounts)*1.2]);
    
    subplot(1,2,2);
    bh = bar(reportStruct.exactHist(1,:), reportStruct.exactHist(2,:), 'b'); hold on;
    bh = bar(reportStruct.correctedHist(1,:), reportStruct.correctedHist(2,:), 'r'); hold on;

    xlabel('Efficiency');
    ylabel('Counts');
    
    PresentationPlot();

    if parameters.saveAndClose
        SaveFigure(reportStruct.figHandles(figCount),'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats);
        close(bitFlipReport.figHandles(figCount));
        bitFlipReport.figHandles(figCount) = -1;
    end

    figCount = figCount + 1;
    
end
