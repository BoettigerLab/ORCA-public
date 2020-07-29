function [reportStruct, parameters] = GenerateMisidentificationReport(words, exactMap, probabilities, varargin)
% ------------------------------------------------------------------------
% [reportStruct, parameters] = GenerateMisidentificationReport(words, exactMap, probabilities, varargin)
% This function uses bit flip probabilities and words to estimate the
% fraction of counts for each species that may have come from
% misidentification.  
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
% February 4, 2015
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'misidentificationRates', 'on'}; %

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
defaults(end+1,:) = {'verbose', 'boolean', false};

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
% Compute misidentification rates
% -------------------------------------------------------------------------
numHybes = size(probabilities,1);
misIdentRates = zeros(length(codewords));
for i=1:length(codewords) % Loop over codewords for genes
    % Generate exact and corrected words for the measured word
    measuredCodewords = {};
    measuredCodewords{1} = codewords{i};
    if ~isempty(parameters.errCorrFunc)
        correctedWords = parameters.errCorrFunc(codewords{i});
        measuredCodewords(2:(length(correctedWords)+1)) = correctedWords;
    end     
    
    % Loop over all possible real words (except the actual one) 
    for k=setdiff(1:length(codewords), i)
        originalCodeword = codewords{k}; % The real word that is misidentified as the measured word
        
        % Initialize accumulator for error probability
        localProb = 0;

        % Loop over all words that are mapped to the given gene
        for j=1:length(measuredCodewords)
            measuredWord = measuredCodewords{j}; % The word that is actually measured
        
            % Find the bits that differ and the direction of the error
            doBitsDiff = find(measuredWord ~= originalCodeword);
            errorDirection = ~originalCodeword(doBitsDiff) + 1; % Create probability inidices: 1= 1->0; 2= 0->1;

            % Extract the specific error probabilities
            localBitFlipErrors = probabilities(...
                sub2ind(size(probabilities), doBitsDiff, errorDirection));
            
            % Compute probability and add to accumulator
            if ~isempty(localBitFlipErrors)
                localProb = localProb + prod(localBitFlipErrors);
            end
        end
        misIdentRates(i,k) = localProb;
    end
    
    if parameters.verbose
        display(['Completed analysis of ' num2str(i) ': ' geneNames{i}]);
    end
    
end
reportStruct.misIdentRates = misIdentRates;

% -------------------------------------------------------------------------
% Estimate number of counts due to misidentification 
% -------------------------------------------------------------------------
reportStruct.misIdentCounts = misIdentRates*geneCounts';

% -------------------------------------------------------------------------
% Estimate misidentification fractions
% -------------------------------------------------------------------------
reportStruct.misIdentFract = reportStruct.misIdentCounts./geneCounts';

% -------------------------------------------------------------------------
% Display report figure
% -------------------------------------------------------------------------
figCount = 1;
reportID = find(strcmp('misidentificationRates', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    reportStruct.figHandles(figCount) = figure('Name',['misidentificationRates'], ...
        'visible', parameters.reportsToGenerate{reportID, 2}, ...
        'Position', [680         558        1188         420]);
    
    % Plot misidentification rates
    subplot(1,4,[1 2]);
    imagesc(reportStruct.misIdentRates);
    colorbar;
    zlabel('Probability');
    
    % Plot misidentification counts and counts
    subplot(1,4,3);
    bh = bar(reportStruct.geneCounts, 1, 'b', 'EdgeColor', 'none'); hold on;
    bh = bar(reportStruct.misIdentCounts, 1, 'r', 'EdgeColor', 'none'); hold on;
    xlabel('Gene ID');
    ylabel('Counts');
    set(gca, 'YScale', 'log');
    xlim([0 length(geneNames)+1]);
    
    % Plot fractions
    subplot(1,4,4);
    bh = bar(reportStruct.misIdentFract, 1, 'c', 'EdgeColor', 'none'); hold on;
    xlabel('Gene ID');
    ylabel('Fraction');
    xlim([0 length(geneNames)+1]);

    PresentationPlot();

    if parameters.saveAndClose
        SaveFigure(reportStruct.figHandles(figCount),'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats);
        close(bitFlipReport.figHandles(figCount));
        bitFlipReport.figHandles(figCount) = -1;
    end

    figCount = figCount + 1;
end
