function [figHandles, parameters] = GenerateHammingSpaceReport(words, exactMap, varargin)
% ------------------------------------------------------------------------
% [figHandles, parameters] = GenerateHammingSpaceReport(words, exactMap, varargin)
% This function creates a variety of 1D projections of hamming space
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
% October 12, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'hamming1DReportAllGenes', 'on'}; %
%defaultReports(end+1,:) = {'hamming1DReportGeneByGene', 'on'}; %

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

defaults(end+1,:) = {'hammingDist', 'positive', 2};
defaults(end+1,:) = {'hammingDisplayMethod', {'nested', 'concentric'}, 'nested'};
defaults(end+1,:) = {'genesToPlot', 'cell', {}};
defaults(end+1,:) = {'numHybs', 'positive', 16};
defaults(end+1,:) = {'gene2geneSpacing', 'positive', 10};
defaults(end+1,:) = {'includeHamming2', 'boolean', true};
defaults(end+1,:) = {'includeSpokes', 'boolean', true};

figHandles = [];

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
% if (nargin < 1 || ...
%         ~isstruct(words) || ...
%         ~isempty(setdiff(fields(CreateWordsStructure(0,0)), fields(words)) ))
%     error('matlabFunctions:invalidArguments', 'Invalid words structure.');
% end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Generate Histogram of Word Counts
% -------------------------------------------------------------------------
n = hist([words.intCodeword], 1:2^parameters.numHybs);

% -------------------------------------------------------------------------
% Find codewords for words
% -------------------------------------------------------------------------
exactCodewords = keys(exactMap);
if ischar(exactCodewords{1})
    exactCodewords = cellfun(@(x)bi2de(fliplr(x=='1')), exactCodewords);
else
    exactCodewords = [exactCodewords{:}];
end
geneNames = values(exactMap);

if isempty(parameters.genesToPlot)
    parameters.genesToPlot = geneNames;
end

% -------------------------------------------------------------------------
% Default display parameters
% -------------------------------------------------------------------------
figCount = 1;

% -------------------------------------------------------------------------
% All Genes
% -------------------------------------------------------------------------
reportID = find(strcmp('hamming1DReportAllGenes', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandles(figCount) = figure('Name',['hammingReportAllGenes'], ...
        'visible', parameters.reportsToGenerate{reportID, 2});

    for i=1:length(exactCodewords)
        % Hamming 0
        geneCodeWord = exactCodewords(i);
        
        [xPos0, yPos0] = ind2sub([ceil(sqrt(length(geneNames))) ceil(sqrt(length(geneNames)))], i);
        xPos0 = parameters.gene2geneSpacing*xPos0; 
        yPos0 = parameters.gene2geneSpacing*yPos0;
        
        % Plot (see function below)
        PlotHamming(n, geneCodeWord, parameters, xPos0, yPos0)
    end
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    xlim([-2.5 2.5]);
    ylim([-2.5 2.5]);
    grid off;
    box off;
    PresentationPlot();
    
    if parameters.saveAndClose
        if parameters.useSubFolderForCellReport
            SaveFigure(figHandles(figCount),'overwrite',parameters.overwrite,...
                'formats',parameters.figFormats);
        else
            SaveFigure(figHandles(figCount),'overwrite',parameters.overwrite,...
                'formats',parameters.figFormats);
        end

        close(figHandles(figCount));
        figHandles(figCount) = -1;
    end
    figCount = figCount+1;

end

% -------------------------------------------------------------------------
% Gene by Gene Report
% -------------------------------------------------------------------------
reportID = find(strcmp('hamming1DReportGeneByGene', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    for g = 1:length(parameters.genesToPlot)
        figHandles(figCount) = figure('Name',['hammingReport_' parameters.genesToPlot{g}], ...
            'visible', parameters.reportsToGenerate{reportID, 2});

        geneInd = find(strcmp(geneNames, parameters.genesToPlot{g}));
        geneCodeWord = exactCodewords(geneInd);
        
        % Plot (see function below)
        PlotHamming(n, geneCodeWord, parameters, 0, 0);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', {});
        set(gca, 'YTickLabel', {});
        xlim([-2.5 2.5]);
        ylim([-2.5 2.5]);
        grid off;
        box off;
        PresentationPlot();
        
        if parameters.saveAndClose
            if parameters.useSubFolderForCellReport
                SaveFigure(figHandles(figCount),'overwrite',parameters.overwrite,...
                    'formats',parameters.figFormats, 'subFolder', 'Genes\');
            else
                SaveFigure(figHandles(figCount),'overwrite',parameters.overwrite,...
                    'formats',parameters.figFormats);
            end
            close(figHandles(figCount));
            figHandles(figCount) = -1;
        end
        figCount = figCount+1;
    end 
end

end

function PlotHamming(n, geneCodeWord, parameters, xPos0, yPos0)
    r1 = 1;
    r2 = 2;
    
    % Hamming 0    
    stem3(xPos0,yPos0, n(geneCodeWord), 'b.'); hold on;

    % Hamming 1
    surroundingWords1 = cellfun(@bi2de, GenerateSurroundingCodewords(de2bi(geneCodeWord, parameters.numHybs)==1,1));
    theta = linspace(0, 2*pi, length(surroundingWords1)+1);
    theta = theta(1:(end-1));
    xPos1 = r1*sin(theta) + xPos0;
    yPos1 = r1*cos(theta) + yPos0;
    stem3(xPos1, ...
        yPos1, ...
        n(surroundingWords1), 'r.');
    plot(r1*sin(linspace(0, 2*pi, 100)), r1*cos(linspace(0, 2*pi, 100)), ...
        'k');

    % Hamming 2
    if parameters.includeHamming2
        data = [];
        pad = -ones(1,5);
        surroundingWords2 = [];
        for j=1:length(surroundingWords1)
            tempWords = cellfun(@bi2de, GenerateSurroundingCodewords(de2bi(surroundingWords1(j), parameters.numHybs)==1,1));
            ind = find(tempWords == geneCodeWord);
            surroundingWords2(j,:) = tempWords(setdiff(1:parameters.numHybs, ind));
            localData = n(tempWords);
            localData(ind) = -1;
            data = [data pad localData]; % Zero adds space
        end
        plot(r2*sin(linspace(0, 2*pi, 100)), r2*cos(linspace(0, 2*pi, 100)), 'k');

        theta = linspace(0, 2*pi, length(data)+1) - 2*pi/32;
        theta = theta(1:(end-1));
        goodInd = find(data~=-1);
        xPos2 = r2*sin(theta(goodInd)) + xPos0;
        yPos2 = r2*cos(theta(goodInd)) + yPos0;
        stem3(xPos2, ...
            yPos2, ...
            data(goodInd), 'g.');
    end

    % Include Spokes
    if parameters.includeSpokes
        for j=1:length(surroundingWords1)
            if surroundingWords1(j) < geneCodeWord
                plot([xPos0 xPos1(j)], [yPos0 yPos1(j)], 'k');
            end
            if parameters.includeHamming2
                for k=1:length(surroundingWords2(j,:))
                    if surroundingWords2(j,k) < surroundingWords1(j)
                        plot([xPos1(j) xPos2(k+(parameters.numHybs-1)*(j-1))], ...
                            [yPos1(j) yPos2(k+(parameters.numHybs-1)*(j-1))], 'k');
                    end
                end
            end
        end
    end

end 

