function figHandles = GenerateCombinedRNASpatialDistribution(words, imageData, geneGroups, varargin)
% ------------------------------------------------------------------------
% [figHandles, parameters] = GenerateSpatialDistributionReport(words, imageData, geneGroups, varargin)
% This function creates figures that display the spatial distribution of
% the provided words
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A structure array with an element for each word. See
%   CreateWordsStructure for information about fields. A stripped word
%   structure is also fine.  All that is required are the following fields: 
%   'geneName', 'wordCentroidX', 'wordCentroidY'
% imageData/ a structure array of imageData. CreateImageDataStructure for
%    information on fields. 
% geneNames/A cell array containing the desired geneNames to display in the
%    order in which they will be displayed
%--------------------------------------------------------------------------
% Outputs
% figHandles/Handles for the generated figures. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% October 25, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'RNAspatialReport', 'on'}; %

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
defaults(end+1,:) = {'colormap', 'function', @jet};
defaults(end+1,:) = {'numColumns', 'nonnegative', 10};
defaults(end+1,:) = {'cellsToKeep', 'array', []};
defaults(end+1,:) = {'displayName', 'boolean', true};
defaults(end+1,:) = {'exactMatchOnly', 'boolean', false};
defaults(end+1,:) = {'edges', 'cell', {linspace(0,256,25), linspace(0,256,25)}};
defaults(end+1,:) = {'contourColormap', 'array', jet(20)};
defaults(end+1,:) = {'displayContour', 'boolean', false};
defaults(end+1,:) = {'displayGroupsTogether', 'boolean', false};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);        

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 2 || ...
        ~isstruct(imageData) || ...
        ~isempty(setdiff(fields(CreateImageDataStructure(0)), fields(imageData)) ))
    error('matlabFunctions:invalidArguments', 'Invalid imageData structures.');
end

% -------------------------------------------------------------------------
% Determine gene group properties
% -------------------------------------------------------------------------
numGroups = size(geneGroups,1);

% -------------------------------------------------------------------------
% Generate Gene Distribution Report Cell By Cell
% -------------------------------------------------------------------------
figCount = 1;
reportID = find(strcmp('RNAspatialReport', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    % Find cells
    cellIDs = sort(unique([words.cellID]));
    if ~isempty(parameters.cellsToKeep)
        cellIDs = intersect(cellIDs, parameters.cellsToKeep);
    end
    for i=1:length(cellIDs)
        figHandles(figCount) = figure('Name',['combinedRNADistributionReport_Cell_' num2str(cellIDs(i))], ...
            'visible', parameters.reportsToGenerate{reportID, 2}, ...
            'Position', [1 41 1920 969]);
        set(figHandles(figCount), 'Color', 'w');

        % Derive composite image
        compositeImage = squeeze(CombineHybImages(imageData, cellIDs(i), ...
            'parameters', parameters));
        [imageH, imageW, depth] = size(compositeImage);

        % Find localWords
        localWords = words([words.cellID] == cellIDs(i));
        
        if parameters.displayGroupsTogether
            numColumns = 1;
            numRows = 1;
        else
            % Replicate composite image in figure
            numColumns = parameters.numColumns;
            numRows = ceil(numGroups/numColumns);
        end
        
        if numRows == 1
            if numGroups < numColumns
                numColumns = numGroups;
            end
        end
        
        if parameters.displayContour
            numDisplayRows = 2*numRows;
        else
            numDisplayRows = numRows;
        end
        
        % newImage
        bigImage = repmat(compositeImage, [numRows numColumns]);
        
        ax1 = subplot(numDisplayRows,numColumns,[1:numRows*numColumns]);
        imshow(bigImage, []); hold on;
        
        axesHandles = [];
        for j=1:numGroups
            % Identify genes
            geneNames = geneGroups{j,1};
            [ind1, ind2] = ind2sub([numRows numColumns], j);

            if parameters.displayGroupsTogether
                xOffset = 0;
                yOffset = 0;
            else
                xOffset = imageW*(ind2-1);
                yOffset = imageH*(ind1-1);
            end
            
            % Find genes and plot
            inds = ismember({localWords.geneName}, geneNames);
            if parameters.exactMatchOnly
                inds = inds & [localWords.isExactMatch];
            end
            if sum(inds) == 0
                continue;
            end
            
            % Compile positions
            pos = [[localWords(inds).wordCentroidX]; [localWords(inds).wordCentroidY];];
            
            % Plot
            axes(ax1);
            PlotCircles(pos(1,:) + xOffset, ...
                 pos(2,:) + yOffset, 2*ones(1,length(pos(1,:))), 'commonPatchOptions', geneGroups{j,2});
            
            % Display name if needed
            if ~parameters.displayGroupsTogether & parameters.displayName & ~isempty(geneGroups{j,3})
                text(xOffset+imageW/20,yOffset+imageH/10, geneGroups{j,3}, 'Color', 'g');
            end
            
            if parameters.displayContour
                % Compute histogram
                N = hist3(pos', 'Edges', parameters.edges);
                subplot(numDisplayRows,numColumns,j+numRows*numColumns);
                contour(N'/sum(N(:))); hold on;
                axis square
                set(gca, 'YDir', 'reverse');
            end
        end
        if parameters.saveAndClose
            if parameters.useSubFolderForCellReport
                SaveFigure(figHandles(figCount),'overwrite',parameters.overwrite,...
                    'formats',parameters.figFormats, ...
                    'subFolder', ['Cell_' num2str(cellIDs(i))]);
            else
                SaveFigure(figHandles(figCount),'overwrite',parameters.overwrite,...
                    'formats',parameters.figFormats);
            end
            close(figHandles(figCount));
            figHandles(figCount) = -1;
        end
        figCount = figCount + 1;
    end
end


