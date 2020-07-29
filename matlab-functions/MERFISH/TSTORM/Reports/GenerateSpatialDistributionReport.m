function figHandles = GenerateSpatialDistributionReport(words, imageData, geneNames, varargin)
% ------------------------------------------------------------------------
% [figHandles, parameters] = GenerateSpatialDistributionReport(words, imageData, geneNames, varargin)
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
defaultReports(end+1,:) = {'geneDistributionReport', 'on'}; %

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
% Generate Gene Distribution Report Cell By Cell
% -------------------------------------------------------------------------
figCount = 1;
reportID = find(strcmp('geneDistributionReport', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    % Find cells
    cellIDs = sort(unique([words.cellID]));
    if ~isempty(parameters.cellsToKeep)
        cellIDs = intersect(cellIDs, parameters.cellsToKeep);
    end
    for i=1:length(cellIDs)
        figHandles(figCount) = figure('Name',['geneDistributionReport_Cell_' num2str(cellIDs(i))], ...
            'visible', parameters.reportsToGenerate{reportID, 2}, ...
            'Position', [1 41 1920 969]);
        set(figHandles(figCount), 'Color', 'w');

        % Derive composite image
        compositeImage = squeeze(CombineHybImages(imageData, cellIDs(i), 'colormap', @gray, ...
            'parameters', parameters));
        
        [imageH, imageW, depth] = size(compositeImage);
        
        % Find localWords
        localWords = words([words.cellID] == cellIDs(i));
        
        % Replicate composite image in figure
        numColumns = parameters.numColumns;
        numRows = ceil(length(geneNames)/numColumns);
        
        if length(geneNames) < numColumns
            numColumns = length(geneNames);
        end
        
        % newImage
        bigImage = repmat(compositeImage, [numRows numColumns]);
        
        imshow(bigImage, [], 'InitialMagnification', 'fit'); hold on;
        
        axesHandles = [];
        for j=1:length(geneNames)
            [ind1, ind2] = ind2sub([numRows numColumns], j);

            xOffset = imageW*(ind2-1);
            yOffset = imageH*(ind1-1);
            % Find genes and plot
            inds = strcmp({localWords.geneName}, geneNames{j});
            if sum(inds) > 0
                circHandles = PlotCircles([localWords(inds).wordCentroidX] + xOffset, ...
                    [localWords(inds).wordCentroidY] + yOffset, 3*ones(sum(inds),1));
                set(circHandles, 'FaceColor', 'g', 'FaceAlpha', 0.5);
            end
            if parameters.displayName
                text(xOffset+imageW/20,yOffset+imageH/10, geneNames{j}, 'Color', 'g');
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
        end
        figCount = figCount + 1;
    end
end



