function [figHandles, parameters] = GenerateTiledImage(imageData, varargin)
% ------------------------------------------------------------------------
% [figHandles, parameters] = GenerateTiledImage(imageData, varargin)
% This function generates a tiled image composed of all hybs for each cell
% in the imageData structure. 
%--------------------------------------------------------------------------
% Necessary Inputs
% imageData/A structure array with an element for each image used to create
%   the elements in words.  See CreateWordsStructure for information on
%   field names. 
%--------------------------------------------------------------------------
% Outputs
% hybReport/A structure containing information from the report
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% October 2, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'compositeHybImage', 'on'}; % {report name, 'off'/'on' do not/do display figure}

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate','cell', defaultReports};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'showNames', 'boolean', false};
defaults(end+1,:) = {'embedNames', 'boolean', true};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'numImageColumns', 'positive', 4};
defaults(end+1,:) = {'displayHybLabel', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~isstruct(imageData) || ...
        ~isempty(setdiff(fields(CreateImageDataStructure(0)), fields(imageData)) ))
    error('matlabFunctions:invalidArguments', 'Invalid words or imageData structures.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Disable warnings
% -------------------------------------------------------------------------
warning('off', 'images:initSize:adjustingMag');

% -------------------------------------------------------------------------
% Generate total report
% -------------------------------------------------------------------------
reportID = find(strcmp('compositeHybImage', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    cellIDs = unique([imageData.cellNum]);
    
    for i=1:length(cellIDs)
        % Create figure handle. 
        figHandles(i) = figure('Name', ['CompositeImage_Cell_' num2str(cellIDs(i))], ...
            'visible', parameters.reportsToGenerate{reportID,2});
        
        % Index image Data and sort by hyb number
        localImageData = imageData([imageData.cellNum] == cellIDs(i));
        hybNums = [localImageData.hybNum];
        [~, sind] = sort(hybNums, 'Ascend');
        localImageData = localImageData(sind);
    
        % Create composite figure
        imageDim = [localImageData(1).imageW localImageData(1).imageH];
        numColumns = parameters.numImageColumns;
        numRows = ceil(length(localImageData)/parameters.numImageColumns);
        
        newImage = zeros([numColumns numRows].*imageDim);
        
        labelPos = [];
        for j=1:length(localImageData)
            [rowInd, colInd] = ind2sub([numRows numColumns], j);
            indsX = (0:(imageDim(1)-1)) + 1 + (rowInd-1)*imageDim(1); 
            indsY = (0:(imageDim(2)-1)) + 1 + (colInd-1)*imageDim(2);
            dax{i} = ReadDax(localImageData(j).infFilePath, 'startFrame', 1, 'endFrame', 1, ...
                'verbose', false);
            newImage(indsY, indsX) = dax{i};
            labelPos(j,:) = [indsX(end) - imageDim(1)/10, indsY(1) + imageDim(2)/10];
        end
        
        imshow(newImage, []);
        
        if parameters.displayHybLabel
            for j=1:length(localImageData)
                text(labelPos(j,1), labelPos(j,2), num2str(j), 'Color', 'r'); hold on
            end
        end
    
        colorbar;
        PresentationPlot();
        
        if parameters.saveAndClose
                if parameters.useSubFolderForCellReport
                    SaveFigure(figHandles(i), 'overwrite', parameters.overwrite, ...
                        'formats', parameters.figFormats, ...
                        'subFolder', ['Cell_' num2str(cellIDs(i))]);
                else
                    SaveFigure(figHandles(i), 'overwrite', parameters.overwrite, ...
                        'formats', parameters.figFormats);
                end
            close(figHandles(i));
            figHandles(i) = -1;
        end   
    end
end

% -------------------------------------------------------------------------
% Reenable Warnings
% -------------------------------------------------------------------------
warning('on', 'images:initSize:adjustingMag');
