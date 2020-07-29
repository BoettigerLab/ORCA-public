function [imageData, parameters] = GenerateFocusLockQuality(imageData, varargin)
% ------------------------------------------------------------------------
% [imageData, parameters] = GenerateFocusLockQuality(imageData, varargin)
% This function analyzes the camera-based focus lock images for a set of
%   imageData and generates figures that summarize the observed lock images. 
%--------------------------------------------------------------------------
% Necessary Inputs
% imageData/A structure array with the following fields
%   --filePath: The full path to the file. The path to the focus lock
%   image will be built from this path
%
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% September 15, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for generic display
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};

% Parameters for building focus lock file names
defaults(end+1,:) = {'focusLockSuffix', 'string', '_lock_cam.png'};
defaults(end+1,:) = {'focusLockNameGenFun', 'function', ...
    @(x)x.filePath(1:(regexp(x.filePath,x.binType))-5)};

% Parameters for determinity focus locking quality
defaults(end+1,:) = {'focusLockQualityFunc', 'function', @(x,y)corr2(x,y)};

% Parameters for controlling report generation
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'focusLockReportByCell', 'off'};
defaultReports(end+1,:) = {'focusLockReportSummary', 'on'};
defaults(end+1,:) = {'reportsToGenerate', 'cell', defaultReports};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~isstruct(imageData) || ~ismember('filePath', fields(imageData))
    error('matlabFunctions:invalidArguments', 'An imageData structure array is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Determine number of cells to load
% -------------------------------------------------------------------------
cellIDs = unique([imageData.cellNum]);

% -------------------------------------------------------------------------
% Loop over cells
% -------------------------------------------------------------------------
for i=cellIDs
    % ---------------------------------------------------------------------
    % Generate printed updates
    % ---------------------------------------------------------------------
    if parameters.printedUpdates
        display('-------------------------------------------------------------');
        display(['Analyzing focus lock quality for cell ' num2str(i)]);
    end
    
    % ---------------------------------------------------------------------
    % Identfy imageData corresponding to this cell
    % ---------------------------------------------------------------------
    localInds = find([imageData.cellNum] == i);
    
    % ---------------------------------------------------------------------
    % Sort images by hyb number
    % ---------------------------------------------------------------------
    [~, sind] = sort([imageData(localInds).hybNum]);
    localInds = localInds(sind);
    
    % ---------------------------------------------------------------------
    % Load camera images
    % ---------------------------------------------------------------------
    images = {};
    for j=1:length(localInds)
        images{j} = imread([parameters.focusLockNameGenFun(imageData(localInds(j))) ...
            parameters.focusLockSuffix]);
    end
    
    % ---------------------------------------------------------------------
    % Calculate Similarity
    % ---------------------------------------------------------------------
    focusLockQuality = zeros(1, length(localInds));
    for j=1:length(localInds)
        focusLockQuality(j) = parameters.focusLockQualityFunc(images{1}, images{j});
        imageData(localInds(j)).focusLockQuality = focusLockQuality(j);
    end
    
    % ---------------------------------------------------------------------
    % Create figure if desired
    % ---------------------------------------------------------------------
    reportID = find(strcmp('focusLockReportByCell', parameters.reportsToGenerate(:,1)));
    if ~isempty(reportID)
        figHandle = figure('Name',['focusLockReport_cell', num2str(i)], ...
            'visible', parameters.reportsToGenerate{reportID, 2}, ...
            'Position', [10   10   850   650]);
        for j=1:length(localInds)
            subplot(4, ceil(length(localInds)/4), j);
            imagesc(images{j});
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            colormap gray;
            
            title(num2str(focusLockQuality(j), 3));
        end
        PresentationPlot();
        if parameters.useSubFolderForCellReport
            SaveFigure(figHandle,'overwrite',parameters.overwrite,...
                'formats',parameters.figFormats, ...
                'subFolder', ['Cell_' num2str(i)]);
        else
            SaveFigure(figHandle,'overwrite',parameters.overwrite,...
                'formats',parameters.figFormats);
        end
        close(figHandle); 
    end
end

% -------------------------------------------------------------------------
% Generate final report
% -------------------------------------------------------------------------
reportID = find(strcmp('focusLockReportSummary', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandle = figure('Name', ['focusLockQualitySummary'], ...
        'visible', parameters.reportsToGenerate{reportID, 2});
    plot3([imageData.cellNum], [imageData.hybNum], [imageData.focusLockQuality], '.');
    xlabel('Cell Number');
    ylabel('Hyb Number');
    zlabel('Focus Lock Quality');
    PresentationPlot();
    SaveFigure(figHandle,'overwrite',true,'formats',{'fig','png'});
    close(figHandle); 
end
