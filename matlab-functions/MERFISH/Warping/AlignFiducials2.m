function [fiducialData, parameters] = AlignFiducials2(fiducialData,varargin)
% ------------------------------------------------------------------------
% [fiducialData] = AlignFiducials(fiducialData,varargin)
% This function analyzes a series of raw conventional images in the 
%   specified directory and creates a words structure which represents all
%   of the identified words in the data.
%--------------------------------------------------------------------------
% Necessary Inputs
% fiducialData/A structure array with elements equal to the number of
%   images to align. This structure must contain an mList field. 
%--------------------------------------------------------------------------
% Outputs
% fiducialData/The same structure array provided but with two new fields
%   added.
%   --transform. The transform that will bring each image into the
%     reference frame of the first image. 
%   --warpUncertainty. A matrix that describes the quality of the generated
%     transform. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% September 6, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};

defaults(end+1,:) = {'reportsToGenerate','cell',cell(0,1)};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~(strcmp(class(fiducialData), 'struct')) || ...
        (length(fiducialData) < 2) || ...
        ~ismember('mList', fields(fiducialData)) )
    
    error('matlabFunctions:AlignFiducials', 'Invalid fiducial data.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
if parameters.printedUpdates & parameters.verbose
    display('--------------------------------------------------------------');
    display('Analyzing fiducials');
end

% -------------------------------------------------------------------------
% Handling optional variables
% -------------------------------------------------------------------------
if ~isfield(parameters,'numHybs')
    parameters.numHybs = length(fiducialData); 
end

% -------------------------------------------------------------------------
% Optional plotting
% -------------------------------------------------------------------------
reportID = find(strcmp('fiducialReport1', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandle1 = figure('Name',['fiducialReport1_cell', num2str(fiducialData(1).cellNum)], ...
        'visible', parameters.reportsToGenerate{reportID, 2});
    clrmap = jet(length(fiducialData)); 
    r = 5; % radius for color wheels of aligned hybes
end

% -------------------------------------------------------------------------
% Build transforms
% -------------------------------------------------------------------------
for i=1:length(fiducialData)
    [fiducialData(i).tform, shiftedList, residuals] = MLists2Transform(fiducialData(1).mList, ...
        fiducialData(i).mList, 'ignoreFrames',true, 'transpose', true);
    fiducialData(i).warpErrors = mean(sqrt(residuals(:,1).^2 + residuals(:,2).^2));

    fiducialData(i).hasFiducialError = false;
    fiducialData(i).fiducialErrorMessage = [];
    if parameters.verbose
        disp(['Alignment error = ',num2str(fiducialData(i).warpErrors)]);
    end
    
    % -------------------------------------------------------------------------
    % Plot Fiducial Report 1 Data: Bead Positions for Each Hyb
    % -------------------------------------------------------------------------
    if exist('figHandle1')
        set(0, 'CurrentFigure', figHandle1);
        plot(shiftedList.xc,shiftedList.yc,'.','color',clrmap(i,:)); hold on;  
    end
    
end
  
% -------------------------------------------------------------------------
% Finalize report figures and save
% -------------------------------------------------------------------------    
if exist('figHandle1')
    xlabel('pixels'); ylabel('pixels'); 
    
    if parameters.saveAndClose
        if parameters.useSubFolderForCellReport
            SaveFigure(figHandle1, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats, ...
                'subFolder', ['Cell_' num2str(fiducialData(1).cellNum)]);
        else
            SaveFigure(figHandle1, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats);
        end

        close(figHandle1); 
    end
end

end
        
