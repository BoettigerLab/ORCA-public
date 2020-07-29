function [fig_handle] = ImbedNamesInFigure(fig_handle, infoStruct, varargin)
%--------------------------------------------------------------------------
% [fig_handle] = ImbedNamesInFigure(varargin)
% This function embeds information in a figure and overloads the text
% output function of the data cursor to allow additional information to be
% displayed with each point
%--------------------------------------------------------------------------
% Necessary Inputs
% fig_handle: A handle to the figure to modify
%
% infoStruct: An information structure. This structure must have the
% following properties:
%   1. The elements must be in the same order as the data points plotted
%   2. There must be a field with the name 'featureName' that contains the
%   name of each data point.  
%
%--------------------------------------------------------------------------
% Outputs
% cvStruct/structure array: A structure containing counts per bp for both
%   strands of a given reference sequence. Each element in the array
%   corresponds to a different reference sequence. 
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% -allSubPlots/boolean/true: Display information for all subplots
%
% -includeXY/boolean(false): Display the x and y position as well
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% Started: February 11, 2013
% Revised: December 3, 2013
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
flags = {'verbose', 'allSubPlots'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = true;
allSubPlots = true;
global includeXY; %Must be global to be seen by display function
includeXY = false;

%--------------------------------------------------------------------------
% Parse file information
%--------------------------------------------------------------------------
if nargin < 2
    error([' ''' mfilename ''' requires a figure handle and a infoStruct']);
end

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch parameterName
        case 'verbose'
            verbose = CheckParameter(parameterValue,'boolean','verbose');
        case 'allSubPlots'
            allSubPlots = CheckParameter(parameterValue,'boolean','allSubPlots');
        case 'includeXY'
            includeXY = CheckParameter(parameterValue, 'boolean', 'includeXY');
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end

%--------------------------------------------------------------------------
% Imbed infoStruct in figure
%--------------------------------------------------------------------------
set(fig_handle, 'UserData', infoStruct);
datacursormode off
dcm_obj = datacursormode(fig_handle);
set(dcm_obj, 'UpdateFcn', @DisplayFeatureName);

z_data = 1:length(infoStruct);
children = get(gca, 'Child');
if allSubPlots
    for i=1:length(children)
        if isprop(children(i), 'ZData')
            set(children(i), 'ZData', z_data); 
        end
    end
else
    set(children(1), 'ZData', z_data); % Only imbed on top data
end

datacursormode on

end


%--------------------------------------------------------------------------
% Display Functions
%--------------------------------------------------------------------------
function output_txt = DisplayFeatureName(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
    global includeXY;
    pos = get(event_obj,'Position');
    infoStruct = get(gcf, 'UserData');
    switch class(infoStruct)
        case 'struct'
            if isfield(infoStruct, 'name')
                output_txt = {[infoStruct(pos(3)).name]};
            elseif isfield(infoStruct, 'featureName')
                output_txt = {[infoStruct(pos(3)).featureName]};
            end
        case 'cell'
            output_txt = infoStruct(pos(3));
    end
    if includeXY
        output_txt = {['X: ' num2str(pos(1))], ...
            ['Y: ' num2str(pos(2))], ...
            [output_txt{1}]};
    end

end

%--------------------------------------------------------------------------
% Default Function
%--------------------------------------------------------------------------
function output_txt = testfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

display('I am here!');
end

