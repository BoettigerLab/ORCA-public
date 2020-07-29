function patchesHandles = PlotPatches(x, y, r, varargin)
% ------------------------------------------------------------------------
% patchesHandles = PlotPatches(x, y, r, varargin)
% This function plots circles of a defined radius r at positions x,y
%--------------------------------------------------------------------------
% Necessary Inputs
% x/array of x positions
% y/array of y positions
% r/array of radius (defaults to 1)
%--------------------------------------------------------------------------
% Outputs
% figHandles/Handles for the generated figures. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% November 22, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'commonPatchOptions', 'cell', {}}; % A list of patch options to apply to all circles
defaults(end+1,:) = {'variablePatchOptions', 'cell', {}}; % A list of patch options to apply to apply to each circle
defaults(end+1,:) = {'numPointsPerCircle', 'nonnegative', 50};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);        

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 3 || ~(length(x) == length(y) & length(x) == length(r))
    error('matlabFunctions:invalidArguments', 'x, y, and r are required and must have the same length.');
end

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
theta = linspace(0,2*pi,parameters.numPointsPerCircle);
x0 = sin(theta);
y0 = cos(theta); 
for i=1:length(x)
    patchesHandles(i) = fill(r(i)*x0 + x(i), r(i)*y0 + y(i),'b', 'EdgeColor', 'none'); hold on;
    if ~isempty(parameters.commonPatchOptions)
        set(patchesHandles(i), parameters.commonPatchOptions{:});
    end
    if ~isempty(parameters.variablePatchOptions)
        set(patchesHandles(i), parameters.variablePatchOptions{i}{:});
    end
end 
    
    