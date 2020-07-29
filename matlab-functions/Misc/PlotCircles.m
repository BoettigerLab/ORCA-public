function circHandles = PlotCircles(x, y, r, varargin)
% ------------------------------------------------------------------------
% circHandles = PlotCircles(x, y, r, varargin)
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
defaults(end+1,:) = {'z', 'nonnegative', []}; % A list of z positions
defaults(end+1,:) = {'colormap', 'colormap', lines(1)}; %
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);        
cirHandles = [];

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if length(r) == 1
    r = repmat(r,length(x),1);
end
    
if nargin < 3 || ~(length(x) == length(y) & length(x) == length(r))
    error('matlabFunctions:invalidArguments', 'x, y, and r are required and must have the same length.');
end

% -------------------------------------------------------------------------
% Check z if provided
% -------------------------------------------------------------------------
if ~isempty(parameters.z)
    if length(parameters.z) ~= length(x)
        error('matlabFunctions:invalidArguments', 'Provided z must be the same length as x.');
    end
end

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
colormap(parameters.colormap); 
theta = linspace(0,2*pi,parameters.numPointsPerCircle);
x0 = sin(theta);
y0 = cos(theta); 
for i=1:length(x)
    if isempty(parameters.z)
        circHandles(i) = fill(r(i)*x0 + x(i), r(i)*y0 + y(i),1, 'EdgeColor', 'none'); hold on;
    else
        circHandles(i) = fill3(r(i)*x0 + x(i), r(i)*y0 + y(i), parameters.z(i)*ones(1, length(theta)),1, 'EdgeColor', 'none'); hold on;
    end
    if ~isempty(parameters.commonPatchOptions)
        set(circHandles(i), parameters.commonPatchOptions{:});
    end
    if ~isempty(parameters.variablePatchOptions)
        set(circHandles(i), parameters.variablePatchOptions{i}{:});
    end
end 
    
    