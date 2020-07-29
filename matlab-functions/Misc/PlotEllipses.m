function circHandles = PlotEllipses(x, y, r, varargin)
% ------------------------------------------------------------------------
% circHandles = PlotEllipses(x, y, r, varargin)
% This function plots ellipses of a defined radi r at positions x,y
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
if nargin < 3 || ~(size(x,1) == size(y,1) & size(x,1) == size(r,1))
    error('matlabFunctions:invalidArguments', 'x, y, and r are required and must have the same length.');
end

if size(r,2) ~= 2
    error('matlabFunctions:invalidArguments', 'r must be an Nx2 matrix');
end

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
theta = linspace(0,2*pi,parameters.numPointsPerCircle);
x0 = sin(theta);
y0 = cos(theta); 
for i=1:length(x)
    circHandles(i) = fill(r(i,1)*x0 + x(i), r(i,2)*y0 + y(i),'b', 'EdgeColor', 'none'); hold on;
    if ~isempty(parameters.commonPatchOptions)
        set(circHandles(i), parameters.commonPatchOptions{:});
    end
    if ~isempty(parameters.variablePatchOptions)
        set(circHandles(i), parameters.variablePatchOptions{i}{:});
    end
end 
    
    