function vol = PolymerBinnedVolume(xyz,varargin)


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'binSize', 'positive', 40};
defaults(end+1,:) = {'zbinSize', 'positive', 80};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a nx3 data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);

% parameters
binSize = parameters.binSize; % 40;
zbinSize = parameters.zbinSize; % 80;


% main function
ymin = min(xyz(:,2));
ymax = max(xyz(:,2));
xmin = min(xyz(:,1));
xmax = max(xyz(:,1));
zmin = min(xyz(:,3));
zmax = max(xyz(:,3));
xBins = xmin:binSize:xmax;
yBins = ymin:binSize:ymax;
zBins = zmin:zbinSize:zmax;
mapIdx = hist4(xyz(:,2),xyz(:,1),xyz(:,3),'bins',cellfun(@length, {yBins,xBins,zBins}));  

bw = mapIdx > 0; 
vol = sum(bw(:))*(binSize*binSize*zbinSize);
