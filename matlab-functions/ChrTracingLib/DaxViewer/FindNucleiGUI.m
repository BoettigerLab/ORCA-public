function cellBorderTable = FindNucleiGUI(imIn,varargin)
% a simple joiner function for FindCirclesApp and FindCellBordersApp

% parameters
defaults = cell(0,3);
% shared
defaults(end+1,:) = {'figHandle','handle',200};
defaults(end+1,:) = {'contrastMin','fraction',.3};
defaults(end+1,:) = {'contrastMax','fraction',.999};
defaults(end+1,:) = {'colormap','colormap','gray'};
% FindCircles
defaults(end+1,:) = {'sensitivity','fraction',.9};
defaults(end+1,:) = {'minRadius','integer',6};
defaults(end+1,:) = {'maxRadius','integer',15};
defaults(end+1,:) = {'downsample','integer',4};
% FindCellBorders
defaults(end+1,:) = {'nAngles','integer',21};
defaults(end+1,:) = {'segLength','integer',20};
defaults(end+1,:) = {'smoothFindEdge','integer',3};
defaults(end+1,:) = {'smoothFinalBorder','integer',3};
defaults(end+1,:) = {'imageFilter','array',fspecial('gaussian',10,20)};
defaults(end+1,:) = {'dilate','integer',1};
defaults(end+1,:) = {'xyNucWindow','positive',2.25}; % search window is the nuc radius x this muliplier
defaults(end+1,:) = {'convex','boolean',false};
defaults(end+1,:) = {'showNum','boolean',true};
defaults(end+1,:) = {'showOverlayFig','integer',3};
defaults(end+1,:) = {'showExtraFig','integer',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% function

fc = FindCirclesApp(imIn,'parameters',pars);
uiwait(fc.pars.stayOpen);
cntrs = evalin('base','cntrs'); % fc.data.cntrs*fc.pars.downsample;
radii = evalin('base','radii'); %  fc.data.radii*fc.pars.downsample;
% cntrs, radii, should be created in caller. doesn't seem to work. 
fb = FindCellBordersApp(imIn,cntrs,radii,'parameters',pars);
% cellBorders = fb.data.cellBorders;
uiwait(fb.pars.figHandle);
try
    cellBorders = evalin('base','cellBorders');
catch
    cellBorders = [];
end

% write cell borders as a table
nCells = size(cellBorders,1); % nCells = 0
cellNum = repmat(1:nCells, [pars.nAngles+1,1]);
cellNum = cellNum(:);
xy = cat(1,cellBorders{:});
x = xy(:,1); y=xy(:,2);
cellBorderTable = table(x,y,cellNum);

