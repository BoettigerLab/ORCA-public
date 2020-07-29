function [I,vlist,imaxes] = Polymer2STORMimg(B,varargin)
%--------------------------------------------------------------------------
% Generate a STORM image from a set of polymer locations
% 
%--------------------------------------------------------------------------
% Optional Input
% {'zm', 'positive', 1};
% {'dotsize', 'positive', 16};
% {'N', 'positive', 1};
% {'buffer', 'nonnegative', .5};
% {'scalebar', 'nonnegative', 0};
% {'showPlots', 'boolean', true};

%% Default Parameters
% folder = 'C:\Users\Alistair\Documents\Research\Projects\OligoSecondaries\Data\2014-03-18_PolymerSamplingSims\';
% load([folder,'polymer_1.mat']);

defaults = cell(0,3);
defaults(end+1,:) = {'zm', 'positive', 1};
defaults(end+1,:) = {'dotsize', 'positive', 16};
defaults(end+1,:) = {'N', 'positive', 1};
defaults(end+1,:) = {'buffer', 'nonnegative', .5};
defaults(end+1,:) = {'scalebar', 'nonnegative', 0};
defaults(end+1,:) = {'showPlots', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'A Nx3 polymer is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


%% Main Function
B = B*parameters.zm; 
numNodes = size(B,1); 
vlist.xc = B(:,1);
vlist.yc = B(:,2);
vlist.zc = B(:,3);
vlist.a = 5*ones(numNodes,1);   
xSize = max(vlist.xc) - min(vlist.xc);
ySize = max(vlist.yc) - min(vlist.yc);
borderBuffer = max(xSize,ySize)*parameters.buffer;
imaxes.zm = 10; imaxes.scale = 2; 
imaxes.xmin = min(vlist.xc)-borderBuffer; imaxes.xmax =  max(vlist.xc)+borderBuffer;
imaxes.ymin = min(vlist.yc)-borderBuffer; imaxes.ymax = max(vlist.yc)+borderBuffer; 
im = list2img(vlist,imaxes,'dotsize',parameters.dotsize,'N',parameters.N,'scalebar',parameters.scalebar);
I = STORMcell2img(im,'colormap',hot(256));


if parameters.showPlots && nargout < 3
    imagesc(I);
end
