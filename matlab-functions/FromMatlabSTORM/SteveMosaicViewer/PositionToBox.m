function box_coords = PositionToBox(position,mosaicPars,varargin)
% 
% Inputs
% takes input postions, an Nx2 array of x,y stage positions and mosaicPars,
% a structure created by LoadTiles which contains the scaling and
% translation information necessary to align position on mosaic. 
% 
% Outputs
% returns box_coordinates, an Nx4 array to plot with rectangle function.
% data is organized as x_min, y_min, width, height.

defaults = cell(0,3);
defaults(end+1,:) = {'frameSize','positive',1024};
pars = ParseVariableArguments(varargin,defaults,mfilename);

numPts = size(position,1); 

box_cx = mosaicPars.mx*position(:,1)-mosaicPars.xmin-pars.frameSize;
box_cy = mosaicPars.my*position(:,2)-mosaicPars.ymin-pars.frameSize;
box_coords = [box_cx,box_cy,pars.frameSize*ones(numPts,1),pars.frameSize*ones(numPts,1)];

% if showbox
% hold on;
%   rectangle('Position',box_coords,'EdgeColor','w');
% end
