function [imBoxes,skipPts,xyz_out] = BoxSpots3D(imIn,xyz,varargin)
% inputs 
%   imIn - input 3D image, full res,
%   xyz - centroids of points to crop (try FindSpots3D to produce these)
% outputs
%   imBoxes - cell array 
%   skip - positions in original xyz that were skipped due to edge

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'r_xy','integer',6};%  box radius in xy
defaults(end+1,:) = {'r_z','integer',4}; % box radius in z 
defaults(end+1,:) = {'showPlot','boolean',false};
defaults(end+1,:) = {'removeNonFits','boolean',true};

% parse variable arguments
pars = ParseVariableArguments(varargin, defaults, mfilename);

% BoxSpots
nPts = size(xyz,1);
imBoxes = cell(nPts,1); 
% zeros(2*r_xy+1,2*r_xy+1,2*r_z+1,nPts,class(imIn));
% if you don't intend to slice a diminsion 4, its easier and faster to keep
% this as a cell array than asa 4D array. 
[yMax,xMax,zMax]=size(imIn);
skipPts = false(nPts,1);
for n=1:nPts
    yL = xyz(n,2)-pars.r_xy;
    yU = xyz(n,2)+pars.r_xy;
    xL = xyz(n,1)-pars.r_xy;
    xU = xyz(n,1)+pars.r_xy;
    zL = xyz(n,3)-pars.r_z;
    zU = xyz(n,3)+pars.r_z;
    if yL>0 && xL>0 && zL>0 && yU<yMax && xU<xMax && zU < zMax
        imBoxes{n} = imIn(yL:yU,xL:xU,zL:zU);
    else
        skipPts(n) = true; % helpful to get back the edge violations. 
    end
end
xyz_out = xyz;
if pars.removeNonFits
    imBoxes(skipPts) = [];
    xyz_out(skipPts,:) = [];
end