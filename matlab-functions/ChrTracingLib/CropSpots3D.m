function [croppedSpots,spotXYZ] = CropSpots3D(im3D,spotXY,varargin)
%% inputs
% im3D   - a 3D image  (uncropped).
% spotXY - Nx2 array of all spots found in the image in 3D
%% outputs
% croppedSpots - Nx1 cell array of XxYxZ pixels arrays
%              spots near the edge will result in smaller (unpadded) arrays
%              pixel array is set by the optional input parameters 
%               spot is centered on brigthest pixel
%% optional output
% spotXYZ - contains xyz brightest pixel in coordiantes of irignal image

%% Main function

defaults = cell(0,3);
defaults(end+1,:) = {'xySearchWindow','integer',10}; % expand by this sq-radius the area cropped
defaults(end+1,:) = {'cropCenter',{'brightest','given'},'brightest'}; % expand by this sq-radius the area cropped
defaults(end+1,:) = {'xyCropWindow','integer',6}; % expand by this sq-radius the area cropped
defaults(end+1,:) = {'zCropWindow','integer',10}; % only take +/- this many pixels
defaults(end+1,:) = {'troubleshoot','boolean',false}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

% some sizes
nSpots = size(spotXY,1);
[ys,xs,zs] = size(im3D);
% crop the 3D image to record the z-position and do fine fit
croppedSpots = cell(nSpots,1);
x = spotXY(:,1); % a bit more compact notation
y = spotXY(:,2); 
spotXYZ = zeros(nSpots,3);
for s=1:nSpots     
    x1 = max(1 ,floor(x(s)-pars.xySearchWindow+1));
    x2 = min(xs,floor(x(s)+pars.xySearchWindow));
    y1 = max(1 ,floor(y(s)-pars.xySearchWindow+1));
    y2 = min(ys,floor(y(s)+pars.xySearchWindow+1));
    cropIm = im3D(y1:y2,x1:x2,:);     
     [~,id] = sort(cropIm(:),'descend') ; % data_3d = im3D;
     [yi,xi,zi] = ind2sub(size(cropIm),id(1));
    
    % just for troubleshooting
    if pars.troubleshoot
       figure(5); clf; 
       ProjectIm3D(cropIm,'showPlots',true);
       colormap(gray);
       pause(1);
    end

        % get peak in z
    z1 = max(1,floor(zi-pars.zCropWindow+1));
    z2 = min(zs,floor(zi+pars.zCropWindow));
    if strcmp(pars.cropCenter,'brightest')
        % xy center on brightest point
        ym = length(y1:y2);
        xm = length(x1:x2);
        x1 = max(1 ,floor(xi-pars.xyCropWindow+1));
        x2 = min(xm,floor(xi+pars.xyCropWindow));
        y1 = max(1 ,floor(yi-pars.xyCropWindow+1));
        y2 = min(ym,floor(yi+pars.xyCropWindow+1));
        croppedSpots{s} = cropIm(y1:y2,x1:x2,z1:z2);
    else % just crop z
        croppedSpots{s} = cropIm(:,:,z1:z2);
    end
    
    
    % just for troubleshooting
    if pars.troubleshoot
       figure(5); clf; 
       ProjectIm3D(croppedSpots{s},'showPlots',true);
       colormap(gray);
       pause;
    end

    spotXYZ(s,:) = [xi,yi,zi];
end

