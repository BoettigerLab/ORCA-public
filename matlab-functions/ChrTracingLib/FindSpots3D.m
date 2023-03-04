function xyz = FindSpots3D(imIn,varargin)
% input, imIn = a 3D image, may contain multiple spots
% output, xyz, a Nx3 matrix of the nearest-pixel centroid coordinates 
% 
% Notes
% only xy is downsampled for speed and simplicity
% should add in future optional alternative downsampling for z. 

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'showPlot','boolean',false};

% parse variable arguments
pars = ParseVariableArguments(varargin, defaults, mfilename);

% main function
imDS = imresize(imIn,1/pars.autoSelectDownsample); 
bw = imregionalmax(imDS);  
peakVals = imDS(bw);
bMin = quantile(peakVals,pars.autoSelectThreshold);
bw(imDS< bMin) = 0;
[y,x,z] = ind2sub(size(bw),find(bw));
xyz = [round([x,y]*pars.autoSelectDownsample),z];

if pars.showPlot || nargout==0
    ProjectIm3D(imIn,'showPlots',true,'caxis',[0,5*bMin]);
    subplot(1,3,1); hold on; plot(xyz(:,1),xyz(:,2),'ro'); title('xy')
    subplot(1,3,2); hold on; plot(xyz(:,2),xyz(:,3),'ro'); title('yz')
    subplot(1,3,3); hold on; plot(xyz(:,1),xyz(:,3),'ro'); title('xz')
end