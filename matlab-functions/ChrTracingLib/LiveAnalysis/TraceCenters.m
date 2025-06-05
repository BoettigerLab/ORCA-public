function [xy0,ptMap1,cnts] = TraceCenters(xs1,ys1,varargin)
% xs1  cell array of length nFrames of x coords
%   UPDATE this to use a table or a structure for more elegance. 
% 
%  imageSize w x h
% 
% see also: LinkLivePoints

defaults = cell(0,3);
defaults(end+1,:) = {'frameFraction','fraction',1}; % fraction of frames to use
defaults(end+1,:) = {'binResolution','positive',4}; % 
defaults(end+1,:) = {'showPlots','boolean',0}; % 
defaults(end+1,:) = {'minSpotsPerBin','integer',20}; % 
defaults(end+1,:) = {'minSpotsPerTrace','integer',500}; % 
defaults(end+1,:) = {'imSize','array',[]}; % 
pars = ParseVariableArguments(varargin,defaults,mfilename);

% % obsolete pars
% defaults(end+1,:) = {'autoSelectThreshold','fraction',.995}; % 
% defaults(end+1,:) = {'autoSelectDownsample','positive',1}; % 


nFrames = length(xs1);
maxFrame = round(nFrames*pars.frameFraction);
binResolution = pars.binResolution;
if isempty(pars.imSize)
w = round(1.1*max(cat(1,xs1{:})));
h = round(1.1*max(cat(1,ys1{:})));
else
    h=pars.imSize(1);
    w=pars.imSize(2); 
end
all1 = [cat(1,ys1{1:maxFrame}),cat(1,xs1{1:maxFrame})]+1;
ptMap1 = hist3(all1,'Ctrs',{linspace(1,h,h/binResolution),linspace(1,w,w/binResolution)}); % 
% I think we could do better here with a region props and a minimum number
% of localizations per bin filter rather than that autoSelectThreshold

bwMap  = ptMap1 >= pars.minSpotsPerBin;
regs = regionprops(bwMap,ptMap1,'Centroid','Area','PixelValues');
pv = {regs.PixelValues};
cnts = cellfun(@sum,pv);
% figure(3); clf; plot(sort(cnts))
spots1 = cat(1,regs.Centroid);
spots1 = spots1(cnts>=pars.minSpotsPerTrace,:);

% spots1 = FindSpots(uint16(ptMap1),'autoSelectThreshold',pars.autoSelectThreshold,'autoSelectDownsample',pars.autoSelectDownsample);
if pars.showPlots > 0 
    try
    figure(pars.showPlots); clf; imagesc(imresize(ptMap1,binResolution,'nearest')); 
    colorbar;  clim([floor(.9*pars.minSpotsPerBin),quantile(ptMap1(:),.9999)])
    hold on; plot(spots1(:,1)*binResolution,spots1(:,2)*binResolution,'ro'); colormap('parula');
    catch 
    end
end
xy0 = (spots1-1)*binResolution;
