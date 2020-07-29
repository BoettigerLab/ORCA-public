function oa = OpenAreaRatio(xyz,varargin)

% parameters
binSize = 40;
zbinSize = 80;
showImage = false;

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
bwA = imdilate(bw,strel('disk',1)) ;
sA = bwA - bw;

oa = sum(sA(:)) ./ sum(bw(:))^(2/3) ;

if showImage
    vprops = regionprops(bw,'PixelList');
    vpix = cat(1,vprops(:).PixelList);

    aprops = regionprops(sA,'PixelList');
    apix = cat(1,aprops(:).PixelList);

    figure(4); clf;
    plot3(vpix(:,1),vpix(:,2),vpix(:,3),'k.');
    hold on;  plot3(apix(:,1),apix(:,2),apix(:,3),'r.');
end
% 
% figure(1); clf; imagesc(sum(bw,3));
% figure(2); clf; imagesc(sum(bwA,3));
% figure(3); clf; imagesc(sum(sA,3));
% 
