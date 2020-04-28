function [imLabeled, bw] = SegmentCellsByCentroid(spotXY,imSize,varargin)
% [imLabeled, bw] = SegmentCellsByCentroid(spotXY,imSize,varargin)
% 
% % EXAMPLE: compute properties
% [imLabeled, bw] = SegmentCellsByCentroid(spotXY,imSize)
% figure(3); clf; imagesc(imLabeled); colormap(lines);
% regProps = regionprops(imLabeled,imData(:,:,2),'Area','PixelIdxlist');
% regProps = regionprops(imLabeled,imData(:,:,2),'MaxIntensity','MeanIntensity','MinIntensity','PixelValues');
% [regProps.MaxIntensity]

defaults = cell(0,3);
defaults(end+1,:) = {'cellSize','positive',20};
defaults(end+1,:) = {'imData','freeType',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% ======== segment cells =========
% -------- spot centered segmentation -------
% keepBrightestFraction = 0.75; % probably would have made more sense to click through and get these set up correctly the first time. Maybe this is a good upgrade for the ChrTracer2. 
% 
% [~,idx] = unique(dataTable.uniqueID);
% fidTable = dataTable(idx,:);
% idx = ismember(fidTable.uniqueID,spotTable.uniqueID);
% bright = fidTable.fid_h(idx);
% keep = bright > quantile(bright,keepBrightestFraction);
% xSpots = spotTable.mosaicX(keep);
% ySpots = spotTable.mosaicY(keep);
% in = inpolygon(xSpots,ySpots,polyDat(:,1),polyDat(:,2));
% xSpots = xSpots(in);
% ySpots = ySpots(in);
% % show spots and validate seletion
% figure(3); clf; imagesc(dnaImageCnst); 
% hold on;  colormap gray; axis image;
% plot(xSpots,ySpots,'r.');  


% [h,w,nChns] = size(imData);
h = imSize(1); 
w = imSize(2);
xSpots = spotXY(:,1);
ySpots = spotXY(:,2);
cellMap = false(h,w);
idx = sub2ind([h,w],ySpots,xSpots);
cellMap(idx) = true;
cmap = hsv(length(ySpots));
cmap = cmap(randperm(length(ySpots)),:);
colormap(cmap);

% ~~~~~~~ watersheding the cell map
im1 = bwdist(cellMap);
im1 = watershed(im1);
bw0 = im1==0; % cuts between cells
bw1 = imdilate(cellMap,strel('disk',pars.cellSize)); 
bw1 = bw1 & ~bw0;
bw2 = imdilate(cellMap,strel('disk',pars.cellSize+1));
bw = (bw2 & ~bw1); % borders of cells 
% figure(3); clf; imagesc(bw);
imLabeled = uint16(im1).*uint16(bw1);

% show results
if ~isempty(pars.imData)
    im2 = IncreaseContrast(pars.imData,'low',.8,'high',.999); % (:,:,selChns)
    im3 = im2 + repmat(uint16(bw)*2^16,1,1,size(im2,3));
    Ncolor(im3); axis image; set(gcf,'color','w');
end


