function pa = PermAreaRatio(xyz)



binSize = 20;%  sqrt(length(xyz))/2;
binSize = min(binSize,100);
binSize = max(binSize,1);
ymin = min(xyz(:,2));
ymax = max(xyz(:,2));
xmin = min(xyz(:,1));
xmax = max(xyz(:,1));
xBins = xmin:binSize:xmax;
yBins = ymin:binSize:ymax;
mapIdx = hist3([xyz(:,2),xyz(:,1)],{yBins,xBins});  


% mean(mapIdx(:))

bw = imresize(mapIdx>0,4);
perm = bwperim(bw,8);
pa = sum(perm(:))/sqrt(sum(bw(:)));


% figure(2); clf; subplot(1,2,1); imagesc(bw); subplot(1,2,2); imagesc(perm); colormap gray;
% title(num2str(pa,3));


