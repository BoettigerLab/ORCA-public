function [vlist,clustData] = ClusterStats(vlist,varargin)
% vlist = ClusterStats(vlist) - adds clusterCount, clusterArea, and
%                clusterDiameter fields to vlist
% vlist = ClusterStats(vlist,'boxSize',boxSize,'minLoc',minLoc,'showMap',t)              

%--------------------------------------------------------------------------
% Default optional variables
%--------------------------------------------------------------------------
boxSize = 32; % in nm
npp = 160;  % nm per pixel
minLoc = 1; % min # of localization per box to count a box
showMap = true;
showHist = true; 
normBoxSize = 0; % density dependent box-size normalization
normMinLoc = 0; % density dependent min localizations per box
%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'boxSize'
                boxSize = CheckParameter(parameterValue,'positive','boxSize'); 
            case 'npp'
                npp = CheckParameter(parameterValue,'positive','npp'); 
            case 'minLoc'
                minLoc = CheckParameter(parameterValue,'nonnegative','minLoc'); 
            case 'showMap'
                showMap = CheckParameter(parameterValue,'boolean','showMap'); 
            case 'showHist'
                showHist = CheckParameter(parameterValue,'boolean','showHist'); 
            case 'normBoxSize'
                normBoxSize = CheckParameter(parameterValue,'nonnegative','normBoxSize'); 
            case 'normMinLoc'
                normMinLoc = CheckParameter(parameterValue,'nonnegative','normMinLoc'); 
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
%% Main code
%--------------------------------------------------------------------------


x = vlist.xc*npp;
y = vlist.yc*npp;  
H = max(y)+1;
W = max(x)+1;

density = length(x)/(H);
if normBoxSize > 0
  boxSize = round(boxSize*sqrt(normBoxSize/length(x)) );
%  boxSize = round(boxSize*sqrt(1/density));
  boxSize = max(10,boxSize); 
end
if normMinLoc > 0
    minLoc = minLoc*(density/normMinLoc);
end
if normMinLoc > 0 || normBoxSize > 0
    disp([boxSize,minLoc]);
end
 
xbins = 0:boxSize:W;
ybins = 0:boxSize:H;
clustMap = hist3([y,x],'Edges',{xbins,ybins});
[h,w] = size(clustMap);


if sum(clustMap) == 0
    warning('no spots, trying centers method');
    clustMap = hist3([y,x],{xbins,ybins});
end

bwMap = clustMap >= minLoc;    
props = regionprops(bwMap,clustMap,'Area','PixelValues','PixelIdxList','Centroid');
locsPerPixel = {props.PixelValues}; % cell array of length number clusters. Each entry contains a vector recording the # of localizations in each pixel
locsPerCluster = cellfun(@sum,locsPerPixel); % total number of localizations in each cluster
% 
% % % The original method 
% xLocsPerCluster = linspace(0,max(locsPerCluster),max(locsPerCluster));
% allLocs = hist(locsPerCluster,xLocsPerCluster);
% clusterCountsPerLoc = allLocs.*xLocsPerCluster; % simply multiply the number of clusters by the cluster size; 
% figure(2); clf; bar(xLocsPerCluster,clusterCountsPerLoc);
% numSizes = length(xLocsPerCluster);  % then merge clusters into bins
% stp = 25; 
% xBins = 0:stp:max(locsPerCluster);
% % xBins = logspace(0,log2(max(locsPerCluster)),50);
% binnedLocs = zeros(length(xBins),1);
% for i=2:length(xBins)
%     inRange = clusterCountsPerLoc > xBins(i-1) & clusterCountsPerLoc < xBins(i);
%     binnedLocs(i) = sum( xLocsPerCluster( inRange) );
% end
% figure(3); clf; bar(xBins,binnedLocs);
% figure(4); clf; bar(log2(xBins),binnedLocs); xlim([0,1+log2(max(locsPerCluster))]);
% % binnedClusterCountsPerLoc = sum(reshape(clusterCountsPerLoc',stp,numSizes/stp));
% % figure(2); clf; bar(xLocsPerCluster(1:stp:end),binnedClusterCountsPerLoc);
% % xlabel('x (cluster size, in counts)');
% % ylabel('number of localizations belonging to clusters of size x');

% The new method
allIdx = cat(1,props.PixelIdxList); % the linear indices of all the active pixels  
pixelClustId = {props.PixelIdxList};
clustNums = num2cell(1:length(props)); % the cluster ID number of each cluster, rendered as a cell array  of length # clusters
pixelClustId = cellfun(@(x,y) y*ones(length(x),1), pixelClustId,clustNums,'UniformOutput',false); 
pixelClustId = cat(1,pixelClustId{:}); % the cluster ID numbers for every pixel containing a localizations of length # active pixels
allArea = [props.Area]';


clustData.Area = allArea*boxSize^2;
clustData.Diameter = 2*sqrt(allArea*boxSize^2/pi);
clustData.Counts = Column(locsPerCluster);
clustData.Centroid = cat(1,props.Centroid); 

% Put it in a table; 
% [pixelClustId,allIdx,allArea(pixelClustId),locsPerCluster(pixelClustId)];

allAreaPerPixel = allArea(pixelClustId)*boxSize^2;
allLocsPerPixel = Column(locsPerCluster(pixelClustId)); % every pixel now knows how many localizations are its cluster 
% figure(2); clf; hist(allAreaPerPixel,50);
% figure(3); clf; hist(locsPerCluster(pixelClustId),50);

% For each molecule in the vlist, record the area and total counts of the
% associated cluster.  
 binId = sub2indFast([h,w],round((y+boxSize/2)/boxSize),round((x+boxSize/2)/boxSize)); % every molecule knows the linear index of its pixel
 % binId2 = sub2ind([h,w],round((y+boxSize/2)/boxSize),round((x+boxSize/2)/boxSize)); % every molecule knows the linear index of its pixel
[~,b] = ismember(binId,allIdx); % intersect linear indices of molecules with linear indices of active pixels  

% because of our filters, some localizations didn't get assigned to clusters.  
% we'll report these molecules as cluster size 0.  
allLocsPerPixel = [allLocsPerPixel;0]; % cluster size 0 means not-assigned to cluster
allAreaPerPixel = [allAreaPerPixel;0];
b(b==0) = length(allLocsPerPixel); 


vlist.clustArea = allAreaPerPixel(b);
vlist.clustDiameter = 2*sqrt(allAreaPerPixel(b)/pi);
vlist.clustCount = allLocsPerPixel(b);




if showMap
    figure(1); clf; 
    imagesc(clustMap); colorbar;
    
    try
    clustMapSuper = clustMap;
    clustMapSuper(binId)  = 20;
    missedPixels = setdiff(binId,allIdx);
    clustMapSuper(missedPixels) = 10;
    figure(3); clf; imagesc(clustMapSuper); caxis([0,20]);
    catch 
        figure(3); clf; imagesc(clustMap>1); 
    end
        
end

if showHist
    figure(2); clf;
    subplot(1,2,1); hist(vlist.clustCount,20);
    xlabel('cluster count (# localizations)'); 
    subplot(1,2,2); hist(vlist.clustDiameter,20);
    xlabel('cluster diameter (nm)'); 
end
% Confirm binning is done correct
% map2 = false(h,w);
% map2(binId) = true;
% figure(3); clf; imagesc(map2);
% figure(3); clf; imagesc(map2 - bwMap);
% 
% 
% figure(2); clf; plot(locsPerCluster,[props.Area],'k.');
% figure(2); clf; hist(locsPerCluster,50);
% figure(2); clf; hist([props.Area],50);

    