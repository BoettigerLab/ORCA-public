function dataOut = ClusterDistMaps(distMap,varargin)
% dataOut = ClusterDistMaps(distMap)
%
%% Inputs
% distMap   - nReads x nReads x nCells matrix
% 
%% Outputs
% skip      - spots rejected due to too many NaNs (see maxBlankFrac)
% rowTree   - dendogram clustering of distance maps (see 'distanceMetric') 
% rowOrder  - the order of the leaves in the dendogram  
% familyIdx - the cluster assignment of each spot
% mainClusters - cluster assignment of major classes
% newMap    - the moving average filtered map (with skipped spots removed)
% sortMap   - the clustergram-sorted, original distMap
%             (with skipped spots removed)
%
%% Optional Inputs
% defaults(end+1,:) = {'dendogramWalk','integer',0};
% defaults(end+1,:) = {'movingAveBoxWidth','integer',5};
% defaults(end+1,:) = {'maxBlankFrac','fraction',.75};
% defaults(end+1,:) = {'distanceMetric',{'cityblock','minkowski','chebychev','cosine','jaccard','spearman','correlation'},'cityblock'};
% defaults(end+1,:) = {'linkageMetric',{'weighted','average'},'weighted'};
% defaults(end+1,:) = {'maxclust','integer',100};
% defaults(end+1,:) = {'minPerCluster','integer',10};
% defaults(end+1,:) = {'showPlots','boolean',true};
% defaults(end+1,:) = {'showExtraPlots','boolean',false};
% defaults(end+1,:) = {'clim','array',[0,1000]}; % caxis color range
% defaults(end+1,:) = {'verbose','boolean',true};
%
%% Examples
% Examples of processing outputs
% rowTreeHandles = dendrogram(rowTree,nCells);
%
% 

defaults = cell(0,3);
defaults(end+1,:) = {'dendogramWalk','integer',0};
defaults(end+1,:) = {'movingAveBoxWidth','integer',5};
defaults(end+1,:) = {'maxBlankFrac','fraction',.75};
defaults(end+1,:) = {'distanceMetric',{'cityblock','minkowski','chebychev','cosine','jaccard','spearman','correlation'},'cityblock'};
defaults(end+1,:) = {'linkageMetric',{'weighted','average'},'weighted'};
defaults(end+1,:) = {'maxclust','integer',100};
defaults(end+1,:) = {'minPerCluster','integer',10};
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'showExtraPlots','boolean',false};
defaults(end+1,:) = {'clim','array',[0,1000]}; % caxis color range
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);


%% moving box instead of moving 
% this should be its own function
% and it should be rewitten for speed with less looping 
%   (or maybe par for at least for large distance maps?) 

[newMap,skip] = SmoothDistanceMaps(distMap,'parameters',pars);
 
% remove entries that are purely nan
newMap(:,:,skip) = [];


nSpotsFilt = size(newMap,3);

if pars.verbose
    disp('computing clustergram');
end

vectMaps = reshape(newMap,nReads^2,nSpotsFilt);
M = vectMaps';
rowDistMatrix = pdist(M,pars.distanceMetric);      
rowTree = linkage(rowDistMatrix,pars.linkageMetric);
rowTree(isnan(rowTree)) = .000001; % fixing NaNs
rowDistMatrix(isnan(rowDistMatrix)) = 0; % fixing NaNs
rowOrder = optimalleaforder(rowTree,rowDistMatrix);


% Build groups based on clustergram tree
familyIdx = cluster(rowTree,'maxclust',pars.maxclust);
% figure(4); clf; hist(familyIdx,1:pars.maxclust);
v = hist(familyIdx,1:pars.maxclust);
keepClusts = v > pars.minPerCluster;
familyIDs = zeros(size(familyIdx));
k = 0;
for n=find(keepClusts)
    k=k+1;
    familyIDs(familyIdx==n) = k;
end

sortMap = distMap;
sortMap(:,:,skip) = [];

dataOut.familyIdx = familyIdx;
dataOut.mainClusters = familyIDs;
dataOut.rowTree = rowTree;
dataOut.rowOrder = rowOrder;
dataOut.skip = skip;
dataOut.newMap = newMap;
dataOut.sortMap = sortMap;

if pars.showPlots
% show major classes fof distance maps
    nClust = k;
    p = ceil(sqrt(nClust));
    for n=1:nClust
        subplot(p,p,n); imagesc( nanmedian( sortMap(:,:,familyIDs==n),3) );
        caxis(pars.clim);
        title(sum( familyIDs==n));
    end
    GetColorMap('RedWhiteBlueSat');
end
  


% show plots, (mostly for troubleshooting)
if pars.showExtraPlots
    figure(1); clf; % averaged images, not yet sorted
    labelNums = cellstr( num2str( (1:nSpotsFilt)' ) );
    [im,lO]= TileImage(newMap,'multicolor',false,'numRows',10,'labelNames',labelNums);
    imagesc(im); hold on; 
    text(lO(:,1),lO(:,2),labelNums,'y');
    colorbar;  
    caxis(pars.clim);
    GetColorMap('RedWhiteBlueSat');
    title('unsorted smoothed maps');

    figure(2); clf;  % average images
    labelNums = cellstr( num2str( (1:nSpotsFilt)' ) );
    [im,lO]= TileImage(newMap(:,:,rowOrder),'multicolor',false,'numRows',30,'labelNames',labelNums);
    imagesc(im); hold on; text(lO(:,1),lO(:,2),labelNums(rowOrder),'y');
    size(im)
    colorbar;  caxis(pars.clim);
    GetColorMap('RedWhiteBlueSat');
    title('sorted smoothed maps');

    figure(3); clf; 
    dendrogram(rowTree,nSpotsFilt);

    figure(4); clf; 
    labelNums = cellstr( num2str( (1:nSpotsFilt)' ) );
    [im,lO]= TileImage(sortMap(:,:,rowOrder),'multicolor',false,'numRows',10,'labelNames',labelNums);
    imagesc(im); hold on; 
    text(lO(:,1),lO(:,2),labelNums(rowOrder),'y');
    size(im)
    colorbar;  
    caxis(pars.clim);
    GetColorMap('RedWhiteBlueSat');
    title('original dist-maps, sorted')
end


%% walk down dendogram

for maxclust = 2:pars.dendogramWalk
    familyIdx = cluster(rowTree,'maxclust',maxclust);
    figure(5); clf; hist(familyIdx,1:maxclust);

    p = ceil(sqrt(maxclust));
    figure(5); clf;
    for m=1:maxclust
        subplot(p,p,m); 
        imagesc( nanmedian( sortMap(:,:,familyIdx==m),3) );
        caxis(pars.clim);
        title(sum( familyIdx==m));
    end
    GetColorMap('RedWhiteBlueSat');
    pause(1);
end

