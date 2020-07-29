function [clusterGene, grp,familyIdx] = PlotClusterGenes(aveCorrMatrix,libRealGenes,varargin)


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'desiredClusters', 'nonnegative', 10}; 
defaults(end+1,:) = {'minLargestCluster', 'nonnegative', 20}; 
defaults(end+1,:) = {'maxClusterSize', 'nonnegative', 0};
defaults(end+1,:) = {'numClusters', 'nonnegative', 0}; 
defaults(end+1,:) = {'minClusterSize', 'nonnegative', 3}; 
defaults(end+1,:) = {'maxClusterSize', 'nonnegative', true};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);



% -------------------------------------------------------------------------
% Main Figure
% -------------------------------------------------------------------------

% Plot Clustergram and define groups
clusterGene = CreateClustergram(aveCorrMatrix,...
    'RowLabels',libRealGenes,'ColumnLabels',libRealGenes,'clustObj',0,...
    'columnClusters',[],'linkageMethod','average','distance','correlation');
caxis([-.4 ,.5]);

%----------------------------------------------------------------------

% Computationally group genes based on tree
desiredClusters = parameters.desiredClusters; 
minLargestCluster = parameters.minLargestCluster;
maxClusterSize = parameters.maxClusterSize;
numClusters = parameters.numClusters;
minClusterSize = parameters.minClusterSize;

numRealGenes = length(libRealGenes);

i=0;
while numClusters < desiredClusters || maxClusterSize < minLargestCluster
    i=i+1;
    clustThresh = numRealGenes-i;
    familyIdx = cluster(clusterGene.rowTree,'maxclust',clustThresh);
    [uniqueVals,numOcc,idx] = occurrences(familyIdx);
    numOccurs = numOcc(idx); 
    medianClustSize = median(numOccurs(numOccurs>= minClusterSize));
    maxClusterSize = max(numOccurs(numOccurs>= minClusterSize));
    
    familyIdx(numOccurs< minClusterSize) = 0;
    nonSingleV = uniqueVals(numOcc >= minClusterSize);
    for j=1:numel(nonSingleV)
       familyIdx(familyIdx==nonSingleV(j)) = j; 
    end
    numClusters = max(familyIdx); 
end
numGroups = max(familyIdx);


% ----------- sort groups linearly along the tree
sortedGenes = libRealGenes(clusterGene.rowOrder); 
groupPos = zeros(numGroups,1);
for g=1:numGroups
    grpCoords = StringFind(sortedGenes,libRealGenes(familyIdx==g),'exactly',1);
    groupPos(g) = mean(grpCoords);
end
[jnk,grpSort] = sort(groupPos);
groupIdx = zeros(numRealGenes,1);
k=0;
for g=grpSort'
    k=k+1;
    groupIdx(familyIdx==g) = k;
end
familyIdx = groupIdx;
familyIdx = familyIdx+1;
numGroups = max(familyIdx);
% -----------------------------------------------------------------------
    

% -------------- color trees
clrMap = hsv(numGroups-1); 
% clrMap = clrMap( randperm(numGroups-1),:);
ColorDendogram(clusterGene,ones(numRealGenes,1),'colormap',.85*ones(1,3));
ColorDendogram(clusterGene,familyIdx-1,'colormap',clrMap);

%--------------- Add color boxes
grp = cell(numGroups,1);
for g=2:numGroups
    grp{g} = libRealGenes(familyIdx==g);
    grpIdx = StringFind(sortedGenes,grp{g},'exactly',1);
    corner = min(grpIdx)-1;
    rectangle('Position',[corner+.5,corner+.5,length(grpIdx),length(grpIdx)],'EdgeColor','k','lineWidth',.5);
end
grp{1} = libRealGenes(familyIdx==1); % ungrouped
