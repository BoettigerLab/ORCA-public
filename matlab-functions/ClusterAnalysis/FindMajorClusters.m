function familyIdx = FindMajorClusters(tree,varargin)
% 
% Example:
% familyIdx  = FindMajorClusters(tree,'desiredClusters',20,...
%              'minLargestCluster',20,'minClusterSize',3);
% Will return at least 20 clusters, all with at least 3 members, the 
% largest of which has at least 20 members.  familyIdx is length numGenes, 
% with the entry for each gene recording the ID number of the cluster it 
% got assigned to.


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'desiredClusters', 'positive', 20};
defaults(end+1,:) = {'minLargestCluster', 'positive', 20};
defaults(end+1,:) = {'minClusterSize', 'positive', 3};
% --------------------------------------------------------

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);

%% Computationally group genes based on tree

numGenes = size(tree,1); 

maxClusterSize = 0;
numClusters = 0;
minClusterSize = parameters.minClusterSize;

i=0;
while numClusters < parameters.desiredClusters || maxClusterSize < parameters.minLargestCluster
    i=i+1;
    clustThresh = numGenes-i;
    familyIdx = cluster(tree,'maxclust',clustThresh);
    [uniqueVals,numOcc,idx] = occurrences(familyIdx);
    numOccurs = numOcc(idx); 
    numClusters = sum(numOccurs>= minClusterSize) ;
    maxClusterSize = max(numOccurs(numOccurs>= minClusterSize));
    
    familyIdx(numOccurs< minClusterSize) = 0;
    nonSingleV = uniqueVals(numOcc >= minClusterSize);
    for j=1:numel(nonSingleV)
       familyIdx(familyIdx==nonSingleV(j)) = j; 
    end
    
%     figure(2); clf; hist(familyIdx,1:numel(nonSingleV));
%     figure(3); clf; plot(n);
%     pause(.1);
end
