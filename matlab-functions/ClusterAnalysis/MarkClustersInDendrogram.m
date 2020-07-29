function MarkClustersInDendrogram(treeHandles,linkageMap, familyIdx, plotProperties, varargin)
% ColorDendogram(clusterGene,familyIdx,'colormap',myColorMap)
% 


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', false};
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

% -------------------------------------------------------------------------
% Determine number of genes, groups, and tree indices
% -------------------------------------------------------------------------
numGenes = length(familyIdx);
numGroups = max(familyIdx);
treeIdx = [linkageMap, (numGenes+(1:numGenes-1))'];

% -------------------------------------------------------------------------
% Color based on groups
% -------------------------------------------------------------------------
for i=1:numGroups
    leafIdx = find(familyIdx == i);
    [~,leafHandles] =intersect( treeIdx(:,1), leafIdx );
    
    set(treeHandles(leafHandles), plotProperties{i}{:}); % color the leaves

    ancestHandles = leafHandles;
    while ~isempty(ancestHandles);
        ancestIdx = treeIdx(ancestHandles,4);
        [~,ancestHandles] = intersect(treeIdx(:,1),ancestIdx);
       set(treeHandles(ancestHandles), plotProperties{i}{:}); % color the parents
    end
end