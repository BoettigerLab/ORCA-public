function ColorCodeDendrogram(rowTree,familyIdx,rowTreeHandles,varargin)
% ColorCodeDendrogram(rowTree,familyIdx,rowTreeHandles,'colormap',myColorMap)
% 
% 'colormap' - a colormap of length the number of groups


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'colormap', 'colormap', []};
defaults(end+1,:) = {'linewidth', 'nonnegative', 2};
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'clusterGene structure,familyIdx vector are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);

%% 
numGenes = length(familyIdx); 
% columnTreeHandles = clusterGene.columnTreeHandles; % 
numGroups = max(familyIdx);
treeIdx = [rowTree, (numGenes+(1:numGenes-1))'];

if isempty(parameters.colormap)
    clrs = jet(numGroups);
    clrs = clrs(randperm(numGroups),:);
else
    clrs = parameters.colormap;
end
    
% % reset
% clrs = [0,0,1]; numGroups = 1; familyIdx = ones(numGenes,1); 
for i=1:numGroups
    leafIdx = find(familyIdx == i);
    [~,leafHandles] =intersect( treeIdx(:,1), leafIdx );
    set(rowTreeHandles(leafHandles),'color',clrs(i,:),'linewidth',parameters.linewidth); % color the leaves
    % set(columnTreeHandles(leafHandles),'color',clrs(i,:),'linewidth',parameters.linewidth); % color the leaves

    ancestHandles = leafHandles;
    while ~isempty(ancestHandles);
        ancestIdx = treeIdx(ancestHandles,4);
        [~,ancestHandles] = intersect(treeIdx(:,1),ancestIdx);
       set(rowTreeHandles(ancestHandles),'color',clrs(i,:),'linewidth',parameters.linewidth); % color the parents
       % set(columnTreeHandles(ancestHandles),'color',clrs(i,:),'linewidth',parameters.linewidth); % color the parents
    end

end