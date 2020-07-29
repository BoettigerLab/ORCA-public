function clusterObj = CreateClustergram(M,varargin)
% 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'ColumnLabels', 'cell', {}};
defaults(end+1,:) = {'RowLabels', 'cell', {}};
defaults(end+1,:) = {'clustObj', 'boolean', false};
defaults(end+1,:) = {'clustFig', 'boolean', true};
defaults(end+1,:) = {'columnClusters', 'nonnegative', []};
defaults(end+1,:) = {'rowClusters', 'nonnegative', []};
defaults(end+1,:) = {'rowGene', 'boolean', false};
defaults(end+1,:) = {'linkageMethod', {'complete','average','single'}, 'average'};
defaults(end+1,:) = {'distance', {'euclidean','seuclidean','correlation','cityblock','mikowski','chebychev','mahalanobis','cosine','spearman','hamming','jaccard'}, 'correlation'};
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

%%
clusterObj = [];

if parameters.clustFig
    [numGenes,numCells] = size(M); 
    if isempty(parameters.columnClusters)
        parameters.columnClusters = numCells;
    end
    if isempty(parameters.rowClusters)
        parameters.rowClusters = numGenes;
    end
    
    rowDistMatrix = pdist(M,parameters.distance);       
    rowTree = linkage(rowDistMatrix,parameters.linkageMethod);
    rowOrder = optimalleaforder(rowTree,rowDistMatrix);
    M1 = M(rowOrder,:);
    
    if ~parameters.rowGene
        columnDistMatrix = pdist(M',parameters.distance);
        columnTree = linkage(columnDistMatrix,parameters.linkageMethod);
        columnOrder = optimalleaforder(columnTree,columnDistMatrix);
        M2 = M1(:,columnOrder);
    else
        M2 = M1(:,rowOrder);
        columnTree = rowTree;
        columnOrder = rowOrder;
    end
    
    rowAxes = subplot(2,2,3); 
    rowTreeHandles = dendrogram(rowTree,parameters.rowClusters,'Reorder',rowOrder,'Orientation','left');
    columnAxes =subplot(2,2,2); 
    columnTreeHandles = dendrogram(columnTree,parameters.columnClusters,'Reorder',columnOrder);
    mapAxes = subplot(2,2,4);
   

    imagesc(M2);
    if ~isempty(parameters.RowLabels);
        set(gca,'YTick',1:numGenes,'YTickLabel',parameters.RowLabels(rowOrder),'Ticklength',[0 0])
    end
    if ~isempty(parameters.ColumnLabels)
        set(gca,'XTick',1:numCells,'XTickLabel',parameters.ColumnLabels(columnOrder),'Ticklength',[0 0]);
    end
    set(gca,'YDir','normal');
    linkaxes([mapAxes,columnAxes],'x'); hold on;
    linkaxes([mapAxes,rowAxes],'y');
%     linkprop([columnAxes, mapAxes],{'XLim'});
%     linkprop([rowAxes, mapAxes],{'YLim'});
    
    PresentationPlot('FontSize',5,'LineWidth',.5);
    set(rowAxes,'Position',[.05 .05 .2 .7],'YTick',[]);
    set(columnAxes,'Position',[.3 .77 .65 .15],'XTick',[]);
    set(mapAxes,'Position',[.3 .05 .65 .7]); 
    % caxis([-.4,.4]);
    
    clusterObj.rowOrder = rowOrder;
    clusterObj.columnOrder = columnOrder; 
    clusterObj.rowDistMatrix = rowDistMatrix;
    clusterObj.columnDistMatrix = columnDistMatrix;
    clusterObj.rowTree = rowTree;
    clusterObj.columnTree = columnTree;
    clusterObj.rowTreeHandles = rowTreeHandles;
    clusterObj.columnTreeHandles = columnTreeHandles;
    clusterObj.columnAxes = columnAxes;
    clusterObj.rowAxes = rowAxes;
    clusterObj.mapAxes = mapAxes; 
    clusterObj.map = M2;  
    if ~isempty(parameters.ColumnLabels);
        clusterObj.columnLabels = parameters.ColumnLabels(columnOrder);
    else
        clusterObj.columnLabels = parameters.ColumnLabels;
    end
    if ~isempty(parameters.RowLabels);
        clusterObj.rowLabels = parameters.RowLabels(rowOrder);
          else
        clusterObj.rowLabels = parameters.RowLabels;
    end
end

if parameters.clustObj
    clustObj = clustergram(M,'RowPDist','correlation','ColumnPDist','correlation','colormap','jet',...
        'RowLabels',parameters.RowLabels,'ColumnLabels',parameters.ColumnLabels);
    if isempty(clusterObj)
        clusterObj = clustObj;
    else
        clusterObj.matlabClusterObj = clustObj;
    end
end


