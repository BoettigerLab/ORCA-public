function clusterObj = MyClustergram(M,varargin)
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
defaults(end+1,:) = {'geneGene', 'boolean', false};
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
    
    geneDistMatrix = pdist(M,'correlation');       
    geneTree = linkage(geneDistMatrix,'average');
    try
        geneOrder = optimalleaforder(geneTree,geneDistMatrix);
    catch er
        disp(er.message);
        disp('surpressing nans');
        geneTree(isnan(geneTree)) = 0;
        geneDistMatrix(isnan(geneDistMatrix)) = 0;
        geneOrder = optimalleaforder(geneTree,geneDistMatrix);
    end
    M1 = M(geneOrder,:);
    
    if ~parameters.geneGene
        cellDistMatrix = pdist(M','correlation');
        cellTree = linkage(cellDistMatrix,'average');
        cellOrder = optimalleaforder(cellTree,cellDistMatrix);
        M2 = M1(:,cellOrder);
    else
        M2 = M1(:,geneOrder);
        cellTree = geneTree;
        cellOrder = geneOrder;
    end
    
    genes1 = subplot(2,2,3); 
    geneTreeHandles = dendrogram(geneTree,parameters.rowClusters,'Reorder',rot90(geneOrder,2),'Orientation','left');
    cells2 =subplot(2,2,2); 
    cellTreeHandles = dendrogram(cellTree,parameters.columnClusters,'Reorder',rot90(cellOrder,2));
    ax3 = subplot(2,2,4);
   

    imagesc(M2);
    if ~isempty(parameters.RowLabels);
        set(gca,'YTick',1:numGenes,'YTickLabel',parameters.RowLabels(geneOrder),'Ticklength',[0 0])
    end
    if ~isempty(parameters.ColumnLabels)
        set(gca,'XTick',1:numCells,'XTickLabel',parameters.ColumnLabels(cellOrder),'Ticklength',[0 0]);
    end
    linkaxes([cells2,ax3],'x'); hold on;
    linkaxes([genes1,ax3],'y');
    
    
    
%     PresentationPlot('FontSize',5);
    set(genes1,'Position',[.05 .05 .2 .7],'YTick',[]);
    set(cells2,'Position',[.3 .77 .65 .15],'XTick',[]);
    set(ax3,'Position',[.3 .05 .65 .7]); 
    % caxis([-.4,.4]);
    
    clusterObj.geneOrder = geneOrder;
    clusterObj.cellOrder = cellOrder; 
    clusterObj.geneDistMatrix = geneDistMatrix;
    clusterObj.cellDistMatrix = cellDistMatrix;
    clusterObj.geneTree = geneTree;
    clusterObj.cellTree = cellTree;
    clusterObj.geneTreeHandles = geneTreeHandles;
    clusterObj.cellTreeHandles = cellTreeHandles;
    clusterObj.cellAxes = cells2;
    clusterObj.geneAxes = genes1;
    clusterObj.mapAxes = ax3; 
    clusterObj.map = M2;  
    if ~isempty(parameters.ColumnLabels);
        clusterObj.cellLabels = parameters.ColumnLabels(cellOrder);
    else
        clusterObj.cellLabels = parameters.ColumnLabels;
    end
    if ~isempty(parameters.RowLabels);
        clusterObj.geneLabels = parameters.RowLabels(geneOrder);
          else
        clusterObj.geneLabels = parameters.RowLabels;
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


