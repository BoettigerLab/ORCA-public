function PlotSpotMatrices(savePath,spotMatrices,libGenes,varargin)
% Plotting SpotMatices (all cells combined).  


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'flag', 'string', ''};
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 3
    error('matlabSTORM:invalidArguments', 'requires savePath,spotMatrices,libGenes');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
flag = parameters.flag; 
%%
% flag = ''; 

numGenes = size(libGenes,1); 
allGenes = 1:numGenes;
negIdx = [StringFind(libGenes,'blank');StringFind(libGenes,'notarget')];
libGenes(negIdx)
realGenes = allGenes; 
realGenes(negIdx) = [];
cntrlGenes = allGenes(negIdx); 
numRealGenes = length(realGenes);
numControls = length(cntrlGenes); 


% Plot spot matrices for top expressed genes
try
k = find( strcmp(libGenes,'THBS1') );
THBS1matrix = figure; clf;  
imagesc(spotMatrices{k}); 
PresentationPlot; colorbar; colormap(hot(256));
xlabel('mRNA #'); ylabel('hybe #'); title(libGenes{k});
export_fig(THBS1matrix,[savePath,'THBS1matrix_',flag,'.png']);
catch
end

% Plot spot matrices for top expressed genes
try
k = find( strcmp(libGenes,'TEAD1') );
TEAD1matrix = figure; clf;  
imagesc(spotMatrices{k}); 
PresentationPlot; colorbar; colormap(hot(256));
xlabel('mRNA #'); ylabel('hybe #'); title(libGenes{k});
export_fig(TEAD1matrix,[savePath,'TEAD1matrix_',flag,'.png']);
catch
end

try
k = find( strcmp(libGenes,'FLNA') );
FLNAmatrix = figure; clf;
imagesc(spotMatrices{k});
PresentationPlot; colorbar; colormap(hot(256));
xlabel('mRNA #'); ylabel('hybe #'); title(libGenes{k});
export_fig(FLNAmatrix,[savePath,'FLNAmatrix_',flag,'.png']);
catch
end


% % Combine all spots to one big matirx
allSpotMatrixFig = figure; clf;
geneCnts = cellfun(@length,spotMatrices(realGenes));
geneLabels = zeros(1,sum(geneCnts));
geneCnts = [1,geneCnts'];
for n=1:numRealGenes
    geneLabels(sum(geneCnts(1:n)):sum(geneCnts(1:n+1))) = n;
end

subplot(4,1,1);
imagesc(geneLabels); 
colormap(lines(numRealGenes));
% cbar = colorbar;  
title('all genes'); 
i = 0;
for n=realGenes;
    i=i+1;
    if i<max(realGenes)
        xposition = mean( [sum(geneCnts(1:i)), sum(geneCnts(1:i+1)) ] );
    else
        xposition = sum(geneCnts(1:i));
    end
    tLabel = text(xposition,1,libGenes{n},'HorizontalAlignment','center','FontSize',16);  % -.4*rand,libGenes{n
    set(tLabel,'rotation',45)
end
freezeColors; % cbfreeze(cbar);
set(gca,'Xtick',[],'Ytick',[],'XtickLabel',{},'YTickLabel',{});

subplot(4,1,2:4);
imagesc(cat(2,spotMatrices{realGenes}));
PresentationPlot; % cbar = colorbar; 
colormap(hot(256));
xlabel('mRNA #'); ylabel('hybe #');
freezeColors; % cbfreeze(cbar);
set(gcf, 'Units','Inches','Position',[0,0,22,6],'PaperUnits', 'Inches', 'PaperSize', [22 6]);
% close all;

cntrlSpotMatrixFig = figure; clf;
geneCnts = cellfun(@length,spotMatrices(cntrlGenes));
geneLabels = zeros(1,sum(geneCnts));
geneCnts = [1,geneCnts'];
for n=1:numControls
    geneLabels(sum(geneCnts(1:n)):sum(geneCnts(1:n+1))) = n;
end
subplot(4,1,1); imagesc(geneLabels);
colormap(lines(numControls));
% colorbar; 
title('not target and blank genes');
i=0;
for n=cntrlGenes % n=4; i = 5;
    i=i+1;
    if i<max(cntrlGenes)
        xposition = mean( [sum(geneCnts(1:i)), sum(geneCnts(1:i+1)) ] );
    else
        xposition = sum(geneCnts);
    end
    tLabel = text(xposition,1,libGenes{n},'HorizontalAlignment','center','FontSize',16); 
    set(tLabel,'rotation',45)
end
set(gca,'Xtick',[],'Ytick',[],'XtickLabel',{},'YTickLabel',{});
freezeColors;

subplot(4,1,2:4); % subplot(8,1,6:8);
imagesc(cat(2,spotMatrices{cntrlGenes}));
PresentationPlot; % cbar = colorbar;  
colormap(hot(256));
xlabel('mRNA #'); ylabel('hybe #'); 
set(gcf, 'Units','Inches','Position',[0,0,22,6],'PaperUnits', 'Inches', 'PaperSize', [22 6]);
freezeColors;%  cbfreeze(cbar); 
spaceplots; 
export_fig(allSpotMatrixFig,[savePath,'allSpotMatrix_',flag,'.png']);
export_fig(cntrlSpotMatrixFig,[savePath,'cntrlSpotMatrix_',flag,'.png']);
