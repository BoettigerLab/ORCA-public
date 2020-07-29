function geneData = ClusterPixels(wordMap,alignedIm,libGenes,c,varargin) 
% Child Function of DecodePixels (step 3)

% ------------------------------------------------------------------------
% Define default values
% ------------------------------------------------------------------------
global figureSavePath;

defaults = cell(0,3);
defaults(end+1,:) = {'minPixels', 'nonnegative', 8}; %  % down scale this if we drop upsampling.
defaults(end+1,:) = {'calcSpotMat', 'boolean', false}; %  
defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; % 
defaults(end+1,:) = {'moreVerbose', 'boolean', true}; % 
defaults(end+1,:) = {'overwrite', 'boolean', true}; % 
defaults(end+1,:) = {'savePath', 'string',figureSavePath}; % 
defaults(end+1,:) = {'exportSpatialImage', 'boolean', false}; % 
defaults(end+1,:) = {'zoomBox', 'array', []}; 
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% get some variables
[yDim,xDim,numBits] = size(alignedIm);
numGenes = length(libGenes);

%% Main Function
geneMap = 1+numGenes*ones(yDim,xDim);
perGeneMat = false(yDim,xDim,numGenes);
geneData = struct();
geneData(numGenes).geneName = '';
for i=1:numGenes
    perGeneMap = false(yDim,xDim);
    perGeneMap(wordMap==i) = true;
    geneMap(wordMap==i) = i;
    regData = regionprops(perGeneMap,sum(alignedIm,3),'Area','Centroid','PixelIdxList','PixelValues');
    validSpot = [regData.Area]>= parameters.minPixels;
   %  validSpot = cellfun(@sum,{regData.PixelValues}) > .1;

    geneData(i).geneName = libGenes(i);
    geneData(i).centroid = cat(1,regData(validSpot).Centroid);
    geneData(i).area = [regData(validSpot).Area]'; 
    geneData(i).pixelIdxList = {regData(validSpot).PixelIdxList};
    geneData(i).pixelValues = {regData(validSpot).PixelValues};

    rejectSpots = cat(1,regData(~validSpot).PixelIdxList);
    geneMap(rejectSpots) = 1+numGenes;
    perGeneMap(rejectSpots) = 1+numGenes;
    perGeneMat(:,:,i) = perGeneMap;
end


% just plotting
if parameters.saveFigures || strcmp(parameters.figVis,'on');

    % plot image of cell with spots labeled
    imagePlot = double(alignedIm - min(alignedIm(:)));
    imagePlot = uint16(imagePlot./max(imagePlot(:))*2^16);    
    genesImage = figure('Name','genesImage','visible',parameters.figVis); clf;  
    rgbIm = Ncolor(uint16(imagePlot),hsv(numBits));
    imagesc( imadjust(rgbIm, [.25 .25 .25; 1 1 1]) ); 
    set(gcf,'Units','Inches','Position',[.5 .5 10 10]);

    % plot map of cell with areas color coded by gene and labeled
    genesMap = figure('Name','genesMap','visible',parameters.figVis); clf; 
    imagesc(geneMap); colormap([jet(numGenes);zeros(1,3)]);
    colorbar; caxis([1,numGenes+1]); hold on; 
    set(gcf,'Units','Inches','Position',[.5 .5 10 10]);

    % speeding up labeling of RNA parti
    geneName = {}; geneX = []; geneY = [];
    for g=1:numGenes
        try
        geneX = [geneX; geneData(g).centroid(:,1)];
        geneY = [geneY; geneData(g).centroid(:,2)];
        geneName = [geneName; repmat(geneData(g).geneName,length(geneData(g).centroid(:,1)),1)];
        catch
        end
    end
     figure(genesImage); hold on; text(geneX,geneY,geneName,'color','w','FontSize',8); %  plot(geneX,geneY,'w.'); 
     figure(genesMap); hold on; text(geneX,geneY,geneName,'color','w','FontSize',8); %  plot(geneX,geneY,'w.'); 

    SaveFigure(genesMap,'name',['genesMap_',num2str(c,'%02d')],'formats',{'png','fig'},...
    'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
        'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);

    SaveFigure(genesImage,'name',['genesIm_',num2str(c,'%02d')],'formats',{'png','fig'},...
    'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
        'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);

    % zoomBox = [500,550,450,500];
    if isempty(parameters.zoomBox);
        parameters.zoomBox = [ceil(xDim*.4),floor(xDim*.6),ceil(yDim*.4),floor(yDim*.6)];
    end    
    zoomBox = parameters.zoomBox;
    subRegionTile = figure('Name','subRegionTile','visible',parameters.figVis); clf;
    try
        TileImage(imagePlot,'zoomBox',zoomBox,'labelLocs',[geneX,geneY],'labelNames',geneName,'gain',3);
        inbox = geneX > zoomBox(1) & geneX < zoomBox(2) & geneY > zoomBox(3) & geneY < zoomBox(4);
        subRegion = figure('Name','subRegion','visible',parameters.figVis); clf; 
        Ncolor(imagePlot(zoomBox(3):zoomBox(4),zoomBox(1):zoomBox(2),:),hsv(numBits+1)); hold on;
        text(geneX(inbox)-zoomBox(1)+1,geneY(inbox)-zoomBox(3)+1,geneName(inbox),'color','c');
        SaveFigure(subRegionTile,'name',['subRegionTile_',num2str(c,'%02d')],'formats',{'png'},...
        'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
            'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
            'savePath',parameters.savePath);     
        SaveFigure(subRegion,'name',['subRegion_',num2str(c,'%02d')],'formats',{'png'},...
        'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
            'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
            'savePath',parameters.savePath);
     catch er
          warning(er.getReport);
    end

     perGeneArea = cellfun(@sum,{geneData.area});
     perGeneCnt = cellfun(@length,{geneData.area});

    %     maxPerSpot = cell(numGenes,1);
    %     for i =1:numGenes
    %     maxPerSpot{i} = cellfun(@max,{geneData(i).pixelValues{:}});
    %     end
    %     figure(4); clf; subplot(2,1,1); hist(cat(2,maxPerSpot{idxBlanks}),20); xlim([0,400]);
    %     subplot(2,1,2); hist(cat(2,maxPerSpot{isGoodGene}),20); xlim([0,400]);

    expPerCell = figure('Name','expPerCell','visible',parameters.figVis); clf; 
    set(gcf,'Units','Inches','Position',[.5 .5 3 12]);
    subplot(1,2,1); barh(perGeneCnt); set(gca,'YTick',1:numGenes,'YTickLabel',libGenes,'FontSize',5); axis(gca,'tight'); set(gcf,'color','w'); title('spot counts','FontSize',12); 
    subplot(1,2,2); barh(perGeneArea); set(gca,'YTick',1:numGenes,'YTickLabel',libGenes,'FontSize',5);  axis(gca,'tight'); set(gcf,'color','w');  title('total area','FontSize',12);
    SaveFigure(expPerCell,'name',['expPerCell_',num2str(c,'%02d')],'formats',{'png','fig'},...
    'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
        'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);

    % sum(perGeneCnt)

    % % Correlation plots per cell
    % figure(50); clf; PlotCorr(libExpect',perGeneArea',libGenes); xlabel('FPKM'); ylabel('area');
    % figure(51); clf; PlotCorr(libExpect',perGeneCnt',libGenes); xlabel('FPKM'); ylabel('count');

    if parameters.exportSpatialImage 
        gMal = StringFind(libGenes,'MALAT1');
        gER = StringFind(libGenes,'THBS1');
        gFLNA= StringFind(libGenes,'FLNA');
        mefFig = figure('Name','mefFig','visible',parameters.figVis); clf; 
        imagesc(cat(3,perGeneMat(:,:,gFLNA),perGeneMat(:,:,gER),perGeneMat(:,:,gMal)));
        SaveFigure(mefFig,'name',['mefFig_',num2str(c,'%02d')],'formats',{'png'},...
            'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
            'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);
    end

end

