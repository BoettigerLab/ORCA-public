function stormIm = RenderSTORMimageFromTable(data,spotNum,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'scaleFactor','positive',20};
defaults(end+1,:) = {'chosenReads','integer', [1,4:6,16:19,35:38]};
defaults(end+1,:) = {'saveData','boolean', false};
defaults(end+1,:) = {'colormap','colormap','hsv'};
defaults(end+1,:) = {'showPlots','boolean', true};
defaults(end+1,:) = {'gain','positive', 1};
defaults(end+1,:) = {'verbose','boolean', true};
defaults(end+1,:) = {'numRows','integer', 3};
defaults(end+1,:) = {'overlayColor','colormap', [1 1 1]};
defaults(end+1,:) = {'saveNameTag','string',''};

parameters = ParseVariableArguments(varargin,defaults,mfilename);

%% Auto extract hybe order

numHybes = max(data.hybN) - min(data.hybN)+1; % sometimes doesn't start at 1.
numReads = max(data.readN) -min(data.readN)+1;
% determine which reads have repeats
readoutHybes = zeros(numReads,2);
for r=1:numReads  % r=4
    isRead_r = data.readN==r;
    rpts = unique(data.hybN(isRead_r));
    readoutHybes(r,1) = rpts(1);
    if length(rpts)==2
        readoutHybes(r,2) = rpts(2);
    end
end
readoutHybes = [(1:numReads)',readoutHybes];

% readout 1:numReads, repeats in order
hybeOrder = [readoutHybes(:,2); nonzeros(readoutHybes(:,3))]   ;


%% plot example spot
% spotNum = 79;
scaleFactor = parameters.scaleFactor;
chosenReads = parameters.chosenReads;
saveData = parameters.saveData;

roiX = unique(data.roi_x); % using roi_x as a unique spot identifier
currSpot = data.roi_x == roiX(spotNum);
currSpotData = data(currSpot,:);
xmin = floor(min(currSpotData.x_nm)/100);
xmax = ceil(max(currSpotData.x_nm)/100);
xBins = xmax-xmin + 1;
ymin = floor(min(currSpotData.y_nm)/100);
ymax = ceil(max(currSpotData.y_nm)/100);
yBins = ymax-ymin + 1;

inHybe = logical(prod(currSpotData.hybN~=Column(hybeOrder(chosenReads))',2)); %  find data from hybeN corresponding to readout h. 
x = currSpotData.x_nm(inHybe)/100;
y = currSpotData.y_nm(inHybe)/100;
z = currSpotData.z_nm(inHybe)/100;

if parameters.showPlots
    figure(1); clf; plot(x,y,'.','color',[.3 .3 .3],'MarkerSize',1); hold on;
end

nReads = length(chosenReads);
storm0 = RenderMList([x,y],'ROI',[ymin,ymax; xmin,xmax],'imageScale',scaleFactor);


cmap = GetColorMap(parameters.colormap,nReads);
stormIm = zeros(yBins*scaleFactor,xBins*scaleFactor,nReads);
k = 0;
com = zeros(nReads,3);
for vDim=chosenReads
    k=k+1;
    inHybe = currSpotData.hybN==hybeOrder(vDim); %  find data from hybeN corresponding to readout h. 
    x = currSpotData.x_nm(inHybe)/100;
    y = currSpotData.y_nm(inHybe)/100;
    z = currSpotData.z_nm(inHybe)/100;
    
    if parameters.showPlots
        plot(x,y,'.','color',cmap(k,:),'MarkerSize',1); hold on;
        plot(nanmedian(x),nanmedian(y),'ko','MarkerSize',20);
        plot(nanmedian(x),nanmedian(y),'+','color',cmap(k,:),'MarkerSize',20);
    end
    com(k,:) = [nanmean(x),nanmean(y),nanmean(z)];
    sm = RenderMList([x,y],'ROI',[ymin,ymax; xmin,xmax],'imageScale',scaleFactor);
    stormIm(:,:,k)= sm;
end

if parameters.showPlots
    newMap = squareform(pdist(com));
    figure(4); clf; imagesc(newMap); colorbar;  caxis([0,12]);
    GetColorMap('redToWhite');
    
    overlayFig = figure(2); clf; Ncolor(IncreaseContrast(stormIm),cmap); axis image;
    SaveFigure(overlayFig,'name',['overlayRegions_spot',num2str(spotNum),parameters.saveNameTag],'formats',{'png','fig'},'overwrite',1,'saveData',saveData);

    tileFig = figure(3); clf; 
    % tileIm = TileImage(IncreaseContrast(stormIm(:,:,:)),'showImage',false,'parameters',parameters); axis image;
    overlayIm = nansum(stormIm,3);
    TileRGBimages(IncreaseContrast(stormIm),'overlayImage',overlayIm,'colormap',cmap,'parameters',parameters); axis image;
    SaveFigure(tileFig,'name',['overlayTiled_spot',num2str(spotNum),parameters.saveNameTag],'formats',{'png','fig'},'overwrite',1,'saveData',saveData);
end



