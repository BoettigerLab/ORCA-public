function subHandles = TileImageStack(imageIn,varargin)
%% Inputs
% imageIn - a HxWxC image which will be plotted as a tile of C separate HxW
% images. The number of rows and columns in the tile can be selected
% (default is 4 rows), columns computed automatically.
% 
% defaults(end+1,:) = {'numRows', 'positive', 4}; 
% defaults(end+1,:) = {'gap', 'nonnegative', 1}; 
% defaults(end+1,:) = {'colorTiles', 'boolean', true}; 
% defaults(end+1,:) = {'colormap', 'colormap', 'hsv'}; 
% defaults(end+1,:) = {'showPanelNumber', 'boolean', true}; 
% defaults(end+1,:) = {'fontSize', 'positive', 10}; 
% defaults(end+1,:) = {'tileLabels','freeType',{}};
%
% Alistair Boettiger
% CC BY  - May 6, 2017
% 


%%
defaults = cell(0,3);
defaults(end+1,:) = {'numRows', 'positive', 4}; 
defaults(end+1,:) = {'gap', 'nonnegative', 1}; 
defaults(end+1,:) = {'colorTiles', 'boolean', true}; 
defaults(end+1,:) = {'colormap', 'colormap', 'hsv'}; 
defaults(end+1,:) = {'autoContrast', 'boolean', true}; 
defaults(end+1,:) = {'showPanelNumber', 'boolean', true}; 
defaults(end+1,:) = {'fontSize', 'positive', 10}; 
defaults(end+1,:) = {'tileLabels','freeType',{}};
defaults(end+1,:) = {'imWhite','freeType',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);
%%

fh = gcf;
[yDim,xDim,numTiles] = size(imageIn);
numRows = pars.numRows; 
numColumns = ceil(numTiles/numRows);

panelWidth = fh.Position(3)/numColumns;
panelHeight = fh.Position(4)/numRows;

if pars.colorTiles
    cmap = GetColorMap(pars.colormap,numTiles);
    cmap = cast(cmap,class(imageIn));
end

if isempty(pars.imWhite)
    imWhite = zeros(yDim,xDim,class(imageIn));
else
    imWhite = pars.imWhite;
end

% Plot the data in separate subpanels, create handles for each subpanel.
subHandles = cell(numTiles,1);
for h=1:numTiles % j = 5 
    figure(fh);
    subHandles{h} = subplot(numRows,numColumns,h,'Units','pixels');   %  slow
    if pars.colorTiles 
        if pars.autoContrast
            imTemp = IncreaseContrast(imageIn(:,:,h));
        else
            imTemp = imageIn(:,:,h);
        end
        im = zeros(yDim,xDim,3,class(imageIn));
        im(:,:,1) = cmap(h,1)*imTemp + imWhite;
        im(:,:,2) = cmap(h,2)*imTemp + imWhite;
        im(:,:,3) = cmap(h,3)*imTemp + imWhite;
    else
        im = imageIn(:,:,h);
    end
    imagesc(im); hold on;
    try
    if pars.showPanelNumber && isempty(pars.tileLabels)
        text(.1*xDim,.1*yDim,num2str(h),'FontSize',pars.fontSize,'color','w');
    elseif ~isempty(pars.tileLabels)
        text(.1*xDim,.1*yDim,pars.tileLabels{h},'FontSize',pars.fontSize,'color','w');
    end
    catch
    end
    axis off;
end
   
% now we adjust the position of the panels. 
% (This needs to be in a separate loop, otherwise matlab randomly deletes
% some of the panels)
for h=1:numTiles
    [colInd, rowInd] = ind2sub([numColumns numRows], h);
     % subHandles{h}.Units = 'pixels';
     subHandles{h}.Position = [(colInd-1)*panelWidth+pars.gap,...
                               fh.Position(4)-(rowInd)*panelHeight+pars.gap,...
                               panelWidth-pars.gap,panelHeight-pars.gap];
     subHandles{h}.Units = 'normalized';    
end

