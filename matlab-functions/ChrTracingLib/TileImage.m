function [newImage,labelOffsets,flatImage] = TileImage(imageIn,varargin)
% TileImage
% Tiles the different images in the 3D stack imageIn into a single image
% 
% defaults(end+1,:) = {'multicolor', 'boolean', true}; 
% defaults(end+1,:) = {'showImage', 'boolean', true}; 
% defaults(end+1,:) = {'showLabels', 'boolean', true}; 
% defaults(end+1,:) = {'numRows', 'positive', 4}; 
% defaults(end+1,:) = {'zoomBox', 'array', []}; 
% defaults(end+1,:) = {'labelLocs', 'positive', []}; 
% defaults(end+1,:) = {'labelNames', 'cell', {}}; 
% defaults(end+1,:) = {'tileLabels', 'cell', {}}; 
% defaults(end+1,:) = {'gain', 'positive', 1}; 
% defaults(end+1,:) = {'imSepColor', 'nonnegative', 2^16}; 
% defaults(end+1,:) = {'colormap','colormap','hsv'};

defaults = cell(0,3);
defaults(end+1,:) = {'multicolor', 'boolean', true}; 
defaults(end+1,:) = {'showImage', 'boolean', false}; 
defaults(end+1,:) = {'showLabels', 'boolean', true}; 
defaults(end+1,:) = {'numRows', 'positive', 4}; 
defaults(end+1,:) = {'zoomBox', 'array', []}; 
defaults(end+1,:) = {'labelLocs', 'positive', []}; 
defaults(end+1,:) = {'labelNames', 'cell', {}}; 
defaults(end+1,:) = {'tileLabels', 'cell', {}}; 
defaults(end+1,:) = {'gain', 'positive', 1}; 
defaults(end+1,:) = {'imSepColor', 'nonnegative', 2^16}; 
defaults(end+1,:) = {'colormap','colormap','hsv'};
defaults(end+1,:) = {'FontSize','positive',7};
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults);

% figure(101); Ncolor(imageIn,hsv(16+1));
%% 


dType = class(imageIn);

[yDim,xDim,numColors] = size(imageIn);
numRows = parameters.numRows; 
numColumns = ceil(numColors/numRows);

% determine if a zoomBox was requested
if isempty(parameters.zoomBox)
    zoomBox = [1,xDim,1,yDim];
else
    zoomBox = parameters.zoomBox;
    imageIn = imageIn(zoomBox(3):zoomBox(4),zoomBox(1):zoomBox(2),:);
    [yDim,xDim,numColors] = size(imageIn);
end



if parameters.multicolor
    numChns = numColors;
else
    numChns = 1;
end

xyDim = [numRows, numColumns].*[yDim+1,xDim+1];
if isa(imageIn,'uint16')
    newImage = parameters.imSepColor*ones(xyDim(1),xyDim(2),numChns,'uint16');
else
    newImage = parameters.imSepColor*ones(xyDim(1),xyDim(2),numChns,dType);    
end

% get label locations and names and trim to fit zoomBox if necessary
if parameters.showLabels && ~isempty(parameters.labelLocs)
    labelLocs = parameters.labelLocs;
    labelNames = parameters.labelNames;
    inbox = labelLocs(:,1) > zoomBox(1) & labelLocs(:,1) < zoomBox(2) & labelLocs(:,2) > zoomBox(3) & labelLocs(:,2) < zoomBox(4);
    labelLocs = labelLocs(inbox,:); 
    labelNames = labelNames(inbox);
    
    numLabels = sum(inbox); 
    labelPrintLocs = zeros(numLabels*numColors,2);
    labelPrintNames = cell(numLabels*numColors,1);
end

labelOffsets = [];
labelPos = [];
gap = 0;
for h=1:numColors % j = 5
    [colInd, rowInd] = ind2sub([numColumns numRows], h);
    indY = (0:(yDim-1)) + 1 + (rowInd-1)*(yDim+1); 
    indX = (0:(xDim-1)) + 1 + (colInd-1)*(xDim+1);
    hybeJ  = imageIn(:,:,h);
    newImage(indY,indX,:) = zeros(length(indY),length(indX),numChns,dType);
    if parameters.multicolor
        newImage(indY,indX,h) = hybeJ;
    else
        newImage(indY,indX,1) = hybeJ;
    end
    labelPos(h,:) = [indX(1) + xDim/10, indY(1) + yDim/10]; %#ok<AGROW,*SAGROW>
    labelOffsets(h,:) = [indX(1)-1, indY(1)-1]; %#ok<AGROW,*SAGROW>
    if parameters.showLabels && ~isempty(parameters.labelLocs)
        labelPrintLocs((h-1)*numLabels+1:h*numLabels,:)  = [labelLocs(:,1)-zoomBox(1)+1 + indX(1),labelLocs(:,2)-zoomBox(3)+1 + indY(1)];
        labelPrintNames((h-1)*numLabels+1:h*numLabels) = labelNames;
    end
end

if ischar(parameters.colormap)
    cmap = GetColorMap(parameters.colormap,numColors);
    flatImage = Ncolor(parameters.gain*newImage,cmap); 
else
    flatImage = max(newImage,[],3);
end

if parameters.showImage || nargout == 0
    imagesc(flatImage); hold on;
    if parameters.showLabels  && ~isempty(parameters.labelLocs)
        text(labelPrintLocs(:,1),labelPrintLocs(:,2),labelPrintNames,'color','w','FontSize',5);
        if isempty(parameters.tileLabels)
            text(labelPos(:,1),labelPos(:,2), cellstr(num2str( (1:numColors)')),'color','w','FontSize',parameters.FontSize);
        end
    end
    if ~isempty(parameters.tileLabels) && parameters.showLabels
        text(labelPos(:,1)-2,labelPos(:,2), parameters.tileLabels,'color','w','interpreter','none' ,'FontSize',parameters.FontSize );
    end
end