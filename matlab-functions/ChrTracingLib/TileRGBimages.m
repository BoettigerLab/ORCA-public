function newImage = TileRGBimages(imageStack,varargin)
% newImage = TileRGBimages(stormIm,...
%           'overlayImage',storm0,'overlayColor',[.25,.25,.25]);

defaults = cell(0,3);
defaults(end+1,:) = {'multicolor', 'boolean', true}; 
defaults(end+1,:) = {'showImage', 'boolean', false}; 
defaults(end+1,:) = {'numRows', 'positive', 4}; 
defaults(end+1,:) = {'overlayColor','colormap',[1,1,1]};
defaults(end+1,:) = {'overlayImage','array',0};
defaults(end+1,:) = {'imSepColor', 'nonnegative', 2^16}; 
defaults(end+1,:) = {'colormap','colormap','hsv'};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
cmap = GetColorMap(parameters.colormap);

[yDim,xDim,nImages] = size(imageStack);
numRows = parameters.numRows; 
numColumns = ceil(nImages/numRows);
xyDim = [numRows, numColumns].*[yDim+1,xDim+1];
newImage = parameters.imSepColor*ones(xyDim(1),xyDim(2),3);
for i=1:nImages
    [colInd, rowInd] = ind2sub([numColumns numRows], i);
    indY = (0:(yDim-1)) + 1 + (rowInd-1)*(yDim+1); 
    indX = (0:(xDim-1)) + 1 + (colInd-1)*(xDim+1);
    newImage(indY,indX,1) = parameters.overlayColor(1,1)*parameters.overlayImage + cmap(i,1)*(imageStack(:,:,i));
    newImage(indY,indX,2) = parameters.overlayColor(1,2)*parameters.overlayImage + cmap(i,2)*(imageStack(:,:,i));
    newImage(indY,indX,3) = parameters.overlayColor(1,3)*parameters.overlayImage + cmap(i,3)*(imageStack(:,:,i));
end

if nargout == 0
    imagesc(newImage); 
end
