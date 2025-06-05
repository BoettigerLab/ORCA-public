function imOut = ChromCorrectImage(imIn,tform3D,varargin)
% make sure the units are correct for the chromatic warp
% currently just a 2D warp using the 3D data
%    should update this for 3D stacks. 

defaults = cell(0,3);
% data loading
defaults(end+1,:) = {'showPlot','boolean',false}; 
defaults(end+1,:) = {'res','integer',2}; % 
defaults(end+1,:) = {'nppXY','positive',108}; % 
pars = ParseVariableArguments(varargin,defaults,mfilename);


res=pars.res;
npp = pars.nppXY; 


[imH,imW,imZ] = size(imIn); % note, dim3 currently ignored 
[Y,X]=meshgrid(1/res:1/res:imH,1/res:1/res:imW); % convert pixels to points
ys = Y(:)*npp;
xs = X(:)*npp; 
nPts = length(xs);

% also get the numerical values of these pixles;  
idM = imresize(reshape((1:(imH*imW)),[imH,imW]),res,'nearest');  % first we need their indices 
id = idM(:);
vs = imIn(id);  % values at each x,y position
% [nPts, length(vs)]  % sanity check


xyc = tforminv(tform3D,ys,xs,zeros(nPts,1)); % apply tform   
xyp = round(xyc/npp*res);

im2c = zeros(res*imH,res*imW,'uint16');
idp = sub2ind(res*[imH,imW],xyp(:,2),xyp(:,1));  % 
im2c(idp) = vs;
imOut = imresize(im2c,1/res); % here we use the default interpolation to downsample and redistribute the sub-pixel values

% plotting for troubleshooting; 
if nargout == 0 || pars.showPlot
    ys = 1:500; xs=1:500;
    % figure(3); clf; 
    % subplot(1,2,1); imagesc(IncreaseContrast( imIn(ys,xs), 'high',.9999,'low',.3) );
    % subplot(1,2,2); imagesc(IncreaseContrast( imOut(ys,xs), 'high',.9999,'low',.3) );
    imC = IncreaseContrast( cat(3, imIn(ys,xs), imOut(ys,xs)), 'high',.9999,'low',.3) ;
    figure(4); clf; Ncolor(imC);
end