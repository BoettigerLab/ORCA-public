function imOut = NColor(imIn,varargin)
% Takes an HxWxC image, returns an HxWx3 image which when rendered as an
% rgb image has C distinct colors.  

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'colormap', 'colormap', 'hsv'};
defaults(end+1,:) = {'contrastRGB', 'boolean', false};
defaults(end+1,:) = {'contrastHigh', 'fraction', .9999};
defaults(end+1,:) = {'contrastLow', 'fraction', 0};
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a 2D image matrix is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------

if length(varargin) == 1 % allow shorthand pass
    % pars.colormap = varargin{1};
    varin = {'colormap',varargin{1}};
else 
    varin = varargin;
end

pars = ParseVariableArguments(varin, defaults, mfilename);


% -------------------------------------------------------------------------
%% Main Function
% -------------------------------------------------------------------------
[vDim,hDim,numColors] = size(imIn);
clrmap = GetColorMap(pars.colormap,numColors); 

if size(clrmap,1) < numColors
    error(['colormap must be a colormap name or a colormap matrix of length at least ',num2str(numColors)]);
end

imOut = zeros(vDim,hDim,3,class(imIn));
for c=1:numColors
    for cc = 1:3
        imOut(:,:,cc) = imOut(:,:,cc) + imIn(:,:,c)*clrmap(c,cc);
    end
end

if pars.contrastRGB
    if ~isa(imOut,'uint16') && ~isa(imOut,'uint8')
        imOut = makeuint(imOut,16);
    end
    imOut = imadjust(imOut,stretchlim(imOut,[pars.contrastLow,pars.contrastHigh]));
end

if nargout == 0
    if ~isa(imOut,'uint16') && ~isa(imOut,'uint8')
        imOut = makeuint(imOut,16);
    end
    imagesc(imOut);
end

