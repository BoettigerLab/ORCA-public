function imOut = NColor(imIn,varargin)
% Takes an HxWxC image, returns an HxWx3 image which when rendered as an
% rgb image has C distinct colors.  

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'colormap', 'colormap', 'hsv'};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a 2D image matrix is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
if length(varargin) == 1
    parameters.colormap = varargin{1};
else
    parameters = ParseVariableArguments(varargin, defaults, mfilename);
end

% -------------------------------------------------------------------------
%% Main Function
% -------------------------------------------------------------------------
[vDim,hDim,numColors] = size(imIn);
clrmap = GetColorMap(parameters.colormap,numColors); 

if size(clrmap,1) < numColors
    error(['colormap must be a colormap name or a colormap matrix of length at least ',num2str(numColors)]);
end

imOut = zeros(vDim,hDim,3,class(imIn));
for c=1:numColors
    for cc = 1:3
        imOut(:,:,cc) = imOut(:,:,cc) + imIn(:,:,c)*clrmap(c,cc);
    end
end

if nargout == 0
    if ~isa(imOut,'uint16') & ~isa(imOut,'uint8')
        imOut = makeuint(imOut,16);
    end
    imagesc(imOut);
end

