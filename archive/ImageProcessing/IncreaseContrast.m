function imOut = IncreaseContrast(imIn,varargin)
% default stretch contast 0,1 
% equivelant to:  imadjust(imIn,stretchlim(imIn,[0,1]))

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'low', 'nonnegative', 0};
defaults(end+1,:) = {'high', 'nonnegative', 1};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a 2D image matrix is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
%% Main Function
% -------------------------------------------------------------------------
[~,~,nChns] = size(imIn);

high = repmat(parameters.high,nChns,1);
low = repmat(parameters.low,nChns,1);  
imOut = imIn;
for i=1:nChns
    imOut(:,:,i) = imadjust(imIn(:,:,i),stretchlim(imIn(:,:,i),[low(i),high(i)])); %   
end
