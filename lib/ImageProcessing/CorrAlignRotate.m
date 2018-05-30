function [xshift,yshift,angle,parameters] = CorrAlignRotate(im1,im2,varargin)
%  [xshift,yshift,parameters] = CorrAlign(Im1,Im2)
% Inputs
% {'region', 'nonnegative', 200}; max number of pixels to use 
% {'showplot', 'boolean', false}; show image of before and after  
% {'upsample', 'positive', 1}; upsample to get subpixel alignment
% 
% Compute xshift and yshift to align two images based on maximizing
% cross-correlation.  

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'region', 'nonnegative', inf};
defaults(end+1,:) = {'showplot', 'boolean', false};
defaults(end+1,:) = {'upsample', 'positive', 1};
defaults(end+1,:) = {'subregion', 'boolean', true};
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'angles','float',-10:1:10};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'two 2D image matrices are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

[H,W] = size(im1);


if parameters.subregion % && parameters.region < H
    h1 = max(1,round(H/2)-parameters.region/2+1);
    h2 = min(H,round(H/2)+parameters.region/2);
    hs = h1:h2;
    w1 = max(1,round(W/2)-parameters.region/2+1);
    w2 = min(W,round(W/2)+parameters.region/2);
    ws = w1:w2;
    im1 = im1(hs,ws);
    im2 = im2(hs,ws);
end
if parameters.upsample ~= 1
    im1 = imresize(im1,parameters.upsample);
    im2 = imresize(im2,parameters.upsample);
end


[H,W] = size(im1);
% normalize data and deal with nans
im1 = single(im1);
im2 = single(im2);
im1 = (im1 - nanmean(im1(:))) / nanstd(im1(:));
im1(isnan(im1(:))) = 0;
im2 = (im2 - nanmean(im2(:))) / nanstd(im2(:));
im2(isnan(im2(:))) = 0;


Hc = min([H,parameters.region*parameters.upsample,parameters.maxShift]);    
Wc = min([W,parameters.region*parameters.upsample,parameters.maxShift]); 
Hc2 = round(Hc/2);
Wc2 = round(Wc/2); 

angles = parameters.angles; % -10:1:10;
numAngles = length(angles);
corrPeak = zeros(numAngles,1);
indmax = zeros(numAngles,1);
for r=1:numAngles
    im2rot = imrotate(im2,angles(r),'bilinear','crop');
    corrM = xcorr2(single(im1),single(im2rot)); % The correlation map
    % Just the center of the correlation map  
    corrMmini = corrM(H-Hc2+1:H+Hc2,W-Wc2+1:W+Wc2);
    [corrPeak(r),indmax(r)] =  max(corrMmini(:));
end
[parameters.corrPeak,a] = max(corrPeak);
angle = angles(a);
[cy,cx] = ind2sub([Hc,Wc],indmax(a) );
xshift = (cx-Wc2);
yshift = (cy-Hc2);

if parameters.showplot
    
   subplot(1,3,1); Ncolor(cat(3,im1,im2));
   im2p = TranslateImage(im2,xshift,yshift);
   im2p = imrotate(im2p,angle,'bilinear','crop');
   subplot(1,3,2); Ncolor(cat(3,im1,im2p)); 
   title(['xshift=',num2str(xshift),' yshift=',num2str(yshift), ' rotCW=',num2str(angle)]);
   freezeColors;
   subplot(1,3,3); imagesc(corrMmini); colormap(jet(256));
   title(['peak=',num2str(parameters.corrPeak/1E9)]);
end

xshift = xshift/parameters.upsample;
yshift = yshift/parameters.upsample;

