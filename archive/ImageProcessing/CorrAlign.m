function [xshift,yshift,parameters] = CorrAlign(im1,im2,varargin)
%  [xshift,yshift,parameters] = CorrAlign(im1,im2)
% Inputs
% {'region', 'nonnegative', 200}; max number of pixels to use 
% {'showplot', 'boolean', false}; show image of before and after  
% {'upsample', 'positive', 1}; upsample to get subpixel alignment
% 
% Compute xshift and yshift to align two images based on maximizing
% cross-correlation.  
% 
% -------------------------------------------------------------------------
% Update to work with NaNs


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'region', 'nonnegative', inf};
defaults(end+1,:) = {'showplot', 'boolean', false};
defaults(end+1,:) = {'upsample', 'positive', 1};
defaults(end+1,:) = {'subregion', 'boolean', true};
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};

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


if parameters.subregion && parameters.region < H
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

% figure(11); clf; 
% subplot(1,2,1); imagesc(im1); 
% subplot(1,2,2); imagesc(im2);



corrM = xcorr2(im1,im2); % The correlation map
Hc = min([H,parameters.region*parameters.upsample,parameters.maxShift*2]);    
Wc = min([W,parameters.region*parameters.upsample,parameters.maxShift*2]); 
Hc2 = round(Hc/2);
Wc2 = round(Wc/2); 
% Just the center of the correlation map  
corrMmini = corrM(H-Hc2+1:H+Hc2,W-Wc2+1:W+Wc2);
[peakCorr,indmaxCorr] =  max(corrMmini(:)); 
if parameters.gradMax
    map = gradient(corrMmini);  % figure(13); clf; imagesc(map); figure(14); clf; imagesc(corrMmini);
    [peakPlus,indmax] =  max(map(:)); 
    [cy1,cx1] = ind2sub(size(map),indmax );
    [peakMinus,indmax] =  min(map(:)); 
    [cy2,cx2] = ind2sub(size(map),indmax );
    cy = mean([cy1,cy2]); 
    cx = mean([cx1,cx2]);
else
    [cy,cx] = ind2sub(size(corrMmini),indmaxCorr );
end

xshift = round(cx-Wc2);
yshift = round(cy-Hc2);

if parameters.showplot
   dispIm1 = makeuint(im1,8);
   dispIm2 = makeuint(im2,8);
   subplot(1,3,1); Ncolor(cat(3,dispIm1,dispIm2));
   dispIm2 = TranslateImage(dispIm2,xshift,yshift);
   subplot(1,3,2); Ncolor(cat(3,dispIm1,dispIm2)); 
   title(['xshift=',num2str(xshift),' yshift=',num2str(yshift)]);
   freezeColors;
   if parameters.gradMax
       subplot(1,6,5); imagesc(map); colormap(jet(256));
       title(['peak=',num2str(peakPlus-peakMinus)]);
       subplot(1,6,6); imagesc(corrMmini); colormap(jet(256));
       title(['peak=',num2str(peakCorr)]);
       parameters.corrPeak = peakPlus-peakMinus;
   else
       subplot(1,3,3); imagesc(corrMmini); colormap(jet(256));
    title(['peak=',num2str(peakCorr)]);
    parameters.corrPeak = peakCorr;
   end
   
end

xshift = xshift/parameters.upsample;
yshift = yshift/parameters.upsample;

