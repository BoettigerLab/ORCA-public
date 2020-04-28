function [xshift,yshift,sc,parameters] = CorrAlignScale(im1,im2,varargin)
%  [xshift,yshift,parameters] = CorrAlign(Im1,Im2)
% Inputs
% {'showplot', 'boolean', false}; show image of before and after  
% {'upsample', 'positive', 1}; upsample to get subpixel alignment
% 
% Compute xshift and yshift to align two images based on maximizing
% cross-correlation.  
% 
% removed 'region' and 'subregion' parameters 
% CorrAignFast now handles acceleration by downsampling much better

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'showplot', 'boolean', false};
defaults(end+1,:) = {'upsample', 'positive', 1};
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'minGrad', 'float', -inf};
defaults(end+1,:) = {'scales','float',1}; % .94:.01:1.02 
defaults(end+1,:) = {'verbose','boolean',false};

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


Hc = min([H,parameters.maxShift*parameters.upsample]);    
Wc = min([W,parameters.maxShift*parameters.upsample]); 
Hc2 = round(Hc/2);
Wc2 = round(Wc/2); 


numScales = length(parameters.scales);
corrPeak = zeros(numScales,1);
peakXY = zeros(numScales,2);
mapPlot = cell(numScales,2);
for s=1:numScales
    if parameters.verbose && numScales > 1
        disp(['scanning angles... fraction complete: ' num2str(s/numScales)]);
    end
    if parameters.scales==1
        im2s = im2;
    else
        sc = parameters.scales(s); % 0.97;
        T = [sc      0      0;
             0       sc     0;
             0      0       1];
         tscale = affine2d(T);
         im2s = imwarp(im2,tscale,'OutputView',imref2d(size(im2)));
    end
    corrM = xcorr2(single(im1),single(im2s)); % The correlation map
    % Just the center of the correlation map  
    corrMmini = corrM(H-Hc2+1:H+Hc2,W-Wc2+1:W+Wc2);
    [corrPeak(s),indmax] =  max(corrMmini(:));
    % use gradient maximum in place of correlation maximum
    if parameters.gradMax
%         % first derivative
%         map = gradient(corrMmini);  % figure(13); clf; imagesc(map); figure(14); clf; imagesc(corrMmini);
%         [peakPlus,indmax] =  max(map(:)); 
%         [cy1,cx1] = ind2sub(size(map),indmax );
%         [peakMinus,indmax] =  min(map(:)); 
%         [cy2,cx2] = ind2sub(size(map),indmax );
%         cy = mean([cy1,cy2]); 
%         cx = mean([cx1,cx2]);
%         corrPeak(r) = peakPlus - peakMinus;
%         peakXY(r,:) = [cx,cy];
      % second derivative
        map = gradient(gradient(corrMmini));
        [corrMin,indmax] =  min(map(:));
        corrPeak(s) = -corrMin;
        [cy,cx] = ind2sub(size(corrMmini),indmax);
        peakXY(s,:) = [cx,cy];
    else
        [cy,cx] = ind2sub(size(corrMmini),indmax);
        peakXY(s,:) = [cx,cy];
    end
    if parameters.showplot
       mapPlot{s,1} =  corrMmini;
       mapPlot{s,2} =  map;
    end
end


[parameters.corrPeak,a] = max(corrPeak);
sc = parameters.scales(a);
xshift = (peakXY(a,1)-Wc2);
yshift = (peakXY(a,2)-Hc2);

if parameters.showplot  
   subplot(2,2,1); Ncolor(cat(3,im1,im2)); 
   T = [sc  0   0;
         0  sc  0;
         0  0   1];
   im2s = imwarp(im2,affine2d(T),'OutputView',imref2d(size(im2)));
   im2s = TranslateImage(im2s,xshift,yshift);
   subplot(2,2,2); Ncolor(cat(3,im1,im2s)); 
   title(['xshift=',num2str(xshift),' yshift=',num2str(yshift), ' scale=',num2str(sc)]);
   freezeColors;
   subplot(2,2,3); imagesc(mapPlot{a,1}); colormap(jet(256));
   title(['peak=',num2str(parameters.corrPeak)]);
   subplot(2,2,4); imagesc(mapPlot{a,2}); colormap(jet(256));
   title(['peak=',num2str(parameters.corrPeak)]);
end
xshift = xshift/parameters.upsample;
yshift = yshift/parameters.upsample;

