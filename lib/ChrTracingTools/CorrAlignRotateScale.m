function [alignValues,parameters] = CorrAlignRotateScale(im1,im2,varargin)
%  [xshift,yshift,angle,rescale,parameters] = CorrAlignRotateScale(im1,im2,varargin)
%  [tform,parameters] = CorrAlignRotateScale(im1,im2,varargin)
%  Note: tform matrix is dependent on image size for rotation, as affine
%  rotation is around 0,0 and needs to be shifted to the center of the
%  image to match imrotate. 
%          
% Inputs
% {'showplot', 'boolean', false}; show image of before and after  
% {'upsample', 'positive', 1}; upsample to get subpixel alignment
% 
% Compute xshift and yshift to align two images based on maximizing
% cross-correlation.  Allows for rescaling and rotation if requested.
% 


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'showplot', 'boolean', false};
defaults(end+1,:) = {'upsample', 'positive', 1};
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'minGrad', 'float', -inf};
defaults(end+1,:) = {'angles','float',0}; % in degrees -10:1:10
defaults(end+1,:) = {'scales','float',1}; % fractional expansion/contraction  .97:.01:1.03
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
numAngles = length(parameters.angles);
corrPeak = zeros(numAngles,numScales);
peakXY = zeros(numAngles,numScales,2);
mapPlot = cell(numAngles,numScales,2);
k = 0;
try
    for r=1:numAngles
        for s=1:numScales
            k=k+1;
            if parameters.verbose && (numScales > 1 || numAngles >1)
                disp(['scanning... fraction complete: ' num2str(k/(numScales*numAngles))]);
            end
            if parameters.scales==1
                im2s = im2;
            else
                sc = parameters.scales(s); % 0.97;
                T = [sc      0      0;
                     0       sc     0;
                     0      0       1];
                 im2s = imwarp(im2,affine2d(T),'OutputView',imref2d(size(im2))); % first rescale
            end
            im2sr = imrotate(im2s,parameters.angles(r),'bilinear','crop');  % then rotate
            corrM = xcorr2(single(im1),single(im2sr)); % The correlation map  % then translate
            % Just the center of the correlation map  
            corrMmini = corrM(H-Hc2+1:H+Hc2,W-Wc2+1:W+Wc2);
            [corrPeak(r,s),indmax] =  max(corrMmini(:));
            % use gradient maximum in place of correlation maximum
            if parameters.gradMax
                % second derivative
%                 [dX,dY] = gradient(corrMmini(:,:));
%                 ddX = gradient(dX);
%                 [~,ddY] = gradient(dY);
%                 map = ddX+ddY;
                % figure(29); clf; imagesc((ddX+ddY));
                 map = gradient(gradient(corrMmini));
                
                [corrMin,indmax] =  min(map(:));
                corrPeak(r,s) = -corrMin;
                [cy,cx] = ind2sub(size(corrMmini),indmax);
                peakXY(r,s,:) = [cx,cy];
            else
                [cy,cx] = ind2sub(size(corrMmini),indmax);
                peakXY(r,s,:) = [cx,cy];
                map = zeros(size(corrMmini));
            end
            if parameters.showplot
               mapPlot{r,s,1} =  corrMmini;
               mapPlot{r,s,2} =  map;
            end
        end
    end
catch er
    disp(er.getReport);
    error('something wrong here');
end


[~,i] = max(corrPeak(:));
[a,z] = ind2sub([numAngles,numScales],i);
a=a; z=z;
parameters.corrPeak = corrPeak(a,z);
sc = parameters.scales(z);
theta = parameters.angles(a);
xshift = (peakXY(a,z,1)-Wc2);
yshift = (peakXY(a,z,2)-Hc2);


if parameters.showplot  
   subplot(2,2,1); Ncolor(cat(3,im1,im2));
   T = [sc  0    0;
         0  sc   0;
         0  0   1];
   im2s = imwarp(im2,affine2d(T),'OutputView',imref2d(size(im2)));
   im2p = imrotate(im2s,theta,'bilinear','crop');
   im2p = TranslateImage(im2p,xshift,yshift);
   subplot(2,2,2); Ncolor(cat(3,im1,im2p)); 
   title(['xshift=',num2str(xshift),' yshift=',num2str(yshift), ' scale=',num2str(sc), ' angle=',num2str(theta)]);
   freezeColors;
   subplot(2,2,3); imagesc(mapPlot{a,z,1}); colormap(jet(256));
   title(['peak=',num2str(parameters.corrPeak)]);
   subplot(2,2,4); imagesc(mapPlot{a,z,2}); colormap(jet(256));
   title(['peak=',num2str(parameters.corrPeak)]);
end
xshift = xshift/parameters.upsample;
yshift = yshift/parameters.upsample;

alignValues.xshift = xshift;
alignValues.yshift = yshift;
alignValues.theta = theta;
alignValues.rescale = sc;

