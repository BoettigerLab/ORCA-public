function [im3] = BlendOverlapImages(im1,im2,xshift,yshift,varargin)
% Ramp the overlapping area from 1 to 0 and then add the images
% see slide 52 of https://web.stanford.edu/class/cs231m/lectures/lecture-5-stitching-blending.pdf
% for both x + y shifts this needs some thought


defaults = cell(0,3);
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'displayHigh','nonnegative',.999};
defaults(end+1,:) = {'displayLow','nonnegative',.1};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% keep zeros
im1(isnan(im1)) = 0;
im2(isnan(im2)) = 0;

% extract some image parameters
[h1,w1] = size(im1);
[h2,w2] = size(im2);
rs = 1:h1;
cs = 1:w1;
if h2 < h1
    im2 = padarray(im2,[h1-h2+1,0],'post');
end
if w2 < w1
    im2 = padarray(im2,[0,w1-w2+1],'post');
end
im2 = im2(rs,cs); % images must be same size; 

% blend images
blend = ones(length(rs),length(cs));
blend1 = ones(length(rs),length(cs));
if xshift >= 0
  blendEdge = repmat(linspace(1,0,xshift),length(rs),1);
  blend(:,end-xshift+1:end) =  blend(:,end-xshift+1:end).*blendEdge;
  blendEdge = repmat(linspace(0,1,xshift),length(rs),1); % 
  blend1(:,end-xshift+1:end) =  blend1(:,end-xshift+1:end).*blendEdge;
elseif xshift < 0 
  blendEdge = repmat(linspace(0,1,-xshift),length(rs),1);
  blend(:,1:-xshift) =  blend(:,1:-xshift).*blendEdge;
  blendEdge = repmat(linspace(1,0,-xshift),length(rs),1); % 
  blend1(:,1:-xshift) =  blend1(:,1:-xshift).*blendEdge;
end
if  yshift >= 0
  blendEdge = repmat(linspace(1,0,yshift)',1,length(cs));
  blend(end-yshift+1:end,:) =  blend(end-yshift+1:end,:).*blendEdge;
  blendEdge = repmat(linspace(0,1,yshift)',1,length(cs)); % 
  blend1(end-yshift+1:end,:) =  blend1(end-yshift+1:end,:).*blendEdge;
elseif yshift < 0
  blendEdge = repmat(linspace(0,1,-yshift)',1,length(cs));
  blend(1:-yshift,:) =  blend(1:-yshift,:).*blendEdge; 
  blendEdge = repmat(linspace(1,0,-yshift)',1,length(cs));
  blend1(1:-yshift,:) =  blend1(1:-yshift,:).*blendEdge; 
end



% for places where im2 has no data, don't decrease values of im1
% blend(im2==0) = 1; blend(isnan(im2)) = 1;


blend2 = 1-blend;
% blend(im1==0) = 1;

im3 = nansum( cat(3, im1.*blend2, im2.*blend ), 3);

if pars.showPlots
 
   im1d = IncreaseContrast(uint16(im1),'high',pars.displayHigh,'low',pars.displayLow);
   im2d = IncreaseContrast(uint16(im2),'high',pars.displayHigh,'low',pars.displayLow);
   im3d = IncreaseContrast(uint16(im3),'high',pars.displayHigh,'low',pars.displayLow);
    
    subplot(2,3,1); imagesc(im3d); title('Blended');
    subplot(2,3,2); imagesc(nansum(cat(3,im1d, im2d),3)); colormap(gray); title('Added');
    subplot(2,3,3); Ncolor(cat(3,im1d,im2d)); title('overlay');
    subplot(2,3,4); imagesc(im1d); title('im1');
    subplot(2,3,5); imagesc(im2d); title('im2');
    
    figure(14); clf;
    subplot(1,3,1); imagesc(blend2); colorbar;
    subplot(1,3,2); imagesc(blend); colorbar;
    subplot(1,3,3); imagesc(blend+blend2); colorbar;
    colormap(gray);
end
