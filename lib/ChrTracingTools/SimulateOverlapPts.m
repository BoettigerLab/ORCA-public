function [imStack,xy] = SimulateOverlapPts(varargin)
% simulate points
% imStack = SimulateOverlapPts() returns 
% 
% 
% 

defaults = cell(0,3);
%           {parameter name, parameter type, default value}
defaults(end+1,:) = {'xy','freeType',[]};
defaults(end+1,:) = {'hs','nonnegative',[]};
defaults(end+1,:) = {'sigma','nonnegative',1.2};
defaults(end+1,:) = {'outputSize','array',[30,30]};
defaults(end+1,:) = {'maxSep','nonnegative',5};
defaults(end+1,:) = {'nPts','nonnegative',4};  % number of emitters to simulate   
defaults(end+1,:) = {'precision','nonnegative',.10}; % 
defaults(end+1,:) = {'scale','nonnegative',10}; % 
defaults(end+1,:) = {'nIms','nonnegative',1}; % number of images to simulate
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'pixelNoise','integer',0}; % pixel noise  
pars = ParseVariableArguments(varargin,defaults,mfilename);

% parse parameters
h = pars.outputSize(1);
w = pars.outputSize(2);
p = pars.maxSep;
nPts = pars.nPts;
sc = 1/pars.precision;
nIms = pars.nIms;
im0 = zeros(h,w,'uint16');
psfL = makeuint(fspecial('gaussian',9*sc,pars.sigma*sc),16)/uint16(2^7); % want this to be smooth
% psf = imresize(psfL,1/sc);
% figure(1); clf; imagesc(psfL); colorbar;
% figure(2); clf; imagesc(psf);  colorbar; colormap(gray);
imL0 = imresize(im0,sc);
[hp, wp] = size(psfL);
% random molecule positions
if isempty(pars.xy)
    xy = h/2-p/2 + p*rand(nPts,2,nIms);
    xy = round(sc*xy)/sc; % make it an even interval for exactness
else
    xy = pars.xy;
    nPts = size(xy,1);
end
% random heights
if isempty(pars.hs)
    hs = 1+1*rand(nIms,nPts);
end
% assemble stacks
imStack = zeros(h,w,nIms,'uint16');
for n=1:nIms  
    xyp = xy(:,:,n)-.5;
    xyL = round(xyp*sc);
    imL = imL0;
    for i=1:nPts
        ymin = xyL(i,2)-(hp/2)+1;
        ymax = xyL(i,2)+(hp/2);
        xmin = xyL(i,1)-(wp/2)+1;
        xmax = xyL(i,1)+(wp/2);
        imL(ymin:ymax,xmin:xmax) = imL(ymin:ymax,xmin:xmax) + hs(n,i)*psfL;
    end
    im = imresize(imL,1/sc); % figure(2); clf; imagesc(imL);
    im = im + uint16(pars.pixelNoise*rand(h,w));
    imStack(:,:,n) = im;        
%     figure(1); clf; imagesc(im); hold on;
%     plot(xy(:,1,n),xy(:,2,n),'m*'); pause(.1);
end
% show results
if pars.showPlots
    imagesc(im); hold on;
    plot(xy(:,1,n),xy(:,2,n),'m*'); 
    colorbar;
end
