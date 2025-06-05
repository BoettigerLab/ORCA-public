function [xy,xz,yz] = ProjectIm4D(im4D,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'caxis','array',[]};
defaults(end+1,:) = {'colormap', 'colormap', 'hsv'};
defaults(end+1,:) = {'contrastRGB', 'boolean', false};
defaults(end+1,:) = {'contrastHigh', 'fraction', .9999};
defaults(end+1,:) = {'contrastLow', 'fraction', 0};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[h,w,z,nChannels] = size(im4D);
xy = zeros(h,w,nChannels);
yz = zeros(z,h,nChannels);
xz = zeros(z,w,nChannels);
for c=1:nChannels
    im3D = im4D(:,:,:,c);
    xy(:,:,c) = max(im3D,[],3);
    yz(:,:,c) = max( permute(im3D,[3,1,2]), [],3);
    xz(:,:,c) = max( permute(im3D,[3,2,1]), [],3);
end
if pars.showPlots || nargout==0
   subplot(1,3,1); Ncolor(xy,'parameters',pars);  xlabel('x'); ylabel('y');
   subplot(1,3,2); Ncolor(yz,'parameters',pars);  xlabel('y'); ylabel('z');
   subplot(1,3,3); Ncolor(xz,'parameters',pars);   xlabel('x'); ylabel('z');
end