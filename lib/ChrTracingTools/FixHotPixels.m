function imOut = FixHotPixels(imIn,varargin)
% imOut = FixHotPixels(imIn,varargin)
% Correct hot pixels using calibration data
% defaults(end+1,:) = {'test','boolean',false};
% defaults(end+1,:) = {'threshold','fraction',.99};
% defaults(end+1,:) = {'scope',{'scope1','scope2'},'scope2'};
% defaults(end+1,:) = {'refDax','string','G:\Alistair\2018-07-09_Scope2calibrations\hot_Pixels_0001.dax'}; % at least 10 frames, match size

global HotPixel

defaults = cell(0,3);
defaults(end+1,:) = {'test','boolean',false};
defaults(end+1,:) = {'threshold','fraction',.99};
defaults(end+1,:) = {'scope',{'scope1','scope2'},'scope2'};
defaults(end+1,:) = {'refDax','string','G:\Alistair\2018-07-09_Scope2calibrations\hot_Pixels_0001.dax'}; % at least 10 frames, match size
pars = ParseVariableArguments(varargin,defaults,mfilename); 

if isempty(HotPixel) % for speed
    hotPixels = ReadDax(pars.refDax);
    hot1 = nanmean(hotPixels(:,:,1:5),3);
    HotPixel.hot1 = hot1;
    HotPixel.scope = pars.scope;
else
    hot1 = HotPixel.hot1;
end
renorm = double(quantile(imIn(:),pars.threshold));
zs = size(imIn,3);
imOut = imIn;
for z=1:zs
    imOut(:,:,z) = uint16( double(imIn(:,:,z))./double(hot1)*renorm );
end

if pars.test
    hot2 = nanmean(hotPixels(:,:,6:end),3);
    hotNorm = uint16(double(hot1)./double(hot2)*quantile(hot2(:),pars.threshold));
    figure(1); clf; imagesc(hot1);  colormap(gray); colorbar; caxis([0,300]);
    figure(2); clf; imagesc(hotNorm); colormap(gray); colorbar; caxis([0,300]);
    disp(['original pixel variance: ',std(hot1(:))./mean(hot1(:))]);
    disp(['normalized pixel variance: ',std(hotNorm(:))./mean(hotNorm(:))]);
end