function spots = FindSpots(im,varargin)

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'autoSelectSigma','positive',3};
defaults(end+1,:) = {'bkdFilterScale','nonnegative',10};
defaults(end+1,:) = {'showPlot','boolean',false};
defaults(end+1,:) = {'MarkerSize','positive',10};

% parse variable arguments
pars = ParseVariableArguments(varargin, defaults, mfilename);


% Auto Select ROI  (This is quick and dirty)
%   for improvement: imgaussfilt3 pre-process, field-flattening...
if size(im,3) > 1
    im = uint16(nanmean(im,3)); % average all hybes
end
if pars.autoSelectDownsample ~=1 
    im = imresize(im,1/pars.autoSelectDownsample); 
end
% optionally, flatten field by subtracting background
% background is determined by aggressive downsampling of the image 
if pars.bkdFilterScale > 0
    bkd = imresize(imresize(im,1/pars.bkdFilterScale),size(im));
    im = im - bkd;
end
% threshold the resulting image
if pars.autoSelectThreshold ~= 1
    mask = IncreaseContrast(im,'low',pars.autoSelectThreshold,'high',1);   % figure(11); imagesc(imAve);  % 
else
    mask = im;
end
bw = imregionalmax(mask);                            %    figure(12); imagesc(bw);
[y,x] = ind2sub(size(mask),find(bw>0));
spots = round(pars.autoSelectDownsample*[x,y] - pars.autoSelectDownsample/2);

if nargout == 0 || pars.showPlot
    % imagesc(im); colormap(gray);
    hold on; 
    plot(spots(:,1),spots(:,2),'yo','MarkerSize',pars.MarkerSize);
end