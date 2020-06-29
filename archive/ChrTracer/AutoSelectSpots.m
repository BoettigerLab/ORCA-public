function [spots,pars] = AutoSelectSpots(imAve1,varargin)
% [spots,pars] = AutoSelectSpots(im); 
% 
defaults = cell(0,3);
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'border','nonnegative',10};
defaults(end+1,:) = {'coarseBackgroundScale','nonnegative',0}; % 50
defaults(end+1,:) = {'fineBackgroundScale','nonnegative',0};   % 5
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'gain','positive',1};
defaults(end+1,:) = {'fov','nonnegative',1};
defaults(end+1,:) = {'laplaceFilter','boolean',true};
defaults(end+1,:) = {'edgeThresh','nonnegative',0}; %.0001
defaults(end+1,:) = {'edgeSigma','nonnegative',1.2};
defaults(end+1,:) = {'overlayIm','array',[]};
defaults(end+1,:) = {'removeEdgeStack','nonnegative',20};

% parse variable arguments
pars = ParseVariableArguments(varargin, defaults, mfilename);

%% Main function
% downsample
temp = imresize(imAve1,1/pars.autoSelectDownsample);  % downsample for speed 

% Optional background subtraction
if pars.coarseBackgroundScale > 0
    bkd = imresize(imresize(temp,1/pars.coarseBackgroundScale),size(temp));  %  figure(11); clf; imagesc(bkd);  % 
    temp = 4*temp ./ bkd;  %  figure(12); imagesc(temp); colormap(gray); % 
end
if pars.fineBackgroundScale > 0
    bkd2 = imresize(imresize(temp,1/pars.fineBackgroundScale),size(temp));  % figure(11); clf; imagesc(bkd2);  % 
    temp = 4*temp ./ bkd2; % figure(12); clf; imagesc(temp); 
end

if pars.edgeThresh > 0
    bw = edge(temp,'log',pars.edgeThresh,pars.edgeSigma/pars.autoSelectDownsample);
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw,5);
    temp(~bw) = 0;
end

% Laplacian filter to ID spots
if pars.laplaceFilter
    temp = imfilter(double(temp),fspecial('laplacian'),'symmetric');  
    temp(temp>0) = 0; 
    temp = uint16(max(temp(:)) - temp); % figure(12); clf; imagesc(temp); colorbar;
end

% threshold
mask = IncreaseContrast(temp,'low',pars.autoSelectThreshold,'high',1); %   figure(12); clf; imagesc(mask);  % 
bw = imregionalmax(mask);                                              %  figure(13); imagesc(bw);
[y,x] = ind2sub(size(mask),find(bw>0));
[h,w] = size(mask);
onEdge = x > w-pars.border | x<pars.border | y > h-pars.border | y<pars.border;

% remove points which stack up vertically or horizontally
if pars.removeEdgeStack > 0
    [v,n] = occurrences( round(x));
    bad = v(n>10);
    for b=1:length(bad)
        onEdge = onEdge | round(x)==bad(b);
    end
    [v,n] = occurrences( round(y));
    bad = v(n>10);
    for b=1:length(bad)
        onEdge = onEdge | round(y)==bad(b);
    end
end

x(onEdge) = []; 
y(onEdge) = [];
spots = round(pars.autoSelectDownsample*[x,y] - pars.autoSelectDownsample/2);


if pars.showPlots
    overlayFig = figure(1); clf; 
    if isempty(pars.overlayIm)
        imagesc(pars.gain*imAve1);
    else
        imagesc(pars.gain*pars.overlayIm);
    end
    hold on; plot(spots(:,1),spots(:,2),'yo');
    text(spots(:,1)+2,spots(:,2),cellstr(num2str( (1:size(spots,1))')),'color','w'  );
    if pars.saveData
        SaveFigure(overlayFig,...
            'name',['fov',num2str(pars.fov,'%03d'),'_SpotMapFig'],...
            'formats',{'fig'},...
            'overwrite',true,'saveData',pars.saveData);
    end
end