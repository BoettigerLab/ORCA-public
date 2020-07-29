function [spots,pars] = AutoSelectSpots(imAve1,varargin)
% [spots,pars] = AutoSelectSpots(im); 
% 
defaults = cell(0,3);
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',1};
defaults(end+1,:) = {'border','nonnegative',10};
defaults(end+1,:) = {'coarseBackgroundSize','nonnegative',0}; % 12
defaults(end+1,:) = {'fineBackgroundSize','nonnegative',0};   % 50
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'showExtraPlots','boolean',false};
defaults(end+1,:) = {'numberSpots','boolean',false};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'fov','nonnegative',0};
defaults(end+1,:) = {'laplaceFilter','boolean',false};
defaults(end+1,:) = {'edgeThresh','nonnegative',0}; %.0001
defaults(end+1,:) = {'edgeSigma','nonnegative',1.2};
defaults(end+1,:) = {'overlayIm','array',[]};
defaults(end+1,:) = {'removeEdgeStack','nonnegative',20};
defaults(end+1,:) = {'figHandle','freeType',1};
defaults(end+1,:) = {'extraFigHandle','freeType',12};
% OBSOLETE, presevered for backwards compatibility
defaults(end+1,:) = {'coarseBackgroundScale','nonnegative',0}; % 50,
defaults(end+1,:) = {'fineBackgroundScale','nonnegative',0};   % 5

% parse variable arguments
pars = ParseVariableArguments(varargin, defaults, mfilename);

%% Main function
% downsample
if pars.autoSelectDownsample ~=1
    temp = imresize(imAve1,1/pars.autoSelectDownsample);  % downsample for speed 
else
    temp = imAve1;
end

% Optional background subtraction
if pars.coarseBackgroundSize > 0
    peakH = quantile(temp(:),.9999);
    imClass = class(temp);
    bkd = imresize(imresize(temp,[pars.coarseBackgroundSize,pars.coarseBackgroundSize]),size(temp));   %  figure(11); clf; imagesc(bkd); colorbar; % 
    temp(bkd<.01*peakH) = 0;
    temp = double(temp) - double(bkd);   %  figure(12); imagesc(temp); colormap(gray); colorbar; % 
    temp = temp - min(temp(:)); temp = temp/max(temp(:));
    temp = cast(temp*double(peakH),imClass);%  
    if pars.showExtraPlots
        figure(pars.extraFigHandle); clf; 
        imagesc(temp); colormap(gray); colorbar; 
        title('coarse background subtraction');
        disp('press any key to continue');
        pause;
    end
end
if pars.fineBackgroundScale > 0
    peakH = quantile(temp(:),.9999);
    imClass = class(temp);
    bkd = imresize(imresize(temp,[pars.fineBackgroundScale,pars.fineBackgroundScale]),size(temp));   %  figure(11); clf; imagesc(bkd); colorbar; % 
    temp(bkd<.01*peakH) = 0;
    temp = double(temp) - double(bkd);   %  figure(12); imagesc(temp); colormap(gray); colorbar; % 
    temp = temp - min(temp(:)); temp = temp/max(temp(:));
    temp = cast(temp*double(peakH),imClass);%  figure(12); imagesc(temp); colormap(gray); colorbar; %
    if pars.showExtraPlots
        figure(pars.extraFigHandle); clf;
        imagesc(temp); colormap(gray); colorbar; 
        title('fine background subtraction');
        disp('press any key to continue');
        pause;
    end
end

if pars.edgeThresh > 0
    bw = edge(temp,'log',pars.edgeThresh,pars.edgeSigma/pars.autoSelectDownsample);
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw,5);
    temp(~bw) = 0;
    if pars.showExtraPlots
        figure(pars.extraFigHandle); clf; 
        imagesc(temp); colormap(gray); colorbar; 
        title('edge filter');
        disp('press any key to continue');
        pause;
    end
end

% Laplacian filter to ID spots
if pars.laplaceFilter
    temp = imfilter(double(temp),fspecial('laplacian'),'symmetric');  
    temp(temp>0) = 0; 
    temp = uint16(max(temp(:)) - temp); % figure(12); clf; imagesc(temp); colorbar;
    if pars.showExtraPlots
        figure(pars.extraFigHandle); clf; 
        imagesc(temp); colormap(gray); colorbar; 
        title('Laplace filter');
        disp('press any key to continue');
        pause;
    end
end

% threshold   figure(14); clf; imagesc(temp);  % 
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
    overlayFig = figure(pars.figHandle); clf; colormap(gray);
    if isempty(pars.overlayIm)
        imagesc(imAve1);
    else
        imagesc(pars.overlayIm);
    end
    hold on; plot(spots(:,1),spots(:,2),'yo');
    if pars.numberSpots
        text(spots(:,1)+2,spots(:,2),cellstr(num2str( (1:size(spots,1))')),'color','w');
    end
    if pars.saveData
        try
            SaveFigure(overlayFig,...
            'name',['fov',num2str(pars.fov,'%03d'),'_SpotMapFig'],...
            'formats',{'png','fig'},'overwrite',true,'saveData',pars.saveData);
        catch er
            disp(er.message);
        end
    end
end