function [mosaic,uls] = TilesToMosaic(imTiles,uls,varargin)
% simple function
% takes tiles and positions and creates a blended mosaic
% no shifting, no automatic alignment etc, 
% simplifed from MosaicViewerRender
% imTiles is an Nx1 stack of tiles
% uls is a Nx2 list of upper left positions
defaults = cell(0,3);
defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'};
defaults(end+1,:) = {'padMosaic','nonnegative',1}; % pad on both sides in multiples of tile size
defaults(end+1,:) = {'downsample','positive',1}; % downsample this fold. 1=no downsampling
defaults(end+1,:) = {'mosaicSize','float',[0,0]}; % if uls have aleady been aligned to some common reference frame, we don't want to disrupt that
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'rezero','boolean',true}; % not sure this is useful - may fix value and remove option in future.  
defaults(end+1,:) = {'flipV','boolean',false};
defaults(end+1,:) = {'flipH','boolean',false};
defaults(end+1,:) = {'transpose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

nTiles = length(imTiles);
uls_in = uls;

% ---- flip if requested
if pars.flipV || pars.flipH || pars.transpose
    for m=1:nTiles
        if pars.flipV
            imTiles{m} = flipud(imTiles{m});
        end
        if pars.flipH
            imTiles{m} = fliplr(imTiles{m});
        end
        if pars.transpose
            imTiles{m} = imTiles{m}';
        end
    end
end

% -------- Downsample for fast display if requested

if pars.downsample~=1
    if pars.verbose
        disp('compressing images...');
    end
    uls = ceil(uls/pars.downsample); %  ceil avoids 0 index errors 
    for m=1:nTiles
        imTiles{m} = imresize(imTiles{m},1/pars.downsample);
    end
else
    uls = ceil(uls); % still enforse integer values, ceil avoids 0 index errors 
end


%--- Determint tile dimensions
[h_i,w_i] = size(imTiles{1});
nFOVs = length(imTiles);

% if a mosaic size is passed, use this.
if sum(pars.mosaicSize) ~= 0 
    pars.padMosaic = 0;
    pars.rezero = false;
    h_a = pars.mosaicSize(1);
    w_a = pars.mosaicSize(2);
end

% if no mosiac size is given (size 0) automatically compute the necessary
% size from the coordinate positions and the tile size. 
if sum(pars.mosaicSize) == 0 
    x_max = ceil(max(uls(:,1)));
    y_max = ceil(max(uls(:,2)));
    if pars.rezero
        x_min = floor(min(uls(:,1)));
        y_min = floor(min(uls(:,2)));
    else
        x_min = 0;
        y_min = 0;
    end
    h_a = y_max-y_min+1+h_i + pars.padMosaic*h_i*2; % pad on both sides
    w_a = x_max-x_min+1+w_i + pars.padMosaic*w_i*2;
    rezero = [x_min,y_min];
    uls = uls - repmat(rezero - pars.padMosaic*[h_i,w_i],nFOVs,1); % +1
end




% The real function: just drop these tiles in the UL positions given
boxCoords = zeros(nFOVs,4);
mosaic = zeros(h_a,w_a,class(imTiles{1}));
for m=1:nFOVs
    try
        if sum(uls_in(m,:)) == 0
            if pars.verbose
                disp('skipping panel with no placement data');
            end
        else
    boxCoords(m,:) = round([uls(m,1),uls(m,1)+w_i-1,uls(m,2),uls(m,2)+h_i-1]);
    mosaic(boxCoords(m,3):boxCoords(m,4),boxCoords(m,1):boxCoords(m,2)) = ...
                    CombineImages(mosaic(boxCoords(m,3):boxCoords(m,4),boxCoords(m,1):boxCoords(m,2)),...
                                  imTiles{m},'parameters',pars);
        end
    catch er
        disp(er.getReport)
        disp('place debug here');
    end
end