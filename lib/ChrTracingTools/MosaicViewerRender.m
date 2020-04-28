function [mosaicGray,boxCoords,xyShifts,rezero,moreOutput] = MosaicViewerRender(imageTiles,uls,varargin)
% INPUTS
% imageTiles - a Nx1 cell array of image tiles
% uls - a vector of upper-left coordinates for each of the images
%
% OUTPUTS
% mosaicGray - assembled mosaic image
% boxCoords - the corner points of each mosaic in the new image
%    this is useful if using corrAlign to tweak box positions
% 
% "Interactive"
%   Interactive mode, with or without corrAlign, asks user to approve the
%   current alignment and allows options to adjust position with mouse
%   keyboard command. May also allow user to skip placing a tile now and
%   return to place it later when more context may be available. 
%   

% set global variables
global MVR scratchPath

% ------------ Parse default Parameters
defaults = cell(0,3);
% MosaicViewerRender pars
defaults(end+1,:) = {'interactive','boolean',false};
defaults(end+1,:) = {'corrAlign','boolean',false};
defaults(end+1,:) = {'saveFolder','string',scratchPath}; % xyShifts determined using interactive mode will be saved here 
defaults(end+1,:) = {'saveName','string','regFOV'}; % xyShifts determined selected using interactive mode will be saved with this name
defaults(end+1,:) = {'saveShifts','boolean',false}; % xyShifts saved anyway whether or not using interactive 
defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'};
defaults(end+1,:) = {'showProgress','boolean',false};
defaults(end+1,:) = {'progressRes','fraction',1/5};
defaults(end+1,:) = {'downsample','integer',1};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'reZeroUL','boolean',false};
defaults(end+1,:) = {'minOverlap','fraction',.05};
defaults(end+1,:) = {'corrFig','integer',9}; % figure to use for correlation plot 
defaults(end+1,:) = {'interactFig','integer',30}; % figure to use for interactive plot 
defaults(end+1,:) = {'xyShifts','freeType',[]}; %
defaults(end+1,:) = {'padMosaic','nonnegative',3}; % in multiples of tile size
defaults(end+1,:) = {'stopOnError','boolean',false}; % in multiples of tile size
% flatten background
defaults(end+1,:) = {'flatten','boolean',false};
defaults(end+1,:) = {'background','freeType',[]};
% mosaic tile adjustments
defaults(end+1,:) = {'transpose','boolean',false};
defaults(end+1,:) = {'fliplr','boolean',false};
defaults(end+1,:) = {'flipud','boolean',false};
% parameters for prepping CorrAlign
defaults(end+1,:) = {'padAlign','integer',250};
defaults(end+1,:) = {'corrAlignHigh','boolean',.9999};
defaults(end+1,:) = {'corrAlignLow','boolean',.8};
defaults(end+1,:) = {'sortData','boolean',true};
% fov sort for CorrAlign
defaults(end+1,:) = {'fovSortThreshold','fraction',.7};
defaults(end+1,:) = {'fovSortResize','fraction',.02};
defaults(end+1,:) = {'fovSortDebug','boolean',false};
% parameters forwarded to CorrAlign
defaults(end+1,:) = {'showplot', 'boolean', false};
defaults(end+1,:) = {'saveCorr', 'boolean', false};
defaults(end+1,:) = {'maxSize', 'positive', 300};
defaults(end+1,:) = {'maxShift', 'nonnegative', 0}; % 0= auto to 1/5th image height 
defaults(end+1,:) = {'minGrad', 'float', 50};
defaults(end+1,:) = {'savePath', 'string', 'G:\Alistair\CorrAlign\'};
defaults(end+1,:) = {'showExtraPlot','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%------------------------Main Function-------------------------------%

% if pars.interactive
%     figure(pars.interactFig); close; pause(.01);
%     figH = figure(pars.interactFig);   pause(.01);
%     set(figH,'KeyPressFcn',{@FigKeyPress});
%     MVR.im3 = zeros(10); 
%     MVR.figH = figH;
%     ResetShifts();
% end

if istable(pars.xyShifts)
   pars.xyShifts = pars.xyShifts{:,:}; % convert table to matrix 
end
% % pad matrix to allow for shifts
% if ~isempty(pars.xyShifts)
%    maxXshift = max(pars.xyShifts(:,1)); % still doesn't handle negative shifts
%    maxYshift = max(pars.xyShifts(:,2));
%    minXshift = min(pars.xyShifts(:,1)); % still doesn't handle negative shifts
%    minYshift = min(pars.xyShifts(:,2));
% else
%     maxXshift = 0;
%     maxYshift = 0;
%     minXshift = 0;
%     minYshift = 0;
% end

% -------- Downsample for fast display if requested
nTiles = length(imageTiles);
if pars.downsample~=1
    if pars.verbose
        disp('compressing images...');
    end
    pars.xyShifts = round(pars.xyShifts/pars.downsample);
    uls = ceil(uls/pars.downsample); %  ceil avoids 0 index errors 
    for m=1:nTiles
        imageTiles{m} = imresize(imageTiles{m},1/pars.downsample);
    end
else
    uls = ceil(uls); % still enforse integer values, ceil avoids 0 index errors 
end

%---- Flip tiles if requested
for m=1:nTiles
    if pars.transpose % in same order as DaxToImageTiles
        imageTiles{m} = imageTiles{m}';
    end
    if pars.fliplr
        imageTiles{m} = fliplr(imageTiles{m});
    end
    if pars.flipud
       imageTiles{m} = flipud(imageTiles{m}); 
    end
end

% --------- Flatten images based on average if requested 
bkd = [];
minTiles = 10;
if pars.flatten && nTiles >= minTiles
    if pars.verbose
        disp('flattening images...');
    end
    if isempty(pars.background)
        bkd = nanmedian(double(cat(3,imageTiles{:})),3);
    else
        bkd = pars.background;
    end
    pk = max(bkd(:));
    if pars.showplot
        figure(9); clf; imagesc(bkd); colormap(gray(2E3)); colorbar;
        title('illumination profile');
    end
    cls = class(imageTiles{1});
    for m=1:nTiles
        imageTiles{m} = cast(pk*double(imageTiles{m})./double(bkd),cls);
    end
elseif pars.flatten && nTiles < minTiles
    if pars.verbose
       disp('fewer than 5 tiles, skipping image based flattening.'); 
    end
end


% ------- compute mosaic size 
if ~pars.corrAlign && ~pars.interactive
    pars.padAlign = 0;
end
% uls(:,1) = uls(:,1) - minXshift;
% uls(:,2) = uls(:,2) - minYshift;
xmax = max(uls(:,1));
ymax = max(uls(:,2));
[h_i,w_i] = size(imageTiles{1});
if pars.reZeroUL
    xmin = min(uls(:,1));
    ymin = min(uls(:,2));
    h_m = round(ymax-ymin+1+h_i+pars.padMosaic*h_i+pars.padAlign ); % + maxYshift - minYshift
    w_m = round(xmax-xmin+1+w_i+pars.padMosaic*w_i+pars.padAlign); %  + maxXshift - minXshift
    rezero = -round([-xmin+pars.padMosaic*w_i+1,-ymin+pars.padMosaic*h_i+1]);
    rezero = max([0,0; rezero]);
    uls = [uls(:,1)-rezero(1),uls(:,2)-rezero(2)];
else
    h_m = round(ymax+1+h_i+pars.padMosaic*h_i+pars.padAlign ); % + maxYshift - minYshift
    w_m = round(xmax+1+w_i+pars.padMosaic*w_i+pars.padAlign ); % + maxXshift - minXshift
    rezero = [0,0];
end
mosaicGray = zeros(h_m,w_m);
boxCoords = zeros(nTiles,4);


% -------- Allows precomputed xyShifts to be used. 
if isempty(pars.xyShifts)
    xyShifts = zeros(nTiles,2);
    xyShiftsOut = nan(nTiles,2); % to distinguish no data from no-move
    hasData = false(nTiles,1);
else
    if size(pars.xyShifts,1)==nTiles
        xyShifts = pars.xyShifts;
        xyShiftsOut = xyShifts;
        hasData = true(nTiles,1);
    else
       warning(['MosaicViewerRender recieved ',num2str(nTiles),' tiles and ',num2str(size(pars.xyShifts,1)),' xyShifts']);
       error('any passed xyShifts must by the same size as nTiles'); 
    end
end

% ------- Check if a previous FOV file already exists
if pars.corrAlign
    hasData = false(nTiles,1); 
   if exist([pars.saveFolder,pars.saveName,'.csv'],'file')
        answer = questdlg(['found existing registration map. Load it or overwrite?'], ... % question
                        'Use existing', ...  % pop-up label
                        'Use existing','Overwrite existing','Cancel',... % op1 op2 op3 
                        'Use existing');  % default op
       if strcmp(answer,'Use existing')
           [xyShifts, xyShiftsOut, hasData] = LoadRegFov([pars.saveFolder,pars.saveName,'.csv']);
           % tableData = readtable(['C:\Data\Scratch\shifts.csv']);
       elseif strcmp(answer,'Cancel')
           error('Operation canceled by user.');
       end
   end
end

% xyShifts doesn't allow NaNs. 
%   Make these into 0s and track this by updating hasData
%   This is just a safety check in case this hasn't been cleaned up.
noDataRows = any(isnan(xyShifts),2);
xyShifts(noDataRows,:) = 0;
hasData(noDataRows) = false;        

% ------- Sort data if requested and using corrAlign
sortIndex = 1:nTiles;
if pars.corrAlign
    [h_i,w_i] = size(imageTiles{1});
    if pars.sortData
        if pars.verbose
            disp('sorting FOVs for optimal tiling...');
        end
        [sortIndex,~] = SortFOVsByOverlap([h_m,w_m], [h_i,w_i],...
            uls,'data',imageTiles,...
            'resize',pars.fovSortResize,'dataThreshold',pars.fovSortResize,'debug',pars.fovSortDebug);  
        [~,sortBack] = sort(sortIndex);
        uls = uls(sortIndex,:);
        imageTiles = imageTiles(sortIndex);
        % if we preloaded these, we have to sort them too
        xyShifts = xyShifts(sortIndex,:);
        xyShiftsOut = xyShiftsOut(sortIndex,:);
        hasData = hasData(sortIndex); 
    end
    if pars.maxShift == 0 
       pars.maxShift = h_i/5;
    end
    if pars.saveCorr
        SetFigureSavePath(pars.savePath);
    end
end

% -------- precompute contrast if using show display
if pars.showProgress
   im = imresize(imageTiles{1},pars.progressRes);
   cmax = quantile(im(:),.999);
end



% -------- The main assembly loop, add images into the mosaic
%  Optionally, adjusts tile positions based on 
%      * CorrAlign,
%      * Previously computed adjustments "XYshifts"
%      * Interactive input
% 
if pars.verbose
    disp('assembling tiles...');
end

k=0;
origIndex = 1:nTiles;
origIndex = origIndex(sortIndex);
toPlace = 1:nTiles;
orderPlaced = zeros(1,nTiles);
flagged = false(1,nTiles);
nPlaced = sum(toPlace==0);
while nPlaced < nTiles    
    n = find(toPlace~=0,1,'first');
    m = toPlace(n);
    try
        addToMosaic = true;
        currTile = double(imageTiles{m});
        [h_i,w_i] = size(currTile);
        boxCoords(m,:) = [uls(m,1),uls(m,1)+w_i-1,uls(m,2),uls(m,2)+h_i-1];
        if hasData(m)
           if pars.veryverbose
              disp(['using existing shifts for position ', num2str(origIndex(m))]);
           end
        end
        if ~hasData(m) && (pars.corrAlign || pars.interactive)
            % add rotation and scaling (tform) support later [I think these are added now, just need testing?] 
            [alignValues, ~, valid] = InteractiveAlign(mosaicGray,currTile,...
                'startPosition',boxCoords(m,[1,3]),...
                'refImageCnst',true,...
                'outputTform',false,...
                'combine','off',...
                'imName',[num2str(origIndex(m)),'  tile:',num2str(m),' of ',num2str(nTiles)],...
                'parameters',pars);
            if valid == 1 %% keep %  ~MVR.shifts.skip && ~MVR.shifts.cancel
                xyShifts(m,:) = round([xyShifts(m,1)+alignValues.xshift,xyShifts(m,2)+alignValues.yshift]); % currently only placed to pixel accuracy
                addToMosaic = true;
            elseif valid == 2  %% MVR.shifts.skip
                % move on and come back to this FOV later
                addToMosaic = false;
                toPlace = [toPlace,m]; %#ok<AGROW>
                toPlace(n) = [];   
                if pars.verbose
                   disp(['skipping panel ',num2str(origIndex(m)),' will try again after other tiles are placed.']); 
                end
            elseif valid == 3
                xyShifts(m,:) = round([xyShifts(m,1)+alignValues.xshift,xyShifts(m,2)+alignValues.yshift]);
                addToMosaic = true;
                flagged(m) = true;
                if pars.verbose
                   disp(['flagged panel ',num2str(origIndex(m))]); 
                end
            end
        end
        if addToMosaic
            % add to mosaic
            boxCoords(m,:) = round([boxCoords(m,1:2)+xyShifts(m,1),boxCoords(m,3:4)+xyShifts(m,2)]);
            mosaicGray(boxCoords(m,3):boxCoords(m,4),boxCoords(m,1):boxCoords(m,2)) = ...
                CombineImages(mosaicGray(boxCoords(m,3):boxCoords(m,4),boxCoords(m,1):boxCoords(m,2)),...
                currTile,'parameters',pars);
            toPlace(n) = 0;
            k=k+1;
            orderPlaced(k) = m;
            % save shift data if it was interactively mapped
            %   this is a bit challenging, since we place the images in a
            %   computationally determined order, not necessarily the order
            %   in which they are passed to the function
            % we write them in the order they were passed. 
            % we analyze them in the order the script recommends.
            if pars.corrAlign && (pars.interactive || pars.saveShifts) && ~hasData(m)
                    xyShiftsOut(m,:) = xyShifts(m,:); % 
                    if pars.sortData
                        xyShiftsSave = xyShiftsOut(sortBack,:);
                    else
                        xyShiftsSave = xyShiftsOut;
                    end
                    tableStruct.xShift = xyShiftsSave(:,1);
                    tableStruct.yShift = xyShiftsSave(:,2);
                    tableStruct.tileOrder = sortBack';
                    tableData = struct2table(tableStruct);
                    tableName = [pars.saveFolder,pars.saveName,'.csv'];
                    writetable(tableData,tableName); % will overwrite
            end      
        end
        % show progress
        if pars.showProgress             
            x1 = max(1,boxCoords(m,1)-3*w_i);
            x2 = min(w_m,boxCoords(m,2)+3*w_i);
            y1 = max(1,boxCoords(m,3)-3*h_i);
            y2 = min(h_m,boxCoords(m,4)+3*h_i);
            im = imresize(mosaicGray(y1:y2,x1:x2),pars.progressRes);
            figure(1); clf; imagesc(im); 
            caxis([0,cmax]);
            colormap(gray);
        end
    catch er
        cprintf([1,0,0],['error at m=',num2str(m), ' tile ',num2str(origIndex(m)) '.']);
        if pars.stopOnError
            disp('place debug to stop here');
            error(er.getReport);
        else
             cprintf([1,0,0],'Atttempting to skip and continue...');
             toPlace(n) = 0;
             k=k+1;
             orderPlaced(k) = m;
        end
    end
    nPlaced = sum(toPlace==0);
end

if pars.verbose
   disp('finished placing all tiles'); 
end

% convert back to original format (much better compression)
if isa(imageTiles{1},'uint16')
   mosaicGray = uint16(mosaicGray); 
elseif isa(imageTiles{1},'uint8')
    mosaicGray = uint8(mosaicGray); 
elseif isa(imageTiles{1},'single')
    mosaicGray = single(mosaicGray); 
end

% beyond 4 outputs, aggegate in a labeled structure
moreOutput.flagged = flagged;
moreOutput.bkd = bkd;

% resort data to original order
if pars.corrAlign && pars.sortData
    boxCoords = boxCoords(sortBack,:);
    xyShifts = xyShifts(sortBack,:);
end


