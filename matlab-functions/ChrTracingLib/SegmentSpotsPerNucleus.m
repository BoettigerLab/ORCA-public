function [fidTable,cellID_full,figHandle] = SegmentSpotsPerNucleus(fidIm,varargin)
% take a fiducial 2D image 
% segment cells using cellpose
% returns brightest spots per cell
% 
% Notes
% fidImage is used both to find the spots per cell and to find the nuclear
%       boundaries.
% By default rerunCellpose is false, so that any existing cellpose
%   segmentation maps will be loaded. This accelerates execution. This
%   function can be called in series to 

defaults = cell(0,3);
defaults(end+1,:) = {'nucImage','freeType',[]}; % optional, give an image of just the nucleus without the DNA spots to use for nuclear segmentation. 
defaults(end+1,:) = {'cellID','freeType',[]}; % optional, give a precomputed cell segmentation map for this image
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'f','integer',0};
defaults(end+1,:) = {'saveFolder','string',''};
% pars passed to AutoFitSpots
defaults(end+1,:) = {'border','nonnegative',10}; % exclude spots this distance in pixels from the border 
defaults(end+1,:) = {'autoSelectDownsample','nonnegative',3};  %autoSelectDownsample
defaults(end+1,:) = {'diameter','nonnegative',0};
defaults(end+1,:) = {'imSize','array',[512,512]};
defaults(end+1,:) = {'nucleiThreshold','fraction',.985};
defaults(end+1,:) = {'displayThreshold','fraction',.9999};
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'rerunCellpose','boolean',false};
defaults(end+1,:) = {'figShowResult','freeType',10};
defaults(end+1,:) = {'spotsPerCell','integer',2};
defaults(end+1,:) = {'figName','string','fidSpots'};
defaults(end+1,:) = {'tableRoot','string','fidTable'};
defaults(end+1,:) = {'figFormats','cell',{'png','fig'}};
defaults(end+1,:) = {'env','string','mlab_cellpose'}; % cellpose env
pars = ParseVariableArguments(varargin,defaults,mfilename);

figHandle = [];
analysisFolder = pars.saveFolder;
f = pars.f;
if ~isempty(analysisFolder)
    figName = ['fov',num2str(f,'%03d'),'_',pars.figName];
    figFile = [analysisFolder,figName];
    cellIDsFile = [analysisFolder,'fov',num2str(f,'%03d'),'_cellIDs.png'];
    fidTableFile =[analysisFolder,'fov',num2str(f,'%03d'),'_',pars.tableRoot,'.csv'];
else
    fidTableFile = '';
end
run = true;

% if data exists, we just load that and move on.
if exist(fidTableFile,'file')  && ~pars.overwrite
    run=false;
    fidTable = readtable(fidTableFile);
    if nargout > 1
        cellID_full = imread(cellIDsFile);
    end
    % if asked for a figure file, we will need to load one
    if nargout > 2
        if exist([figFile,'.fig'],'file')
            uiopen([figFile,'.fig'],1);
            pause(.01); 
            figHandle = gca;
        elseif exist([figFile,'.png'],'file') && pars.figShowResult
            oldFig = imread([figFile,'.png']);
            figHandle = figure(pars.figShowResult); 
            imagesc(oldFig);
        else
            run = true;
            warning('reruning')
        end
    end
    % otherwise we run analysis
end
if run
    
    if isempty(pars.nucImage)
        nucImage = fidIm;
    else
        nucImage = pars.nucImage;
    end
    % downsample and increase contrast
    nucSmall = imresize(nucImage,pars.imSize);
    nucSmall = makeuint(nucSmall,8);  
    nucContrast = IncreaseContrast(nucSmall,'high',pars.nucleiThreshold);
    % fidSmall = imresize(fidIm,pars.imSize);

    % ---- Cellpose based segmentation
    if isempty(pars.cellID)
        % RunCellpose creates a new save folder and cellpose will by analyze
        % all images in that folder regardless of their name.   If passed a
        % cell array of images, all the images will be converted to compressed
        % pngs and then analyzed.
        allMasks = FindFiles([analysisFolder,'*_cp_masks.tif']);
        if length(allMasks)>1 && pars.overwrite 
            cprintf([1 0 0],['error in saveFolder: ',analysisFolder]);
            error('found multiple existing masks in this folder. SegmentSpotsPerNucleus may confuse idenitity');
        end
    
        % Optionally, use a separate image for the nuclear segmentation from
        % the spot selection image
        cellID = RunCellpose(nucContrast,'model','nuclei','diameter',pars.diameter,'imSize',size(nucContrast),...
            'figShowLoadedImages',0,'saveFolder',analysisFolder,'overwrite',pars.rerunCellpose,'env',pars.env); 
    else
        cellID = imresize(pars.cellID,size(fidIm));
    end
    % % just for troubleshooting -- show results of cellpose
%     mask = boundarymask(cellID); 
%     im = labeloverlay(IncreaseContrast(nucSmall,'high',pars.displayThreshold),mask,'Transparency',0);
%     figure(4); clf; imagesc(im); colormap(gray); colorbar;  hold on;

    %%
    % fit lots of points, score by brightness (we'll take the brightest N per cell, generally brightests 2) 
    xy = AutoSelectSpots(fidIm,'autoSelectThreshold',.75,'laplaceFilter',false,...
        'autoSelectDownsample',pars.autoSelectDownsample,'border',pars.border,...
        'showPlots',false,'removeEdgeStack',0);  
    % Assign spots to cells at full size
    spot_xyLin = sub2ind(size(fidIm),xy(:,2),xy(:,1));
    spot_brightness = fidIm(spot_xyLin);

    % just for troubleshooting
%  figure(1); clf; 
%  cellID_full = imresize(cellID,size(fidIm),'nearest');
%  mask = boundarymask(cellID_full); 
%  im = labeloverlay(IncreaseContrast(fidIm,'high',.9999),mask,'Transparency',0);
%  sc = size(fidIm)/size(fidSmall,1);
%  y = sc(1)*xy(:,2) - sc(1)/2;
%  x = sc(2)*xy(:,1) - sc(2)/2;
%  imagesc(im); hold on; 
%  plot(x,y,'yo'); colormap(gray);
%  id = spot_brightness > quantile(spot_brightness,.95)
%  plot(sc(2)*xybc(id,1),sc(1)*xybc(id,2),'ro'); colormap(gray);

    cellProps = regionprops(cellID,'PixelIdxList');
    nSpots = size(xy,1);
    nCells = length(cellProps);
    spotsPerCell = pars.spotsPerCell;
    xybc = cat(2,xy,spot_brightness,zeros(nSpots,2));
    for c=1:nCells
        cell_xyLin = cellProps(c).PixelIdxList; % linear index of cell xy position
        [pixel_idx,spot_id]  = intersect(spot_xyLin,cell_xyLin); % in cell
        [spotBright,bright_id] = sort(spot_brightness(spot_id),'descend'); % brightness of all spots in cell;
        if length(bright_id)>spotsPerCell+1
            next_bright = spotBright(spotsPerCell+1);
        else
            next_bright = 0;
        end
        if length(bright_id)>spotsPerCell
            bright_id = bright_id(1:spotsPerCell);
        end
        xybc(spot_id(bright_id),5) = c;
        xybc(spot_id(bright_id),4) = next_bright;
    end

    [~,sortIdx] = sort(xybc(:,4));
    spotMatrix = xybc(sortIdx,:);
    % plot(spotMatrix(:,1),spotMatrix(:,2),'yo');
    spotMatrix(spotMatrix(:,4)==0,:) =  [];
    % plot(spotMatrix(:,1),spotMatrix(:,2),'ro');

    % full size
    cellID_full = imresize(cellID,size(fidIm),'nearest');
    mask = boundarymask(cellID_full); 
    im = labeloverlay(IncreaseContrast(fidIm,'high',.9999),mask,'Transparency',0);
    % sc = size(fidIm)/size(fidSmall,1);
    % y = sc(1)*spotMatrix(:,2) - sc(1)/2;
    % x = sc(2)*spotMatrix(:,1) - sc(2)/2;
    % spotMatrix(:,2) = y; % update with new values
    % spotMatrix(:,1) = x;
    % add more data columns and convert to table
    nTraces = size(spotMatrix,1);
    spotMatrix = [spotMatrix, (1:nTraces)',f*ones(nTraces,1)];
    fidTable = array2table(double(spotMatrix));
    fidTable.Properties.VariableNames = {'x_pix','y_pix','brightness','contrast','cellID','traceID','fov'};
    
    % plot segmentation results and spots
    if pars.figShowResult
        figHandle = figure(pars.figShowResult); 
        clf; imagesc(im); 
        colormap(gray); colorbar;  
        hold on;
        plot(spotMatrix(:,1),spotMatrix(:,2),'ro');
    end
    % %     add columns 
    %% save results
  
    
    if pars.overwrite
        imwrite(cellID_full,cellIDsFile);
        writetable(fidTable,fidTableFile);
        if  pars.figShowResult
            SetFigureSavePath(analysisFolder);
            SaveFigure(figHandle,'name',figName,'formats',pars.figFormats,'overwrite',pars.overwrite);
        end
    end
    
end 
