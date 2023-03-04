function [spotTablePerRead,cellProj,cropSpotsPerRead] = SpotsInCellsFromDax(cellBorders,fov,eTableXLS,varargin)
% Inputs
%
%  Outputs
%     spotTablePerRead  - table of fitted spot locations in FOV
%     cellProj          - cell array, nCells x nHybes x nRperHyb x 2
%                         filtered xy and xz projections if each nuc
%    cropSpotsPerRead   - zoom in 3D crop of the brightest 2 spots from
%                         each cell

%
% maybe this should be a class / app where the child function parameters
% could be handled differently.

defaults = cell(0,3);
% self parameters
defaults(end+1,:) = {'dataFolder','string',''};
defaults(end+1,:) = {'analysisFolder','string',''};
defaults(end+1,:) = {'hybFolderRoot','string','Hyb'};
defaults(end+1,:) = {'daxRoot','string','ConvZscan_'}; % all dax files start with this root
defaults(end+1,:) = {'padFrames','integer',10}; % frames to remove from start and end of dax file 
defaults(end+1,:) = {'skipFirst','integer',0}; % skip first frame 0=no, 1=yes
defaults(end+1,:) = {'nChns','integer',3}; % number of channels in dax. alternating channel assumed
defaults(end+1,:) = {'idxFid','integer',3}; % index of fiducial channel, if used. 0 if no fiducial data  
defaults(end+1,:) = {'showRawImageFig','integer',1}; % handle for a figure in which to plot the raw image. set to 0 to supress figure display.
defaults(end+1,:) = {'xyTileFig','integer',1}; % handle for a figure in which to plot the raw image. set to 0 to supress figure display.
defaults(end+1,:) = {'xzTileFig','integer',2}; % handle for a figure in which to plot the raw image. set to 0 to supress figure display.
defaults(end+1,:) = {'vecImageFig','integer',3}; % handle for a figure in which to plot the raw image. set to 0 to supress figure display.
defaults(end+1,:) = {'xyNucWindow','positive',1}; % cropping area = radius x this for xy and xz Tile Figs
defaults(end+1,:) = {'radi','positive',50}; % cropping area = radi x xyNucWindow for xy and xz Tile Figs
defaults(end+1,:) = {'maxHybs','nonnegative',0}; % 0=auto detect. Or optionally set the max number of hybs (good for troubleshooting to just test a few). 
% parameters passed to FindSpots
defaults(end+1,:) = {'autoSelectThreshold','fraction',.99}; 
defaults(end+1,:) = {'bkdFilterScale','positive',5}; 
% parameters passed to SpotsInNucleus
defaults(end+1,:) = {'maxSpotsPerCell','integer',2}; 
defaults(end+1,:) = {'brightWindow','integer',1}; % expand by this sq-radius the area used for brightness peak
defaults(end+1,:) = {'xyCropWindow','integer',8}; % expand by spot-fitting area by this sq-radius the area cropped for PSF fitting
defaults(end+1,:) = {'zCropWindow','integer',14}; % expand by spot-fitting z-stack by this sq-radius the area cropped for PSF fitting
defaults(end+1,:) = {'troubleshoot','boolean',false};
defaults(end+1,:) = {'refineFit','boolean',true}; % highly recommend - Gaussian fitting for sub-pixel accuracy and "spot-like" data 
    % fit PSF pars
    defaults(end+1,:) = {'keepBrightest', 'integer', 1};
    defaults(end+1,:) = {'maxFitWidth', 'positive', 8}; 
    defaults(end+1,:) = {'maxFitZdepth', 'positive', 14};
    defaults(end+1,:) = {'initSigmaXY','positive',1.25};
    defaults(end+1,:) = {'initSigmaZ','positive',2.5};
    defaults(end+1,:) = {'minHBratio','nonnegative',1.2}; % peak value over background value
    defaults(end+1,:) = {'minAHratio','nonnegative',.25}; % fitted height over background vs peak value
    defaults(end+1,:) = {'maxUncert','nonnegative',2}; % pixels
pars = ParseVariableArguments(varargin,defaults,mfilename);

%
% This is parent routine for the following sub routines
%  FindSpots 
%  SpotsInNucleus
%% step 4
% dataFolder = 'N:\Derek\20200929_L10-b2-E14mESC_im6\DNA_Expt\'
% eTableXLS = [dataFolder,'ExpLayout_Corrected.xlsx'];
% analysisFolder = [dataFolder,'AnalysisAB\'];
% load experiment table
eTable = readtable(eTableXLS);
nHybs = height(eTable);
nDatChns = length(str2num(eTable.Readouts{1}));  %#ok<*ST2NM>
%%

%---- deal with table input 
if istable(cellBorders)
    cellBorderTable = cellBorders;
    nCells = max(cellBorderTable.cellNum);
    cellBorders = cell(nCells,1);
     for c=1:nCells 
        isC = cellBorderTable.cellNum == c;
        cellBorders{c} = [cellBorderTable.x(isC),cellBorderTable.y(isC)];
    end
end

%-- optional data
nCells = size(cellBorders,1);
cellProj = cell(nCells,nHybs,nDatChns,2);
pars.xyNucWindow = 1;
yTiles = 4;
xTiles = ceil(nCells/yTiles);

nReads = nDatChns*nHybs;
spotTablePerRead = cell(nReads,1);
cropSpotsPerRead = cell(nReads,1);
% start loop
r = 0;
for h = 1:nHybs
    daxInH = LoadDaxFromEtable(eTable,...
            'dataFolder',pars.dataFolder,...
            'fov',fov,'hybNumber',h,...
            'fixDrift',true,'driftFolder',pars.analysisFolder,...
            'dataType','data',...
            'maxProject',false);

    chns = squeeze(daxInH(h,fov,:)); % 
    nChns = length(chns);
    % -------- OPTIONAL, Display image, mostly for troubleshooting 
    if pars.showRawImageFig
        figure(pars.showRawImageFig); clf; 
        for c=1:nChns
            im = max(chns{c},[],3);
            im = IncreaseContrast(im,'high',.99999); 
            subplot(1,nChns,c); imagesc(im);
        end
    end
    
    % ------- Find spots
    for n=1:nChns
        r=r+1;
        im3D = chns{n};
        [xy,~] = ProjectIm3D(im3D); % starting with 2D is faster
        spotXY = FindSpots(xy,'autoSelectThreshold',pars.autoSelectThreshold,'bkdFilterScale',pars.bkdFilterScale,'showPlot',false);
        [spotTablePerRead{r},cropSpotsPerRead{r}] = SpotsInNucleus(im3D,spotXY,cellBorders); % 10% of tot time, mostly on lsqn fitting as expected   
    end
    % 
    % ----- OPTIONAL, for display only  
    for n=1:nChns
        for c=1:nCells
            im3D = chns{n};
            [ys,xs,~] = size(im3D);
            cntrs = mean(cellBorders{c},1); % 
            xc1 = max(1, floor( cntrs(1,1)- pars.xyNucWindow*pars.radi)+1);
            xc2 = min(xs, floor( cntrs(1,1)+ pars.xyNucWindow*pars.radi) );
            yc1 = max(1, floor( cntrs(1,2)- pars.xyNucWindow*pars.radi)+1);
            yc2 = min(ys, floor( cntrs(1,2)+ pars.xyNucWindow*pars.radi) );
            nucCrop = im3D(yc1:yc2,xc1:xc2,:);
            [xy,xz,~] = ProjectIm3D(nucCrop);
            im0 = xy;
            fs = -fspecial('log',8,2);
            imF = imfilter(im0,fs,'conv','symmetric'); % default is correlation not convolution 
            cellProj{c,h,n,1} = imF;
            im0 = xz;
            fs = -fspecial('log',8,2);
            imF = imfilter(im0,fs,'conv','symmetric'); 
            cellProj{c,h,n,2} = imF;
        end
    end
end

% display cell images (filtered overlay data)
figure(pars.xyTileFig);
for c=1:nCells
    im = cat(3,cellProj{c,:,1,1});
    im = Ncolor(im,'contrastRGB',true,'contrastHigh',.99);
    figure(pars.xyTileFig); 
    subplot(yTiles,xTiles,c); 
    imagesc(im);
end

figure(pars.xzTileFig);
for c=1:nCells
    im = cat(3,cellProj{c,:,1,2});
    im = Ncolor(im,'contrastRGB',true,'contrastHigh',.99);
    figure(pars.xzTileFig); 
    subplot(yTiles,xTiles,c); 
    imagesc(im);
end

% Vector image of cell outlines and spots
figure(pars.vecImageFig); clf; % figure(3); clf;
for c=1:length(cellBorders)
    plot(cellBorders{c}(:,1),cellBorders{c}(:,2),'color','k'); hold on;
    text(cellBorders{c}(1,1)+5,cellBorders{c}(1,2),num2str(c),'color','k','fontSize',15); hold on;
end
cmap = GetColorMap('hsvCut',nReads);
for r=1:nReads
    x = spotTablePerRead{r}.x(spotTablePerRead{r}.fineFit);
    y = spotTablePerRead{r}.y(spotTablePerRead{r}.fineFit);
    plot(x,y,'.','color',cmap(r,:),'MarkerSize',10);
end

    