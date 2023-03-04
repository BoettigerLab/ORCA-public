function [spotCodeTable,mapData,polyData] = DemultiplexByPixel(expFolder,analysisFolder,varargin)
%% required inputs
% expFolder - this is  string filepath, which ends in a filesep symbol, 
%  for the folder where the data files are stored. The barcode images for
%  demultiplexing should be in this folder.  The default settings assume
%  these folders begin with "Bar*", e.g. Barcode01, Barcode02. If you named
%  them something different you can specify the "barcodeRoot" flag as
%  something different. 
% analysisFolder - a string filepath (ending with a filesep) to the folder
%  containing the _AllFits.csv output tables from ChrTracer3
% 
%% outputs
% spotCodeTable -- a data table recording for every fiducial spot found by
%   ChrTracer3, what barcode was assigned to it. Also tracks fiducial x,y
%   position and fov number.
% mapData -- a matrix of all the distance maps for all fiducial spots (aka
%   cells).  The matrix size is nR x nR x nC, where nB is the number of
%   readouts in the experiment (or readouts + repeats) and nC is the number
%   of fiducial spots, and notably is equal to the number of rows in the
%   spotCodeTable.
%  polyData -- a matrix of all the polymer trajectories, nB x 3 x nC. 
% 
% 
%% Optional inputs
% 'dataChns' - if you have 3 color images (2 barcodes + 1 fiducial) update
%     this. For example 1:2 would say the first 2 channels in the hyb image
%     are the barcode data.  2:3 would say the second two, etc.  
%     Note: Channels are always sorted red to blue.
% 'fidChnBar' - which channel has the fiducial data in the barcode images. 
% 'fidChnRef' - which channel has the fiducial data in the reference image
% 'bins' - if the code has trouble concatinating your images, it will help
%       tell it how many bins to expect for the distance map. If you
%       specify a number fewer than you use, the map will truncate. If you
%       specify a larger number, the extra rows will be padded with NaN. 
%  'downsample' - speed up processing by downsampling the barcode images by
%      this factor
% 'barHybs' - allows you to specify the full file path where to find the
%       barcode images. Use this to bypass autodetection from the
%       expFolder.
% 'barcodeRoot' - change the root of the barcode folders in the expFolder
% 'daxRoot' - change the root of the daxfiles. default is 'ConvZscan'.
% 
%
%%  Example use
% [spotCodeTable,mapData]= DemultiplexByPixel('path/to/myBarcodeData/','path/to/AllFits/');
% for b=1:4
%   isB = spotCodeTable.spotGroup==b;
%   cMap = ContactFrac(mapData(:,:,isB));
%   subplot(1,4,b); imagesc(cMap);
% end

%% parse optional inputs
defaults = cell(0,3);
defaults(end+1,:) = {'datChns','integer',1}; 
defaults(end+1,:) = {'fidChnBar','integer',2};
defaults(end+1,:) = {'fidChnRef','integer',2};
defaults(end+1,:) = {'bins','integer',0}; 
defaults(end+1,:) = {'downsample','positive',1}; 
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'refFolder','string','Hyb_001'}; 
defaults(end+1,:) = {'daxRoot','string','ConvZscan'};
defaults(end+1,:) = {'barHybs','cell',{}};
defaults(end+1,:) = {'barcodeRoot','string','Bar'};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'overwrite','boolean',false}; % overwrite spotCodeTable
defaults(end+1,:) = {'fovs','nonnegative',inf};
% figs
defaults(end+1,:) = {'saveFigs','boolean',true};
defaults(end+1,:) = {'figMask','integer',0};
defaults(end+1,:) = {'figAlign','integer',30}; 
defaults(end+1,:) = {'figShowBars','integer',31}; 
defaults(end+1,:) = {'figCDF','integer',32}; 
defaults(end+1,:) = {'figFracNA','integer',33};
defaults(end+1,:) = {'plotResults','boolean',true};
defaults(end+1,:) = {'figOverlay','integer',34}; 
defaults(end+1,:) = {'figTile','integer',35}; 
% 
defaults(end+1,:) = {'nppXY','freeType',[]}; % add to over-ride auto detect of pixel size
defaults(end+1,:) = {'codebook','freeType',[]};


% use full 3D stack or just max project.
pars = ParseVariableArguments(varargin,defaults,mfilename);


%% main function
if isempty(pars.saveFolder)
    pars.saveFolder = SetFigureSavePath([analysisFolder,'deMultiplex',filesep],'makeDir',true);
else
    pars.saveFolder = SetFigureSavePath(pars.saveFolder,'makeDir',true);
end

sc = pars.downsample; % downsample rescale images for faster processing


%% check if data already exists
if exist([pars.saveFolder,'spotCodeTable.csv'],'file') && pars.overwrite == false
    skip = true;
    spotCodeTable = readtable([pars.saveFolder,'spotCodeTable.csv']);
else
    skip = false;
end

%%
if skip  % just load polys and maps
    [polys,maps] = CombineAllFits(analysisFolder,'byFOV',true); 
else
    refDaxFiles = FindFiles([expFolder,pars.refFolder,filesep,pars.daxRoot,'*.dax']);
    imProps = LoadDax([refDaxFiles{1}],'justImProps',true);
    if isempty(pars.nppXY)
        nppXY = imProps.xy2um*1e3;
    else
        nppXY = pars.nppXY;
    end
    if isempty(pars.barHybs)
        barHybs = FindFiles([expFolder,pars.barcodeRoot,'*'],'onlyFolders',true);
    else
        barHybs = pars.barHybs;
    end
    B = length(barHybs);
    barDaxFiles = cell(B,1);
    for b=1:B % b=2
        barDaxFiles{b} = FindFiles([barHybs{b},filesep,pars.daxRoot,'*.dax']);
    end

    if length(barDaxFiles{1}) ~= length(refDaxFiles)
        error('different number of FOVs in refHyb and barHyb');
    end

    if isinf(pars.fovs)
        nFOV = length(barDaxFiles{1});
        fovs = 1:nFOV;
    else
        fovs = pars.fovs;
        nFOV = length(fovs);
    end
    spotCodes = cell(nFOV,1);
    polys = cell(nFOV,1);
    maps = cell(nFOV,1);
    imBarStk = cell(nFOV,1); 
    for f=fovs
        % --- load ORCA AllFits table
        try
            orcaTable = readtable([analysisFolder,'fov',num2str(f,'%03d'),'_AllFits.csv']);
            [polys{f},maps{f},spotData] = TableToPolymer(orcaTable,'bins',pars.bins);         
            spotXY = round(spotData/nppXY/sc); 
        catch
            warning(['failed to load AllFits for fov ',num2str(f)]);
            continue
        end
        % -- Load data, ref first
        refDax = LoadDax([refDaxFiles{f}],'maxProject',true);
        refDax = imresize(refDax,1/sc);
        barIms = cell(B*length(pars.datChns),1);
        barMask= cell(B*length(pars.datChns),1);
        k = 0; % a counter for the total barcode images
        for b=1:B  
            % -- now load each barcode hyb
            [barDax,imBarData] = LoadDax([barDaxFiles{b}{f}],'maxProject',true,'verbose',pars.veryverbose);
            barDax = imresize(barDax,1/sc);
            % -- align images to ref Hyb
            % ---- check if previous alignValues have been saved  (if exist, default to skip)
            if pars.figAlign
                f30 = figure(pars.figAlign); clf;
            end
            alignV = CorrAlignFast(refDax(:,:,pars.fidChnRef),barDax(:,:,pars.fidChnBar),'showplot',pars.figAlign>0,'showExtraPlot',false); pause(.01);
            % ---- save alignValues (or load them if uses previous data)
            %   ADD CODE HERE
            % ---- apply alignments 
            for d=pars.datChns
                k=k+1;
                barIms{k} = ApplyReg(barDax(:,:,d),alignV);
                % build mask  [Not currently used. The idea was to use the mask
                % to group pixels outside beyond just the overlapped fid spot
                % to classify better. Neighboring pixels contribute brightness
                % in any case, and in the embryo the seperation is better also,
                % so this appears not needed.]
                barMask{k} = b*uint8(imbinarize(barIms{k})); % uses Otsu' method. 
                % show mask
                im = labeloverlay(barIms{k},logical(barMask{k}),'Transparency',.5);
                if pars.figMask
                    figure(pars.figMask); clf; imagesc(im);
                    hold on; plot(spotXY(:,1),spotXY(:,2),'yo');
                    pause(.1);
                end
            end

            if pars.saveFigs && pars.figAlign
                SaveFigure(f30,'name',['alignFig_f',num2str(f,'%03d'),'_b',num2str(b)],'formats',{'png'},'overwrite',pars.overwrite);
            end
            imBarStk{f} = barIms;
        end
        if pars.figMask % mostly for troubleshooting
            f31 = figure(pars.figMask);  clf; 
            Ncolor(cat(3,barIms{:})); 
            hold on; plot(spotXY(:,1),spotXY(:,2),'yo');
            pause(.1);
            if pars.saveFigs
                SaveFigure(f31,'name',['figMask_f',num2str(f,'%03d')],'formats',{'png'},'overwrite',pars.overwrite);
            end
        end



        nSpots = size(spotXY,1);
        spotCodes{f} = zeros(nSpots,3+B);
        nB = size(barIms,1);
        for b=1:nB        
             % record spot brightness
             [nRows,nCols] = size(barIms{b});
             for s=1:nSpots
                 xi = spotXY(s,1);
                 yi =  spotXY(s,2);
                 xi = max([min([xi,nRows]),1]);
                 yi = max([min([yi,nCols]),1]);
                 if ~isnan(xi)
                     spotCodes{f}(s,1) = xi;
                     spotCodes{f}(s,2) = yi;
                     spotCodes{f}(s,3) = f; 
                     spotCodes{f}(s,3+b) = barIms{b}(yi,xi);
                 end
             end
        end
    end


    
    %% Balance barcodes
    % after we go through all FOVs, we have a better idea the typical ranges of
    % the barcodes used. We then 

    codeData = cat(1,spotCodes{:});
    codeScore = codeData(:,4:end);
    if pars.figCDF
        f32 = figure(pars.figCDF); clf;
        for b=1:nB
            cdfplot(codeScore(:,b)); hold on;
        end
        legend( cellstr(num2str( (1:nB)' )));
        if pars.saveFigs
            SaveFigure(f32,'name',['barcodeCDFs_f',num2str(f,'%03d')],'formats',{'png','fig'},'overwrite',pars.overwrite);
        end
    end

    totCells = size(codeScore,1);
    mL = quantile(codeScore,.01,1);
    mB = quantile(codeScore,.99,1);
    codeScore = codeScore./ repmat(mB,totCells,1);
    % codeScore = codeScore.* repmat(mB-mL,totCells,1);
    codeScore = codeScore./ repmat( max(codeScore,[],2),1,nB); 
    if pars.figShowBars
        f33 = figure(pars.figShowBars); clf; 
        imagesc(codeScore(randi(totCells,200,1),:)); 
        colorbar; 
        if pars.saveFigs
            SaveFigure(f33,'name',['exampleBarcodes_f',num2str(f,'%03d')],'formats',{'png','fig'},'overwrite',pars.overwrite);
        end
    end

    %% decode

    if ~isempty(pars.codebook)
        if ischar(pars.codebook)
            codebook = readtable(pars.codebook);
        else
            codebook = pars.codebook;
        end
        codebookBarcodes = codebook{1,2:end}; % record the barcode numbers used
        codebookNames = codebook{2:end,1};
        codeMatrix = codebook{2:end,2:end};
    else
        warning('no codebook provided. Assuming each barcode is a unique cell type');
        codeMatrix = eye(nB);
        codebookBarcodes = 'unknown';
        codebookNames = cellstr(repmat('unknown',nB,1));
    end


    spotGroup = zeros(totCells,1);
    contrast = zeros(totCells,1);
    groupName = cell(totCells,1);
    nearestGroup = zeros(totCells,1);
    for c=1:totCells % c = 200
        [idx,dis] = knnsearch(codeMatrix,codeScore(c,:),'K',2);
        spotGroup(c) = idx(1);
        nearestGroup(c) = idx(2);
        contrast(c) = (1/dis(1)) ./ (1/dis(2));
        groupName{c} = codebookNames{idx};
    end

    x = codeData(:,1)*sc;
    y = codeData(:,2)*sc;
    fov = codeData(:,3);
    spotCodeTable = table(spotGroup,nearestGroup,contrast,x,y,fov);

    writetable(spotCodeTable,[pars.saveFolder,'spotCodeTable.csv']);
    disp(['wrote ',pars.saveFolder,'spotCodeTable.csv']);
    
    if pars.figFracNA
        f34 = figure(pars.figFracNA); clf; 
        histogram(contrast,0:.1:10);
        fracNonAmbig = num2str(sum(contrast>1.5)./ sum(contrast>0)*100,3);
        title(['non-ambiguous spots =',fracNonAmbig,'%'])
        if pars.saveFigs
            SaveFigure(f34,'name',['fracNonAmbig_f',num2str(f,'%03d')],'formats',{'png','fig'},'overwrite',pars.overwrite);
        end
    end
    
    
   %% Plot results

    if pars.plotResults
        nG = length(unique(spotCodeTable.spotGroup));
        figF = figure(pars.figOverlay); clf;
        subF = gca;
        figT = figure(pars.figTile); clf;
        for fov = fovs
            if pars.saveFigs
               figName = [pars.saveFolder,'fov',num2str(fov,'%03d'),'_Demultiplex.png'];
               if exist(figName,'file') && ~pars.overwrite || isempty(imBarStk{fov})
                  continue 
               end
            end        
            subF.NextPlot = 'replace';
            cMap = GetColorMap('hsvCut',nB); % barcodes
            gMap = GetColorMap('hsvCut',nG); % groups
            ncImage = imresize(cat(3, imBarStk{fov}{:}),1/sc);
            for b=1:nB
                ncImage(:,:,b) =uint16( 2^15*double(ncImage(:,:,b))./double(mB(b)) );
            end
            % the overlay figure
            figure(figF);
            Ncolor(ncImage,'colormap',cMap);
            colormap(cMap);
            colorbar;
            isIn = spotCodeTable.fov==fov;
            xyA = [spotCodeTable.x(isIn),spotCodeTable.y(isIn)]; % spotCodeTable{isIn,1:2};
            hold on; plot(xyA(:,1),xyA(:,2),'k.','MarkerSize',20)
            hold on; plot(xyA(:,1),xyA(:,2),'w.')
            for g=1:nG  
                isIn = spotCodeTable.fov==fov & spotCodeTable.spotGroup==g;
                isLow = isIn & spotCodeTable.contrast<1.5;
                xy = [spotCodeTable.x(isIn),spotCodeTable.y(isIn)]; % {isIn,1:2};
                xyL =[spotCodeTable.x(isLow),spotCodeTable.y(isLow)]; %  spotCodeTable{isLow,1:2};
                subF.NextPlot = 'add';
                plot(subF,xy(:,1),xy(:,2),'.','color',gMap(g,:),'MarkerSize',16);
                plot(subF,xyL(:,1),xyL(:,2),'+','color',gMap(g,:),'MarkerSize',16);
            end
            % the tile figure
            figure(figT); clf;
            subHandles = TileImageStack(ncImage,'colormap',cMap,'numRows',2);
            for sH=1:length(subHandles)    
                hold on; plot(subHandles{sH},xyA(:,1),xyA(:,2),'k.','MarkerSize',14)
                hold on; plot(subHandles{sH},xyA(:,1),xyA(:,2),'w.')
                subHandles{sH}.NextPlot = 'add';
                 for g=1:nG  
                    isIn = spotCodeTable.fov==fov & spotCodeTable.spotGroup==g;
                    isLow = isIn & spotCodeTable.contrast<1.5;
                    xy = [spotCodeTable.x(isIn),spotCodeTable.y(isIn)]; % {isIn,1:2};
                    xyL =[spotCodeTable.x(isLow),spotCodeTable.y(isLow)]; %  spotCodeTable{isLow,1:2};
                    subHandles{sH}.NextPlot = 'add';
                    plot(subHandles{sH},xy(:,1),xy(:,2),'.','color',gMap(g,:),'MarkerSize',10);
                    plot(subHandles{sH},xyL(:,1),xyL(:,2),'x','color',gMap(g,:),'MarkerSize',5);
                end
            end

            if pars.saveFigs
               figName = ['fov',num2str(fov,'%03d'),'_Demultiplex'];
               SaveFigure(figF,'name',figName,'formats',{'png'},'overwrite',pars.overwrite); 
                figName = ['fov',num2str(fov,'%03d'),'_DemultiplexTiled'];
               SaveFigure(figT,'name',figName,'formats',{'png'},'overwrite',pars.overwrite); 
            end
            pause(.01);
        end 
    end
    
    
end

%% flatten polymer and map data
try
polyData = CatPolys(polys); %  cat(3,polys{:});
mapData = CatMaps(maps);
catch
    polyData = polys;
    mapData = maps;
    warning('failed to concatinate maps');
end
%%

