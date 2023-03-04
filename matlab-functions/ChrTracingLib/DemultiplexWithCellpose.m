function [spotCodeTable,mapData,polyData,codeData] = DecodeSampleBarcodes(eTableXLS,analysisFolder,varargin)

 % Refs
 % https://cellpose.readthedocs.io/en/latest/command.html
 
defaults = cell(0,3);
defaults(end+1,:) = {'scope','string','auto'};
defaults(end+1,:) = {'nppXY','freeType',[]}; % will read from table eTableXLS info file by default
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'veryverbose','boolean',false}; 
defaults(end+1,:) = {'minFracOn','fraction',.05}; 
defaults(end+1,:) = {'overwrite','boolean',false}; 
defaults(end+1,:) = {'recordQuality','boolean',false}; 
defaults(end+1,:) = {'codebook','string',''};
defaults(end+1,:) = {'bins','integer',0}; % 0 for autodetect
defaults(end+1,:) = {'saveFolder','string',''}; % will create a default folder if blank. This is where the cellpose output and processed table will be saved.
defaults(end+1,:) = {'nucImSize','array',[512,512]}; % for cellpose, we downsample images to this size to run faster
defaults(end+1,:) = {'nucDiameter','nonnegative',0}; % 0 = autodetect  a fixed diameter is substantially faster, and may improve uniformity. 
defaults(end+1,:) = {'refHyb','integer',0}; % By default will try to guess using the regData.csv.  This is the hyb in expTable.xlsx used as alignment reference.  We will use this hyb for matching the spots into their segmented cells
defaults(end+1,:) = {'maxProject','boolean',true}; % use full 3D stack or just max project.
defaults(end+1,:) = {'nucleiImageRoot','string','fidMax_ConvZscan'}; % filename root of 2D images for cellpose
defaults(end+1,:) = {'contrastHigh','fraction',.995}; % contrast for cellpose
defaults(end+1,:) = {'contrastLow','fraction',.0}; % contrast for cellpose

% plotting
defaults(end+1,:) = {'saveFigure','boolean',true}; 
defaults(end+1,:) = {'plotResults','boolean',true}; 
defaults(end+1,:) = {'showExtraPlots','boolean',false}; 
defaults(end+1,:) = {'figShowLoadedImages','integer',1};
defaults(end+1,:) = {'figContrastHist','integer',2};
defaults(end+1,:) = {'showFOVfig','integer',0}; 
defaults(end+1,:) = {'overlayFig','integer',10}; 
defaults(end+1,:) = {'tileFig','integer',11}; 
% python path stuff
defaults(end+1,:) = {'cellpose_env','string','mlab_cellpose'}; 

pars = ParseVariableArguments(varargin,defaults,'DecodeSampleBarcodes');


if isempty(pars.saveFolder) 
    saveFolder = SetFigureSavePath([analysisFolder,'cellposeFitsGPU\'],'makeDir',true);
else
    saveFolder = SetFigureSavePath(pars.saveFolder,'makeDir',true);
end
% check if results already exist, if so, we load these, unless 'overwrite'
% is specified.  If 'overwrite' is true, we repeat the analysis and
% overwrite any existing previous analysis.
% 
spotCodeTableCsv = [saveFolder,'spotCodeTable.csv'];
if exist(spotCodeTableCsv,'file') && ~pars.overwrite
    spotCodeTable = readtable(spotCodeTableCsv);
    runDemultiplex = false;
    if pars.verbose
       disp('detected existing spotCodeTable.csv, Demultiplex will use this file. Specify overwrite=1 to overide'); 
    end
else
    runDemultiplex = true;
end

polys = {};
maps = {}; 

if runDemultiplex
    % determine if cellpose has already been run 
    cellPoseFiles = FindFiles([saveFolder,'*_cp_masks.tif']);
    if isempty(cellPoseFiles) || pars.overwrite
        runCellpose = true;
    else
        runCellpose = false;
    end
    % load the experiment table
    eTable = readtable(eTableXLS);


    if runCellpose
        if pars.verbose
           disp('saving images to run in cellpose'); 
        end
        %%  Save downsampled png images to feed to cellpose for segmentation
        % this version uses 2D projections for spped 
        
        % determine which hyb was the refHyb
        if pars.refHyb == 0
            regData = FindFiles([analysisFolder,'fov*_regData.csv']);
            regTable = readtable(regData{1});
            pars.refHyb = find(regTable.xshift==0 & regTable.yshift==0  & regTable.xshift2==0  & regTable.yshift2==0);  
            pars.refHyb = pars.refHyb(1); % enforce single
        end
        
        root = pars.nucleiImageRoot;
        dataFolder = [fileparts(eTableXLS),'\'];
        nucleiFolder = [dataFolder,eTable.FolderName{pars.refHyb}];
        fidMaxFiles = FindFiles([nucleiFolder,'\',root,'*.dax']);
        if isempty(fidMaxFiles)
            error(['no dax files found in ', nucleiFolder]);
        end
        % hybFolders = FindFiles([dataFolder,'Hyb*']);
        % fidMaxFiles = FindFiles([hybFolders{1},'\fidMax_ConvZscan*.dax'])
        % fidMaxFiles = FindFiles([hybFolders{1},'\dat647Max_ConvZscan*.dax'])
        nFOV = length(fidMaxFiles);
        for f=1:nFOV
            im1 = ReadDax(fidMaxFiles{f},'verbose',pars.veryverbose);
            [~,imName] = fileparts(fidMaxFiles{f});
            imOut = makeuint(imresize(im1,pars.nucImSize),8);
            imOut = IncreaseContrast(imOut,'high',pars.contrastHigh,'low',pars.contrastLow);
           if pars.figShowLoadedImages
                figure(pars.figShowLoadedImages); clf;
                imagesc(imOut); pause(.01); colormap(gray);
           end
          imwrite(imOut,[saveFolder,imName,'.png']);
        end
        sc = size(im1,1)./pars.nucImSize(1); % rescaling factor

        %%  Run Cellpose!
        if pars.verbose
            disp('running cellpose to ID nuclei, please wait...');
        end
        % this is 2D, since we're running on 2D images based on above 
        %  a fixed diameter is substantially faster, and may improve uniformity. 
        tic
        diameter = num2str(pars.nucDiameter);
        callPython = ['python -m cellpose --dir ',saveFolder,' --pretrained_model nuclei --diameter ',diameter,' --save_tif --no_npy --use_gpu'];
        setEnv = ['! activate ',pars.cellpose_env, ' '];% '! activate mlab_cellpose '; % '! activate cellpose ';
        cmdOut = [setEnv,' && ',callPython];
        if pars.verbose
            disp(cmdOut);
        end
        eval(cmdOut); % should move to system Run 
        toc
    end

    %%
    if pars.verbose
        disp('loading barcode data...');
    end
    % this handles the dift correction assuming the eTableXLS includes the
    % barcode hybes and that this eTableXLS was used in running ChrTracer
    % (rather than a different one ommitting the barcode hybs). 
    [imMax,imInfo] = LoadDaxFromEtable(eTableXLS,...
        'hybNumber',inf,...
        'fov',inf,...
        'hybType','B',...
        'dataType','data',...
        'fixDrift',true,...
        'maxProject',pars.maxProject,...
        'driftFolder',analysisFolder,...
        'verbose',pars.verbose);
    imMax0 = imMax; % backup for debugging only.
    % remove empty spaces for skipped hybs
    %     if barcodes are the last N hybs, there will be H blank cells which
    %     were skipped for the readout hybs prior to the cellular barcodes. 
    noData = cellfun(@isempty,imMax);
    imMax(noData(:,1),:,:) = [];
    imInfo(noData(:,1),:,:) = [];
    % stack barcodes
    %   This step flattens out the multi-color data (if there were multiple
    %   barcodes read out per hyb).  Barcodes are stacked in the order they
    %   were acquired, following the convention used in ChrTracer3(). 
    %   Example: red-1, green-1, red-2, green-2, ...
    %   This makes our nHyb x nFOV x nDatChns cell array into a new array which
    %   has dimensions nHyb*nDatChns x nFOV.
    [nH,nFOVs,nDatChns] = size(imMax);
    if nDatChns>1
       imMaxStk = cell(nH*nDatChns,nFOVs);
       k=0;
       for h=1:nH
           for c=1:nDatChns
               k=k+1;
              imMaxStk(k,:) = imMax(h,:,c);
           end
       end
       imMax = imMaxStk;
    end
    if size(imMax,1) <= 1
       warning('barcode folders must be designated as dataType "B" in the experiment table');
       error('Did not find multiple barcode folders, check your experiment table'); 
    end
    % determine scope settings from info file of first image
    if isempty(pars.nppXY)
        scopePars = GetScopeSettings(imInfo{1,1,1});
        pars.nppXY = scopePars.nmPixXY;
    end

    isB = strcmp(eTable.DataType,'B');
    datPropTable = DataChnsFromTable(eTable(isB,:));
    reads = datPropTable.readout;
    idZero = find(reads==0);
    imMax(idZero,:) = []; %#ok<FNDSB>

    [nB,nFOV] = size(imMax);
    %% Load the mask and score barcode values by cell
    % First we load the masks and we load the ChrTracer data
    % Then we match spots & traces to Cell IDs
    % Then we match cell IDs with barcode images and record the value in each
    % cell.  We'll wait until all FOV have been processed before we try
    % correcting barcode brightness variation and assigning cells based on the
    % codebook.
    if pars.verbose
        disp('loading ORCA data and cellpose results...');
    end
    tic
    allMasks = FindFiles([saveFolder,'*_cp_masks.tif']);
    polys = cell(nFOV,1);
    maps = cell(nFOV,1);
    cellMasks = cell(nFOV,1);
    traceTables = cell(nFOV,1);
    cellTables = cell(nFOV,1);
    currCellTotal = 0;
    currTraceTotal = 0;
    for f=1:nFOV %    f=2
        if pars.veryverbose
           disp(['loading FOV ',num2str(f), ' of ',num2str(nFOV)]); 
        end
        mask1 = imread(allMasks{f}); % cell ID#s
        cellMasks{f} = mask1;
        cellProps = regionprops(mask1,'Centroid','PixelIdxList');
        sc = scopePars.imHeight./pars.nucImSize(1); % rescaling factor.  (redefined here, in case runCellpose was not recalled)
        nCells = length(cellProps);
        if pars.figShowLoadedImages
            figure(pars.figShowLoadedImages); clf; 
            imagesc(mask1); GetColorMap('distColors',nCells);
            title(['FOV=',num2str(f), ' of ',num2str(nFOV)]);
        end
         % --- load ORCA AllFits table    
         orcaTable = readtable([analysisFolder,'fov',num2str(f,'%03d'),'_AllFits.csv']);
         [polys{f},maps{f},spotData] = TableToPolymer(orcaTable,'bins',pars.bins);  % why is this hanging sometimes? this was supposed to be faster
         % assign cell ID based on the ID of the underlying pixel in the
         % segmented map. 
         traceXY = round(spotData/pars.nppXY/sc); 
         nSpots = size(traceXY,1);
         [nRows,nCols] = size(mask1);
         cellIDfov = zeros(nSpots,1);
         for s=1:nSpots
             xi = traceXY(s,1);
             yi = traceXY(s,2);
             xi = max([min([xi,nRows]),1]);
             yi = max([min([yi,nCols]),1]);
             cellIDfov(s) = mask1(yi,xi);
         end
        % store fov data in fov table for traces
        traceID = (1:nSpots)' + currTraceTotal;
        traceX = traceXY(:,1); % position is relative to fixed 512x512 barcode image
        traceY = traceXY(:,2); % position is relative to fixed 512x512 barcode image
        cellUID = cellIDfov + currCellTotal;
        traceTables{f} = table(traceID,traceX,traceY,cellIDfov,cellUID);
        currTraceTotal = currTraceTotal + nSpots;
        
        % get the full barcode pattern in each cell  
        cellBarcodeAves = zeros(nCells,nB);
        for c=1:nCells
            for b=1:nB 
                imSmall = imresize(imMax{b,f},1./sc) ;
                cellBarcodeValues = imSmall(cellProps(c).PixelIdxList);
                cellBarcodeAves(c,b) = median(cellBarcodeValues);
            end
        end       
        % store cell data in table for cells
        if nCells > 0
            cellIDfov_c = (1:nCells)';
            cellUID_c = cellIDfov_c + currCellTotal; % need unique names,
            cellXY = cat(1,cellProps(:).Centroid);        
            cellX = cellXY(:,1); % position is relative to fixed 512x512 barcode image
            cellY = cellXY(:,2);  % position is relative to fixed 512x512 barcode image
            fov = f*ones(nCells,1);
            barcodeValues = array2table(cellBarcodeAves);
            cellTable_f = table(cellUID_c,cellIDfov_c,cellX,cellY,fov);
            cellTable_f = cat(2,cellTable_f,barcodeValues);
            cellTables{f} = cellTable_f;
            currCellTotal = currCellTotal + nCells;
        end
        % combine info into tables
    end
    toc
    
    cellTable  = cat(1,cellTables{:});
    traceTable = cat(1,traceTables{:});
    barcodeIntensity = cellTable{:,6:end};
    totCells = size(barcodeIntensity,1);
    unassignedSpots = traceTable.cellIDfov==0;

    %% Decode
    % minFracOn = .05;% minimum fraction observed in each barcode
    % normalize the brightness of the relative barcodes, then compare locally
    % within the cell.  We need to normalize against when this barcode is 'ON'.
    % We assume all barcodes are ON in at least 1% of cells, and the variation
    % in the ON fraction is substantially less than between the ON and OFF
    % populations, so as long as 1% lands us in the ON fraction (whether 98%
    % were on or only 1% were on), it is an appropriate norm.
    mB = quantile(barcodeIntensity,1-pars.minFracOn,1);
    cellNormBarcodeBright = barcodeIntensity./repmat(mB,totCells,1);  % balance barcode brightnesses
    figure(2); clf; imagesc(cellNormBarcodeBright(randi(totCells,100,1),:)); colorbar; 
    % codeScore = cellNormBarcodeBright./ repmat( sum(cellNormBarcodeBright,2),1,nB); % scale all brightnesses to sum to 1
    codeScore = cellNormBarcodeBright./ repmat( max(cellNormBarcodeBright,[],2),1,nB); 
      figure(2); clf; imagesc(cellNormBarcodeBright(randi(totCells,100,1),:)); colorbar; 
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
    
    % save some additional information we may want to keep track of. 
    codeData.codebookBarcodes = codebookBarcodes;
    codeData.codebookNames = codebookNames;
    codeData.codeMatrix = codeMatrix;
    codeData.scopePars = scopePars;
    codeData.scaleFactor = sc;

    spotGroup = zeros(totCells,1);
    contrast = zeros(totCells,1);
    groupName = cell(totCells,1);
    nearestGroup = zeros(totCells,1);
    for c=1:totCells
    %         c = randperm(totCells,1);
    %         codeScore(c,:)
        [idx,dis] = knnsearch(codeMatrix,codeScore(c,:),'K',2);
        spotGroup(c) = idx(1);
        nearestGroup(c) = idx(2);
        contrast(c) = (1/dis(1)) ./ (1/dis(2));
        groupName{c} = codebookNames{idx};
    end
    
    figure(pars.figContrastHist); clf; 
    histogram(contrast,0:.1:10);
    fracNonAmbig = num2str(sum(contrast>1.5)./ sum(contrast>0)*100,3);
   title(['non-ambiguous spots =',fracNonAmbig,'%'])
    %% combine tables
    newCellCols = table(spotGroup,nearestGroup,contrast,groupName);
    cellTableFull = cat(2, cellTable,newCellCols);

    % we add a blank row to the cell table, to map all the spots which fall
    % outside annotated cells. 
    blankRow = table('Size',[1,size(cellTableFull,2)],...
        'VariableTypes',varfun(@class,cellTableFull,'OutputFormat','cell'),...
        'VariableNames',cellTableFull.Properties.VariableNames);
    cellTableFull = cat(1,cellTableFull,blankRow);
    
    totCells = height(cellTableFull); 
    traceTable.cellUID(traceTable.cellIDfov==0) = totCells;
    cell2traceTable = cellTableFull(traceTable.cellUID,:);
    spotCodeTable = cat(2,traceTable,cell2traceTable);
    

    % display a randomly selected 10 rows of the spotCodeTable
    totTraces = height(spotCodeTable);
    spotCodeTable(randi(totTraces,10,1),:)
    
    % remove sanity check columns
    checkCol1 = strcmp(spotCodeTable.Properties.VariableNames,'cellUID_c');
    checkCol2 = strcmp(spotCodeTable.Properties.VariableNames,'cellIDfov_c');
  %   spotCodeTable(:,checkCol1 | checkCol2) = [];

    spotCodeTable(randi(totTraces,10,1),:)
   
    
    %% save results
    writetable(spotCodeTable,spotCodeTableCsv);
    save([saveFolder,'codeData.mat'],'codeData');
    if pars.verbose
        disp(['saved ',spotCodeTableCsv]);
    end
    %% Plot results

    if pars.plotResults
        nG = length(unique(spotCodeTable.spotGroup));
        figF = figure(pars.overlayFig); clf;
        subF = gca;
        figT = figure(pars.tileFig); clf;
        for fov = 1:nFOV
            if pars.saveFigure
               figName = [saveFolder,'fov',num2str(fov,'%03d'),'_Demultiplex.png'];
               if exist(figName,'file') && ~pars.overwrite
                  continue 
               end
            end        
            subF.NextPlot = 'replace';
            cMap = GetColorMap('hsvCut',nB); % barcodes
            gMap = GetColorMap('hsvCut',nG); % groups
            ncImage = imresize(cat(3, imMax{:,fov}),1/sc);
            for b=1:nB
                ncImage(:,:,b) =uint16( 2^15*double(ncImage(:,:,b))./double(mB(b)) );
            end
            % the overlay figure
            figure(figF);
            imOut = Ncolor(ncImage,'colormap',cMap);
            mask = boundarymask(cellMasks{fov}); 
            imOut =labeloverlay(imOut,mask,'Transparency',0);
            imagesc(imOut);
            colormap(cMap);
            colorbar;
            isIn = spotCodeTable.fov==fov;
            xyA = [spotCodeTable.traceX(isIn),spotCodeTable.traceY(isIn)]; % spotCodeTable{isIn,1:2};
            hold on; plot(xyA(:,1),xyA(:,2),'k.','MarkerSize',20)
            hold on; plot(xyA(:,1),xyA(:,2),'w.')
            for g=1:nG  
                isIn = spotCodeTable.fov==fov & spotCodeTable.spotGroup==g;
                isLow = isIn & spotCodeTable.contrast<1.5;
                xy = [spotCodeTable.traceX(isIn),spotCodeTable.traceY(isIn)]; % {isIn,1:2};
                xyL =[spotCodeTable.traceX(isLow),spotCodeTable.traceY(isLow)]; %  spotCodeTable{isLow,1:2};
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
                    xy = [spotCodeTable.traceX(isIn),spotCodeTable.traceY(isIn)]; % {isIn,1:2};
                    xyL =[spotCodeTable.traceX(isLow),spotCodeTable.traceY(isLow)]; %  spotCodeTable{isLow,1:2};
                    subHandles{sH}.NextPlot = 'add';
                    plot(subHandles{sH},xy(:,1),xy(:,2),'.','color',gMap(g,:),'MarkerSize',10);
                    plot(subHandles{sH},xyL(:,1),xyL(:,2),'x','color',gMap(g,:),'MarkerSize',5);
                end
            end
            
            if pars.saveFigure
               figName = ['fov',num2str(fov,'%03d'),'_Demultiplex'];
               SaveFigure(figF,'name',figName,'formats',{'png'},'overwrite',pars.overwrite); 
                figName = ['fov',num2str(fov,'%03d'),'_DemultiplexTiled'];
               SaveFigure(figT,'name',figName,'formats',{'png'},'overwrite',pars.overwrite); 
            end
            pause(.01);
        end
    end
elseif nargout > 1  % files already exist, we just need to load
    orcaTables = FindFiles([analysisFolder,'fov*_AllFits.csv']);
    nFOV = length(orcaTables);
    polys = cell(nFOV,1);
    maps = cell(nFOV,1);
    for f=1:nFOV
         orcaTable = readtable(orcaTables{f});
         [polys{f},maps{f}] = TableToPolymer(orcaTable,'bins',pars.bins);
    end
    try 
        load([saveFolder,'codeData.mat'],'codeData');
    catch
        warning('error loading codeData.mat');
        codeData = [];
    end
end


% flatten polymer and map data
 try
    polyData = cat(3,polys{:});
    mapData = cat(3,maps{:});
 catch er
    warning(er.message);
    warning('polyData and mapData are cell arrays by fov');
    polyData = polys;
    mapData = maps; 
    
 end
