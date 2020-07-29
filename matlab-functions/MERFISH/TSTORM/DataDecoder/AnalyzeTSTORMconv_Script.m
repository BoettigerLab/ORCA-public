
% TSTORM
% Analyze multi-day STORM images of sequential staining


%% Parsing filepaths:
startup;
 
codebookPath = '\\cajal\TSTORMdata\Library7\Lib7_140808\';
TSTORMdata = '\\cajal\TSTORMdata\' ;
codebookFile = [codebookPath,'E1_codebook.fasta'];
dataPath = [TSTORMdata,'140826_L7E1_conv\'];
fpkmPath = [TSTORMdata,'GenomeData\IMR90\sequencing data from Steven\'];
listType = 'alist.bin';

reverseBits = false;
showCellImages = true;
saveCellImages = true;

savePath = [dataPath,'Analysis3\'];
if ~exist(savePath,'file')
    mkdir(savePath);
end

hybe1Data = dir([dataPath,'STORM_01','_*_','c1_',listType]);
codebook = fastaread(codebookFile);

numGenes = length(codebook);  % number of genes (including neg controls) 
numHybes = length(regexprep(codebook(1).Header,' ','')); % number of stains 
numCells = length(hybe1Data);

if numCells < 1
    error(['No cells found in ',dataPath,' of type ',listType]);
end

rnaCnts = zeros(numGenes,numCells);
perCnts = zeros(numGenes,numCells);
spotData = cell(numCells,1); 
spotDaxData = cell(numCells,1); 

colordef white;

%% Filters

%-----------------------------------------------------------------------
% choose a colormap 

clrmap = jet(numHybes);

% Drift Correction Parameters:
xshift = zeros(1,numHybes);
yshift = zeros(1,numHybes);

% Parameters for aligning hybes
maxD = 8; % max distance over which to attempt to match beads
warp2Hybe1 = false; 

% Parameters for Filtering STORM data
minPhotons =   1; % updated  ;  % min # photons per spot
STORMstartingFrame = 1; 

% Parameters for Clustering Localizations
binSize = 50;  % size of bins in nm
minDotPerBin = 1;  % min number of localizations to call a bin occupied
minLocsPerDot = 1; % min number of localization in all bins assigned to a cluster to be called an mRNA
minArea = 0;   % min area in bins to be called a cluster of localization
maxArea = 10; % max area in bins to be called a cluster of localization

% Parameters for Assigning stains to mRNA centroids
% min number of localizations within MaxDtoCentroid to be called 'detected in stain d'
minPhotsPerStain =  ones(1,numHybes);  %  [ones(1,3),3E3*ones(1,7),8E3*ones(1,2) ]; % 
maxDtoCentroid = 1; % max distance you can be from the centroid of the mRNA

for c = 1:numCells; % c=3

    try
        if numCells>=100
            cellnum = sprintf('%03d',c-1);
        elseif numCells>=10 
            cellnum = sprintf('%02d',c-1);
        else
            cellnum = sprintf('%01d',c-1);
        end

        mRNAdata = dir([dataPath,'*STORM*','_',cellnum,'_','c1_',listType]);
        mRNAbin = strcat(dataPath,{mRNAdata.name}');

        beaddata = dir([dataPath,'*STORM*',cellnum,'_','c2_list.bin']);
        beads = strcat(dataPath,{beaddata.name}');

            % HYBES DONE IN REVERSE, FLIP THESE
        if reverseBits
            mRNAbin = flipud(mRNAbin);
            beads = flipud(beads); 
        end

    %     % Remove bit 0 from list of hybes 
    %     bit0 = StringFind(mRNAbin,'STORM_0_');
    %     mRNAbin(bit0) = []; 
    %     beads(bit0) = []; 

        if isempty(beads) || length(beads) > numHybes ||  length(beads) < numHybes
            warning('some bead data is missing: ');
            disp(beads);
            error(['error in automatic detection of bead data for cell ',cellnum,'!']);
        end
        if isempty(mRNAbin) || length(mRNAbin) > numHybes ||  length(mRNAbin) < numHybes
            warning('some STORM data is missing: ');
            disp(mRNAbin);
            error(['error in automatic detection of STORM data for cell ',cellnum,'!']);   
        end

        beadImage = regexprep(beads{1},'_list.bin','.dax'); % we plot the bead localizations ontop of the bead image from hybe 1.  



        %% ------------------------------  Read Data -----------------------

        imLists = cell(numHybes,1); 
        fedPos = cell(numHybes,1); 
        hybeLocs = zeros(numHybes,1); 
        hybeBrightness = zeros(numHybes,1); 
        for d=1:numHybes;
            imLists{d}   = ReadMasterMoleculeList(mRNAbin{d}); 
            hybeLocs(d) = length(imLists{d}.xc);
            hybeBrightness(d) = median(imLists{d}.a);

            beadList = ReadMasterMoleculeList(beads{d}); 
            frame1 = beadList.frame==1;
            fedPos{d} = double([beadList.xc(frame1),beadList.yc(frame1)]);
        end

        %% Align separate stains

        if showCellImages
            matchFig = figure(1); clf;
            alignFig = figure(2); clf;
        else
            matchFig = []; 
            alignFig = []; 
        end
        [tform,warpUncert] = AlignHybes(fedPos,...
            'maxD',maxD,...
            'warp2Hybe1',false,...
            'alignFig',alignFig,...
            'matchFig',matchFig,...
            'daxFile',beadImage,...
            'rawData',beads,...
            'showPlots',showCellImages); 
        if saveCellImages
            export_fig(alignFig,[savePath,'hybe_alignment_cell',cellnum,'.png']);               
            saveas(alignFig,[savePath,'hybe_alignment_cell',cellnum,'.fig']); 
        end
        %% Plot AlignDax
         
        if showCellImages; 
           snapshots = figure(3); close; figure(3); clf;  
           set(gcf,'Units','Inches','Position',[1,1,7.1,9.5],'PaperUnits', 'Inches', 'PaperSize', [7.1, 9.5]);
        else
            snapshots = [];
        end

        alignedIm = PlotAlignedDax(mRNAbin,listType,tform,...
                    'showPlots',showCellImages,...
                    'fighandle',snapshots);

        if showCellImages
            figure(3); FluorImage;
            alignedImFig = figure(4); clf; 
            Ncolor(uint16(alignedIm),jet(numHybes));
        end
            
        if saveCellImages
            export_fig(alignedImFig,[savePath,'alignedImFig_cell',cellnum,'.png']);    
            export_fig(snapshots,[savePath,'hybe_snapshots_cell',cellnum,'.png']);    
        end
        
        
        colordef white; % change colordef back to white for furture images
        
        %% Apply stain-alignment optional filters to STORM data

        % slow step -- lots of knn search calls.  
        if showCellImages; 
            clusterFig = figure(4); hold on;
        else
            clusterFig = []; 
        end

        [xf,yf,flists] = FilterConvData(imLists,tform,...
            'minPhotons',1,...
            'startframe',1,...
            'endframe',1,...
            'showPlots',showCellImages);

        
%%  Plot spots ontop of aligned dax image

        % plot dots on top of image;
        if showCellImages
                snapshots = figure(3); close; figure(3); clf;
            for d=1:numHybes
                subplot(4,numHybes/4,d); imagesc(alignedIm(:,:,d)); colormap gray; hold on; 
                plot(xf{d},yf{d},'.','color',clrmap(d,:));
            end
            set(gcf,'Units','Inches','Position',[1,1,7.1,9.5],'PaperUnits', 'Inches', 'PaperSize', [7.1, 9.5]);
            FluorImage;
        end     

        if saveCellImages
            export_fig(clusterFig,[savePath,'alignedImFigOverlay_cell',cellnum,'.png']);   
            export_fig(snapshots,[savePath,'hybe_snapshotsOverlay_cell',cellnum,'.png']);  
        end
        
        colordef white; % change colordef back to white for furture images

        %% Cluster localizations to find centroids for each mRNA
        if showCellImages
            histFig = figure(5); clf;
        else
            histFig = [];
        end;

        [mRNAcents,stainCents,stainCounts,clusterStats] = FindMRNA(xf,yf,...
                                        'binSize',binSize,...
                                        'minDotPerBin',minDotPerBin,...
                                        'minLocsPerDot',minLocsPerDot,...
                                        'minArea',minArea,...
                                        'maxArea',maxArea,...
                                        'showPlots',showCellImages,...
                                        'clusterFig',clusterFig,...
                                        'histFig',histFig);
        if saveCellImages
            export_fig(clusterFig,[savePath,'clusters_cell',cellnum,'.png']);               
            saveas(clusterFig,[savePath,'clusters_cell',cellnum,'.fig']);
        end

        %% Extract brightness from dax files     
        mRNApixels = sub2indFast([256,256],round(mRNAcents(:,2)),round(mRNAcents(:,1)));
        mRNApixels(mRNApixels<1) = 1;
        daxBrightness = zeros(numHybes,size(mRNApixels,1));
        for h=1:numHybes
            hybeIm = alignedIm(:,:,h);
            daxBrightness(h,:) = hybeIm(mRNApixels);
        end
        
        % figure(1); clf; imagesc(daxBrightness);
        
%         % Some test plots to ensure we grabbed the right pixels
%         % 
%         d = 1;
%         mRNApixels = sub2indFast([256,256],round(yf{d}),round(xf{d}));  % per hybe spots  
%         hybeIm = alignedIm(:,:,d)
%         daxBrightness = hybeIm(mRNApixels)
%         figure(1); clf; PlotCorr(flists{d}.a,daxBrightness);
%         xlabel('DaoSTORM "a"'); ylabel('dax file pixel value');
%         PresentationPlot();
%         
%         % validate that we chose the right pixels
%         hybeIm(mRNApixels) = 0; 
%         figure(2); imagesc(hybeIm); hold on;
%         plot(xf{d},yf{d},'go');

      
      
        %% Assign localizations from each stain to mRNA centroids
        % Method 2: map using just the filtered localizations straight
        [wordsDetected,brightnessPerSpot] = AssignConvSpotsToCentroids(flists,mRNAcents,...
            'maxDtoCentroid',maxDtoCentroid,...
            'minPhotsPerStain',minPhotsPerStain);
       

        %%  Analyze clusters

        [libExpect,libExpectPerHybe,libGenes,libCodes] = ...
                                    ExpectedLocsPerHybe(codebook,...
                                    'fpkmPath',fpkmPath,...
                                    'codebookPath',codebookFile);

        [perCnt,allCnt,spotMatrix,spotIDs,spotMatrixDax] = DecodeConvImage(wordsDetected,codebook,...
            'brightnessPerSpot',brightnessPerSpot,'daxValuePerSpot',daxBrightness');


        if showCellImages
            % Plot Counts per Gene
            cntsFig = figure(5); clf;
            subplot(2,1,1); BarGraphAllGenes(allCnt,libGenes,'codebook',codebook) ;        
            subplot(2,1,2); BarGraphAllGenes(perCnt,libGenes) ;
            set(gcf, 'Units','Inches','Position',[0,0,14,8],'PaperUnits', 'Inches', 'PaperSize', [14 8]);

            % Plot Distribution of ClusterCounts
            hybehistFig = figure(6); clf; 
            hist( sum(wordsDetected,2),0:numHybes); xlim([0,numHybes]);
            xlabel('number of hybes per cluster');
            PresentationPlot();

            % Plot Gene locations
            mRNALocFig = figure(7); clf;
            PlotConvGeneLocations(alignedIm,mRNAcents,libGenes,libCodes,spotIDs,'mlists',flists,'MarkerSize',1,'showNames',true);
            set(gcf, 'Units','Inches','Position',[0,0,12,10],'PaperUnits', 'Inches', 'PaperSize', [12 10]);

            % Plot FPKMs
            fpkmFig = figure(8); clf; 
            subplot(1,2,1); PlotCountsVsFPKM(libExpect,allCnt,libGenes);
            subplot(1,2,2); PlotCountsVsFPKM(libExpect,perCnt,libGenes);
            set(gcf, 'Units','Inches','Position',[0,0,14,5],'PaperUnits', 'Inches', 'PaperSize', [14 5]);
            expectHybeFig = figure(9); clf; PlotExpectHybe(libExpectPerHybe,hybeLocs);            
            
            % Plot pairwise correlations between bits
            bitCorrFig = figure(10); clf; 
            [corrM,bitPairs] =  BitCorrelation(wordsDetected);  
        end

        if saveCellImages
            export_fig(cntsFig,[savePath,'geneCounts_cell',cellnum,'.png']); 
            export_fig(hybehistFig,[savePath,'hybesPerSpot_cell',cellnum,'.png']);  
            export_fig(bitCorrFig,[savePath,'bitCorrelation_cell',cellnum,'.png']);
            export_fig(mRNALocFig,[savePath,'mRNALocations_cell',cellnum,'.png']);   
            saveas(mRNALocFig,[savePath,'mRNALocations_cell',cellnum,'.fig']);   
            export_fig(fpkmFig,[savePath,'fpkmVScounts_cell',cellnum,'.png']);
            saveas(fpkmFig,[savePath,'fpkmVScounts_cell',cellnum,'.fig']);
            export_fig(expectHybeFig,[savePath,'expectHybeCounts_cell',cellnum,'.png']);
        end

        rnaCnts(:,c) = allCnt;
        perCnts(:,c) = perCnt;
        spotData{c} = spotMatrix; 
        spotDaxData{c} = spotMatrixDax;


    %     wordsWithCorrBits = WordsWithCorrBits(bitPairs,libCodes,libGenes);
    %     [topNwords,topNcnts] = TopNWords(wordsDetected,libGenes,libCodes);
    %
    %     idxNegCntrls = 13:14; % indices of the negative control RNAs
    %     DistanceFromTopWords(idxNegCntrls,topNwords,topNcnts,libCodes,EcCount)
    catch er
        warning(['Problem analyzing data from cell ',num2str(c)]);
        warning(er.getReport);
        disp('skipping this cell');
    end
end

%% All Cell Stats

%% Save Data and script
save([savePath,'cntData.mat'],'rnaCnts','perCnts','spotData','spotDaxData','libGenes','libExpect','libCodes','minPhotsPerStain','maxDtoCentroid','codebook','hybeBrightness');      
copyfile( [mfilename('fullpath'),'.m'],[savePath,mfilename,'.m']);

%% All Cell Statistics

close all;
AllCellStats(savePath);
     

%%

