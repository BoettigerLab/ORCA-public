
% TSTORM
% Analyze multi-day STORM images of sequential staining


%% Parsing filepaths:
% startup;
codebookFile = 'C:\Users\Alistair\Documents\Research\Projects\TSTORM\Lib\Lib5_140615\E1_codebook.fasta';
dataPath = '\\cajal\TSTORMdata\140709_L5E1\';
fpkmPath = 'D:\Data\Genomics\HumanGenome\IMR90\mthread_IMR90_totRNA_rep1\';
listType = 'mlist.bin';

savePath = [dataPath,'Analysis1\'];
if ~exist(savePath,'file')
    mkdir(savePath);
end

hybe1Data = dir([dataPath,'STORM_1','_*_','c1_',listType]);
codebook = fastaread(codebookFile);

numGenes = length(codebook);  % number of genes (including neg controls) 
numHybes = length(regexprep(codebook(1).Header,' ','')); % number of stains 
numCells = length(hybe1Data);

if numCells < 1
    error(['No cells found in ',dataPath,' of type ',listType]);
end

rnaCnts = zeros(numGenes,numCells);

%% Filters

%-----------------------------------------------------------------------
% choose a colormap 

clrmap = jet(numHybes);

% Drift Correction Parameters:
xshift = zeros(1,numHybes);
yshift = zeros(1,numHybes);
maxdrift = 8; % pixels
showDriftPlots = true; % 
minFramesForFeducials = .5; % 
samplingrate = 60; 
startframe = 1; % for beads;

% Parameters for aligning hybes
maxD = 10; % max distance over which to attempt to match beads
warp2Hybe1 = false; 

% Parameters for Filtering STORM data
minPhotons =   1000; % updated  ;  % min # photons per spot
maxDistanceForDensityFilter = .5; % max distance for density filter (in pixels)
minDensity = 3; % min number of localizations within maxD
STORMstartingFrame = 1; 

% Parameters for Clustering Localizations
binSize = 12;  % size of bins in nm
minDotPerBin = 2;  % min number of localizations to call a bin occupied
minLocsPerDot = 80; % min number of localization in all bins assigned to a cluster to be called an mRNA
minArea = 0;   % min area in bins to be called a cluster of localization
maxArea = 300; % max area in bins to be called a cluster of localization

% Parameters for Assigning stains to mRNA centroids
minLocsPerStain = 10; % min number of localizations within MaxDtoCentroid to be called 'detected in stain d'
maxDtoCentroid = .4; % max distance you can be from the centroid of the mRNA

for c = 1:numCells; % c=2

    cellnum = sprintf('%01d',c-1);

    mRNAdata = dir([dataPath,'*STORM*','_',cellnum,'_','c1_',listType]);
    mRNAbin = strcat(dataPath,{mRNAdata.name}');
    
    beaddata = dir([dataPath,'*STORM*',cellnum,'_','c2_list.bin']);
    beads = strcat(dataPath,{beaddata.name}');
    
    % Remove bit 0 from list of hybes 
    bit0 = StringFind(mRNAbin,'STORM_0_');
    mRNAbin(bit0) = []; 
    beads(bit0) = []; 

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
    % slow step, potentionall rough on the memory.  

    imLists = cell(numHybes,1); 
    hybeLocs = zeros(numHybes,1); 
    hybeBrightness = zeros(numHybes,1); 
    for d=1:numHybes;
        imLists{d}   = ReadMasterMoleculeList(mRNAbin{d}); 
        hybeLocs(d) = length(imLists{d}.xc);
        hybeBrightness(d) = median(imLists{d}.a);
    end


    %% Use beads to compute drift correction in each movie
    driftFig = figure(1); 
    [dxc,dyc,fedPos] = CorrectImDrift2(beads,imLists,... % Compute drift
                        'samplingrate',samplingrate,...
                        'integrateframes',samplingrate,...
                        'maxdrift',maxdrift,...
                        'startframe',startframe,...
                        'fmin',minFramesForFeducials,...
                        'fighandle',driftFig);
    saveas(driftFig,[savePath,'drift_cell',cellnum,'.png']);               
    saveas(driftFig,[savePath,'drift_cell',cellnum,'.fig']);

    %% Align separate stains

    matchFig = figure(2); clf;
    alignFig = figure(3); clf;
    [tform,warpUncert] = AlignHybes(fedPos,'maxD',maxD,'warp2Hybe1',false,...
        'alignFig',alignFig,'matchFig',matchFig,'daxFile',beadImage,...
        'rawData',beads); 
    saveas(alignFig,[savePath,'hybe_alignment_cell',cellnum,'.png']);               
    saveas(alignFig,[savePath,'hybe_alignment_cell',cellnum,'.fig']); 

    %% Align frame 1 images
    
    figure(4); clf; 
    snapshots = figure(5); clf; 
    rawIm  = zeros(256,256,numHybes);
    alignedIm = zeros(256,256,numHybes);
    for d=1:numHybes
        daxfile = regexprep(mRNAbin{d},['_',listType],'.dax');
        dax = max(ReadDax(daxfile,'endFrame',20),[],3); % max project the first 50 frames 
        [H,W] = size(dax); 
        tformInv = fliptform(tform{d}); 
        alignedDax = imtransform(dax,tformInv,...
                        'XYScale',1,'XData',[1 W],'YData',[1 H]);
        rawIm(:,:,d) = dax;
        alignedIm(:,:,d) = alignedDax;
        figure(5); subplot(4,2,d); imagesc(uint16(alignedDax)); colormap gray;  % caxis([0,2^13.5]);
    end
    
    saveas(snapshots,[savePath,'hybe_snapshots_cell',cellnum,'.png']);    
    %% Apply drift-correction, stain-alignment, and optional filters to STORM data

    % slow step -- lots of knn search calls.  
    clusterFig = figure(4); clf; 
    [xf,yf] = FilterSTORMdata(imLists,dxc,dyc,tform,xshift,yshift,...
        'minPhotons',minPhotons,...
        'maxDistance',maxDistanceForDensityFilter,...
        'minDensity',minDensity,...
        'startframe',STORMstartingFrame,...
        'showPlots',true);


    %% Cluster localizations to find centroids for each mRNA
    histFig = figure(15); clf;
    [mRNAcents,stainCents,stainCounts,clusterStats] = FindMRNA(xf,yf,...
                                    'binSize',binSize,...
                                    'minDotPerBin',minDotPerBin,...
                                    'minLocsPerDot',minLocsPerDot,...
                                    'minArea',minArea,...
                                    'maxArea',maxArea,...
                                    'showPlots',true,...
                                    'clusterFig',clusterFig,...
                                    'histFig',histFig);
    saveas(clusterFig,[savePath,'clusters_cell',cellnum,'.png']);               
    saveas(clusterFig,[savePath,'clusters_cell',cellnum,'.fig']);
    
   
    x = 0:10:5E3;
    cnts = clusterStats.localizations(clusterStats.goodClusters);
    histCnts = figure(6); clf; hist(cnts,x); 
    ylabel('number of clusters');
    xlabel('number of localizations in cluster'); 
    saveas(histCnts,[savePath,'histCnts_cell',cellnum,'.png']);               
    saveas(histCnts,[savePath,'histCnts_cell',cellnum,'.fig']);
    
   
    
    %% Assign localizations from each stain to mRNA centroids

    % Method 1: map using clusters from individual stains
    % [Mcents,datPos]= AssignLocsToMRNAcentroids(stainCents,stainCounts,mRNAcents,...
    %                                       'minLocsPerStain',minLocsPerStain,...
    %                                       'maxDtoCentroid',maxDtoCentroid,...
    %                                       'method','clusters');

    % Method 2: map using just the filtered localizations straight
    [Mcents,datPos]= AssignLocsToMRNAcentroids(xf,yf,mRNAcents,...
                                          'minLocsPerStain',minLocsPerStain,...
                                          'maxDtoCentroid',maxDtoCentroid,...
                                          'method','dots');


    %% Use mRNA to compute new warp
    % This doesn't help more than a couple nm for most localizations


    % mRNAwarp = ComputeWarpByMRNA(mRNAcents,datPos,xf,yf,...
    %     'maxDtoCentroid',maxDtoCentroid);
    % 
    % % Apply warp to correct the localizations
    % for d=1:length(xf)
    % [xf{d},yf{d}]  = tforminv(mRNAwarp{d},xf{d},yf{d});
    % end
    % 
    % % Recompute clusters
    % [mRNAcents,stainCents,stainCounts] = FindMRNA(xf,yf,'binSize',binSize,...
    %                                 'minDotPerBin',minDotPerBin,...
    %                                 'minLocsPerDot',minLocsPerDot,...
    %                                 'minArea',minArea,...
    %                                 'maxArea',maxArea);
    %                             
    % % Re-assign localizations from each stain to a cluster                            
    % [Mcents,datPos]= AssignLocsToMRNAcentroids(xf,yf,mRNAcents,...
    %                                       'minLocsPerStain',minLocsPerStain,...
    %                                       'maxDtoCentroid',maxDtoCentroid,...
    %                                       'method','dots');


    %%  Analyze clusters
    
    [libExpect,libExpectPerHybe,libGenes,libCodes] = ExpectedLocsPerHybe(codebook,'fpkmPath',fpkmPath,'codebookPath',codebookFile);
    [uniqueMsg,hammingDis,freqUniqueMsg] = DecodeImage(Mcents,codebook);

    % Plot Counts per gene with and without error correction 
    cntsFig = figure(7); clf;
    [PmCount,EcCount] = PlotCountsPerGene(codebook,hammingDis,freqUniqueMsg,'libGenes',libGenes);
    saveas(cntsFig,[savePath,'geneCounts_cell',cellnum,'.png']);               


    % Plot Distribution of ClusterCounts
    ClusCount = zeros(size(uniqueMsg,2),1);
    for n = 1: size(uniqueMsg,2)
        ClusCount(n) = sum(freqUniqueMsg(sum(uniqueMsg,2)==n));
    end   
    hybehistFig = figure(8); clf; bar(ClusCount); 
    xlabel('number of hybes per cluster');
    PresentationPlot();
    saveas(hybehistFig,[savePath,'hybesPerSpot_cell',cellnum,'.png']);               


    % Compute which clusters are true mRNA and Plot gene Locations
    mRNALocFig = figure(clusterFig); clf; 
    showbkd = 0; 
    [correctClusters,spotMatrix] = PlotGeneLocations(mRNAcents,Mcents,uniqueMsg,hammingDis,codebook,xf,yf,libGenes,showbkd);
    nums = cellfun(@num2str,num2cell(1:numHybes),'UniformOutput',false);
   %  legend({nums{:},nums{:},'clusters',codebook.Sequence})
    saveas(mRNALocFig,[savePath,'mRNALocations_cell',cellnum,'.png']);   
    saveas(mRNALocFig,[savePath,'mRNALocations_cell',cellnum,'.fig']);   

     fpkmFig = figure(9); clf; PlotCountsVsFPKM(libExpect,EcCount)
     saveas(fpkmFig,[savePath,'fpkmVScounts_cell',cellnum,'.png']);
    % 
     expectHybeFig = figure(10); clf; PlotExpectHybe(libExpectPerHybe,hybeLocs);
     saveas(expectHybeFig,[savePath,'expectHybeCounts_cell',cellnum,'.png']);

    rnaCnts(:,c) = EcCount;

    % Plot pairwise correlations between bits
    bitCorrFig = figure(11); clf; [corrM,bitPairs] =  BitCorrelation(Mcents);
    saveas(bitCorrFig,[savePath,'bitCorrelation_cell',cellnum,'.png']);

    % 
     wordsWithCorrBits = WordsWithCorrBits(bitPairs,libCodes,libGenes);


    % 
    [topNwords,topNcnts] = TopNWords(Mcents,libGenes,libCodes);

%     %
%     idxNegCntrls = 13:14; % indices of the negative control RNAs
%     DistanceFromTopWords(idxNegCntrls,topNwords,topNcnts,libCodes,EcCount)

end

cellcorr = figure(12); clf;
subplot(1,2,1); PlotCorr(rnaCnts(:,1),rnaCnts(:,2)); xlabel('cell 1'); ylabel('cell 2');
subplot(1,2,2); PlotCorr(rnaCnts(:,3),rnaCnts(:,4)); xlabel('cell 3'); ylabel('cell 4');
PresentationPlot('MarkerWidth',20);
saveas(cellcorr,[savePath,'cellcell_correlation','.png']);

allCnts = sum(rnaCnts,2); disp(allCnts);
fpkmFigAll = figure(13); clf; PlotCountsVsFPKM(libExpect,allCnts);
     saveas(fpkmFig,[savePath,'fpkmVScounts_AllCells','.png']);

copyfile( [mfilename('fullpath'),'.m'],[savePath,mfilename,'.m']);
     