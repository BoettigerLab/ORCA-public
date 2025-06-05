% clc; clear all; startup;


expFolder = '\\169.254.113.81\NAS02_Vol1\Sedona\2023-08-16_96-168hrgastruloids_hoxRNA_try2\';
rnaFolder = [expFolder, 'RNA_Expt\']; 
analysisFolder = [rnaFolder,'Analysis\']; %  
saveFolder = [analysisFolder,'RNA_Fits\'];

% Load RNA names
% (Yes, you should make an experiment table for every dataset).
rnaTable = readtable([rnaFolder,'Hyb_Identities.xlsx'])  
rnaNames = rnaTable.RNA;


% hybFolders = FindFiles([rnaFolder,'Hyb*'],'onlyFolders',true); % autofind 
hybFolders = strcat(rnaFolder,rnaTable.FolderName,'/')  % in this case Hybs 1-5 are membrane barcodes, not RNA

 % only the first 16 are RNA data, the Hybs in the RNA experiment are
 % membrane barcodes and WGA / DAPI. We will address these separately.
hybFolders = hybFolders(1:16)
rnaNames = hybFolders(1:16);
nH = length(hybFolders);

% Find the dax files 
daxFiles = FindFiles([hybFolders{1},'*.dax']); 
dc = (1:nH);  
nFOV = length(daxFiles);
fov_nums = cellstr( num2str( (1:nFOV)' ));

% mosaic tile orientation 
%   (the different scopes have cameras oriented different ways and on
%   different sides of the scopebody, which leads to a scope specific
%   orientation).  
flipV = false; flipH = true; transposeI = true; % scope3 specific transpose

%% Quick view mosaic
% verify that the mosaic assembled correctly
[imTiles,uls,miniMosaics] = DaxToTiles(daxFiles,'verbose',false);   % a compact function 
figure(1); clf; imagesc(miniMosaics{1})
%  mosaicImage =  TilesToMosaic(imTiles(:,c),uls); % helper to convert the imTiles to a mosaic   
% 
[imH,imW] = size(imTiles{1});
fovs =1:nFOV;


%% compute overlap (simple and complete)
% we use this to avoid double-counting the RNA
% (to update - should do this after drift correction). 
boxOverlap = nan(nFOV,4,nFOV);
for f=2:nFOV
    for ff=1:f-1
        boxOverlap(f,:,ff) = [uls(ff,1:2),imW,imH];
    end
end


%% RNA drift correction
% by default this will 
% recomp
driftData = FindFiles([analysisFolder,'alignTable_fov*']);
if length(driftData) < nFOV
    SetFigureSavePath(analysisFolder);
    datImageAlign = cell(nFOV,nH);
    for f=1:nFOV
        [~,daxName]=fileparts(daxFiles{f}); %  (dax names don't change between hybes)  
        loadFiles = strcat(hybFolders,daxName,'.dax'); %  folder names do change 
        [~,xyMax] = ShowDaxSeries(loadFiles,'driftCorrectUsingChn',2,'analysisFolder',analysisFolder,'overwriteDrift',false); 
        for h=1:nH
            datImageAlign{f,h} = xyMax(:,:,h);
        end
    end
end

%% plot all RNA
overwrite = false; 

selFOV = 1:nFOV;
driftCorrectFolder = analysisFolder;



figure(1); clf;
hs = 1:nH   % [3,4,6,10,12] %
rnaNames(hs)
cmap = hsv(length(hs));
rnaSpots = cell(nH,1);
for f=selFOV
    k = 0;
    driftCorrectTable = readtable([driftCorrectFolder, 'alignTable_fov',num2str(f-1, '%03d'), '.csv']);
    if f==1
        driftCorrectTable  % display table for validation
    end

    % loop over RNA hybes
    for h=hs
        isH = strcmp(rnaTable.FolderName{h},driftCorrectTable.hyb);
        fidAlign = driftCorrectTable(isH,:);
     
        k=k+1;
        % % earlier name convention
        % daxNewName = ['3D_Hyb',num2str(h,'%02d'),'_FOV',num2str(f,'%02d')]; 

        % new name convention
         daxNewName = ['3D_Hyb',rnaNames{h,1},'_FOV',num2str(f,'%02d')]; % full 3D
        binName = [saveFolder,daxNewName,'.hdf5'];
        try
            data = LoadDaoFits(binName,'verbose',false);
        catch
            data = [];
        end
        if isempty(data)
            continue
        end
        x = data.x;
        y = data.y;
        x = x + fidAlign.xshift + fidAlign.xshift2;
        y = y + fidAlign.yshift + fidAlign.yshift2;
        if flipV
            y = imH-y;
        end
        if flipH
            x = imW-x;
        end
        if transposeI
            xx =x;
            yy =y;
            x = yy;
            y = xx;
        end 
        x = x + uls(f,1); 
        y= y + uls(f,2);
        % remove overapped points
        % boxOverlap  nFOV x 4 x nFOV
        inBox = false(length(x),1);
        for ff=1:f-1
            x1 = boxOverlap(f,1,ff);
            x2 = boxOverlap(f,1,ff)+boxOverlap(f,3,ff);
            y1 = boxOverlap(f,2,ff);
            y2 = boxOverlap(f,2,ff)+boxOverlap(f,4,ff);
            inBox = inBox | (x > x1 & x < x2 & y > y1 & y < y2);
        end

        % --- filter data
        % sel = data.significance > quantile(data.significance,.2); 
        sel = data.xsigma < 1.9 & ~inBox  & data.significance > 40 & data.significance < 150;
        rna_x = x(sel); rna_y = y(sel);
        figure(1); plot(rna_x,rna_y,'.','color',cmap(k,:),'MarkerSize',1); hold on;
        text(uls(f,1)+imW/2,uls(f,2)+imH/2,num2str(f-1),'Color','k','FontSize',12);


        rna_xy = data(sel,:); % [rna_x,rna_y];
        rna_xy.name = repmat(rnaNames(h),height(rna_xy),1);
        rna_xy.x = rna_x; % in stage coordinates with drift correction
        rna_xy.y = rna_y; %  in stage coordinates with drift correction
        try
            rna_xy.category = [];  % remove unnecessary lines from table
        catch
        end
        rnaSpots{h} = cat(1,rnaSpots{h},rna_xy);      
    end
end

cb = colorbar; colormap(cmap);
set(cb,'YTick',0.5/length(hs) + linspace(0,1,length(hs)+1),'YTickLabel',rnaNames(hs));
set(gcf,'color','w');


rnaSpotTableName = [analysisFolder,'rnaSpots.csv'];
rnaSpotTable = cat(1,rnaSpots{:});

%% Write RNA table
if overwrite || ~exist(rnaSpotTableName,'file')
    writetable(rnaSpotTable,rnaSpotTableName);
    disp(['wrote ',rnaSpotTableName, '  with n spots=',num2str(height(rnaSpotTable))])
end


%% load and align DNA
% DaxToTiles is a bit slow the first time you run it (it has to load all
% the data to create the max-projection files). 

npp = 108; % scope specific xy pixel size
dna_analysis_folder = [expFolder,'DNA_Expt\Analysis\'];
dnaHybs = FindFiles([expFolder,'DNA_Expt\Hyb*'],'onlyFolders',true);
daxFilesDNA = FindFiles([dnaHybs{1},'Conv*.dax']);
nFOV_dna = length(daxFilesDNA);
nBins = length(dnaHybs);
[imTilesDNA,ulsDNA,miniMosaics] = DaxToTiles(daxFilesDNA); 

allFitFiles = FindFiles([dna_analysis_folder,'*AllFits.csv']);

%% compute overlap from DNA
[imW,imH] = size(imTilesDNA{1});
boxOverlap_DNA = nan(nFOV_dna,4,nFOV_dna);
for f=2:nFOV_dna
    for ff=1:f-1
        boxOverlap_DNA(f,:,ff) = [ulsDNA(ff,1:2),imW,imH];
    end
end

%%
figure(7); clf;
npp = 108;
showPlot = true;
nBins = length(dnaHybs);
allMaps = cell(nFOV_dna,1);
dnaSpotTables = cell(nFOV_dna,1);
for f=1:nFOV_dna
[polys,maps,spotXY] = TableToPolymer(allFitFiles{f},'bins',nBins); % load dna data
    dna_x =  spotXY(:,1)/npp;
    dna_y =  spotXY(:,2)/npp;
    if flipV
        dna_y = imH-dna_y;
    end
    if flipH
        dna_x = imW-dna_x;
    end
    if transposeI
        xx =dna_x;
        yy =dna_y;
        dna_x = yy;
        dna_y = xx;
    end
    dna_x = dna_x + ulsDNA(f,1); 
    dna_y= dna_y + ulsDNA(f,2);
        
    rna_cnts = zeros(length(dna_x),nH);
    if showPlot
        figure(7);
        plot(dna_x,dna_y,'r+'); hold on;
        text(ulsDNA(f,1),ulsDNA(f,2),num2str(f),'color','k','FontSize',12);
    end

    dnaSpotTables{f} = [dna_x,dna_y]; %
    allMaps{f} = maps;
    if showPlot
        figure(6); clf; imagesc(nanmedian(maps,3)); 
        colorbar; caxis([100,500]);
        title(['FOV ',num2str(f)])
        pause(.01);
    end
end

%%
allMap = cat(3,allMaps{:});
[cf,nObs] = ContactFrac(allMap,'threshold',250);
xy_dna = cat(1,dnaSpotTables{:});

%% Align RNA to DNA -- step 1 convert the points to images by binning 
sc = 50; % scale to map data
xy_rna = [rnaSpotTable.x,rnaSpotTable.y]; % [rnaSpots{h}.x,rnaSpots{h}.y];
rna_lo = min(xy_rna);
xy_rna = [xy_rna(:,1)-rna_lo(1),xy_rna(:,2)-rna_lo(2)]; % make  coordinates positive  

dna_lo = min(xy_dna);
xy_dna = [xy_dna(:,1)-dna_lo(1),xy_dna(:,2)-dna_lo(2)]; % make  coordinates positive

rna_hi = max(xy_rna);
dna_hi = max(xy_dna);
im_wh = max([rna_hi; dna_hi],[],1)
xs = 0:sc:im_wh(1);
ys = 0:sc:im_wh(2);

imRNA = hist3([xy_rna(:,2),xy_rna(:,1)],'Ctrs',{ys,xs});
imDNA = hist3([xy_dna(:,2),xy_dna(:,1)],'Ctrs',{ys,xs});

% % optional (with multiple embryos on the same slide with large gaps
% inbetwen, in may be better to align them independently)
% im_wh(1) = 12000; % truncate
% imRNA = imRNA(:,1:(12000/sc));
% imDNA = imDNA(:,1:(12000/sc));

figure(6); clf; 
subplot(2,1,1); imagesc(imRNA); clim([0,5]);
subplot(2,1,2); imagesc(imDNA); clim([0,5]);

imRNA_1 = imRNA > 5; % cleanup some noise in the RNA channel
imDNA_1 = imDNA > 0;

%% Align DNA-to-RNA step 2, compute registration -- Key step
% This will go faster if you have a narrower range of angles to scan. You
% might start with a something coarse (-20:2:20) and refine to something
% more precise (e.g. -4:.25:-2) after you have an estimate of what rotation
% is needed. 
figure(30); clf;
 [alignO,imO] = CorrAlignFast(imRNA_1,imDNA_1,'angles',-15:.5:15,'scales',1); 


 %% Apply the alignment
 alignA = alignO;
 xy_dna_shift =xy_dna/sc; 
ptsOut = RotateTranslatePoints(xy_dna_shift,alignA,...
    'center',(im_wh/2)/sc,'upperLeft',false,'invert',false); % 
ptsOut = ptsOut*sc;

inBlock = xy_rna(:,1) < inf; %  12000;
figure(5); clf; 
subplot(1,2,1);
plot(xy_rna(inBlock,1),xy_rna(inBlock,2),'k.'); hold on;
plot(xy_dna(:,1),xy_dna(:,2),'ro'); hold on;
set(gca,'YDir','reverse')
subplot(1,2,2);
plot(xy_rna(inBlock,1),xy_rna(inBlock,2),'k.'); hold on;
plot(ptsOut(:,1),ptsOut(:,2),'ro'); hold on;
set(gca,'YDir','reverse')

x = ptsOut(:,1);
y = ptsOut(:,2);
dna_aligned = table(x,y);

dna_aligned_name = [dna_analysis_folder,'dna_aligned_XY.csv'];
if overwrite || ~exist(dna_aligned_name,'file')
    writetable(dna_aligned,dna_aligned_name);
    disp(['wrote ',dna_aligned_name]);
end

%% 
% MosaicAnalyzer2 uses these dna_aligned table and RNAspots table
