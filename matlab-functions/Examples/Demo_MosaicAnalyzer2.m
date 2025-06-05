%% A mimimal Mosaic Analyzer Script
% 04/17/2024  

% clc; clear all; close all; startup; 


%% filepaths 
% E10.5 dataset
NAS02_Vol1 = '\\169.254.113.81\NAS02_Vol1\';
expFolder = [NAS02_Vol1,'Sedona\2023-08-23_E105_RNA_DNA_replicate\'];
npp = 108;

%% 
rnaFolder = [expFolder, 'RNA_Expt\'];  
analysisFolder = [rnaFolder,'Analysis\'];
saveFolder = [analysisFolder,'RNA_Fits\'];
dna_analysis_folder = [expFolder,'DNA_Expt\Analysis\'];
allFitFiles = FindFiles([dna_analysis_folder,'*AllFits.csv']);
dna_bins = 66;
%% load data
rnaTableOutName = [saveFolder,'rnaSpots.csv'];
rnaTable = readtable(rnaTableOutName); % >9 million (!) RNA molecules loading takes a few minutes 
rnaNames = unique(rnaTable.name,'stable'); rnaNames = rnaNames(1:16)
nR = length(rnaNames);

[polys,maps,spotXYf] = CombineAllFits(dna_analysis_folder,'bins',dna_bins);
dna_aligned = readtable([dna_analysis_folder,'dna_aligned_XY.csv']);
%%

% table size
ymin = min(rnaTable.y)-500;
ymax = max(rnaTable.y)+500;
xmin = min(rnaTable.x)-500;
xmax = max(rnaTable.x)+500;

xy = [dna_aligned.x,dna_aligned.y];

 rs = 1:nR;
K = length(rs);
cmap = hsv(K);
nCells = size(xy,1);
rna_cnts = zeros(nCells,nR);
imRNA = cell(K,1);
figure(9); clf;
k = 0; 
 plot(dna_aligned.x,dna_aligned.y,'k+','MarkerSize',2); hold on;
prev_RNA = [];
for r=rs
    k=k+1;
    disp(rnaNames{r});
    % get current RNA channel from data table. 
    isRNA_r = strcmp(rnaTable.name,rnaNames{r}) ;
    rna_r = rnaTable(isRNA_r,:);

    % some ad-hoc filtering.  It's worth playing around with different
    % values and seeing which filters remove the cross-talk and ectopic
    % patterns better.
    sel = rna_r.xsigma < 1.6 & rna_r.xsigma > 0.2  &  rna_r.significance > 40 & rna_r.significance < 150;
    rna_r = rna_r(sel,:);

    % figure(1); subplot(4,4,r); hist(rna_r.significance,20:150); xlim([20,150]);

    % attempt to remove blead-through spots
    if ~isempty(prev_RNA)
        % [id,dis] = knnsearch([prev_RNA.x,prev_RNA.y,prev_RNA.frame],[rna_r.x,rna_r.y,rna_r.frame],'K',1); % (3D overlap)
        [id,dis] = knnsearch([prev_RNA.x,prev_RNA.y],[rna_r.x,rna_r.y],'K',1); % (2D overlap)
        drop = dis<1e3/npp;      
        % figure(10); clf;
        % figure(10); plot(rna_r.x,rna_r.y,'rx'); hold on;
        % figure(10); plot(prev_RNA.x,prev_RNA.y,'b+'); hold on;
        % figure(10); plot(rna_r.x(drop),rna_r.y(drop),'ko'); hold on;
        % disp('test')
        rna_r = rna_r(~drop,:);
    end
    if height(rna_r) < 1
        disp('no RNA left')
    end
    figure(9); plot(rna_r.x,rna_r.y,'.','color',cmap(k,:),'MarkerSize',2); hold on;
    title(rnaNames{r}); pause(.1);


    % count RNA per DNA spot within a 10um radius
    if ~isempty(rna_r)
        [id,dis] = knnsearch([rna_r.x,rna_r.y],xy,'K',500);
        rna_cnts(:,r) = sum(dis<10e3/npp,2); % convert nm back to pixels
    end

    % bin the RNA spots to create new RNA-images
    % divide embryo into ~10x10 um boxes and count RNA per box (100 108 nm pixels)    
      imRNA{r} =  hist3([rna_r.y,rna_r.x],'Ctrs',{ymin:100:ymax,xmin:100:xmax});

    figure(11); clf; 
    imRNAc = Ncolor( IncreaseContrast( uint16(cat(3,imRNA{r})),'high',.9995) );
    imagesc(imRNAc);
    figure(11); clf;
    plot(rna_r.x,rna_r.y,'.','color',cmap(k,:),'MarkerSize',2); hold on;
    title(rnaNames{r}); pause(.1);

      % record previous RNA to remove bleadthrough
      prev_RNA =   cat(1,prev_RNA,rna_r);
end 
figure(9); legend('dna',rnaNames{rs});
box off; set(gcf,'color','w')
title('RNA + DNA mosaic')

%% Show the new RNA images from binning the RNA localizations 
figure(10); clf;
rs =1:5;
imRNAc = Ncolor( IncreaseContrast( uint16(cat(3,imRNA{rs})),'high',.9995) );  
subplot(1,4,1); imagesc(imRNAc); set(gca,'YDir','normal');
colormap(hsv(length(rs)))
cb2 = colorbar; set(cb2,'YTick',0.5/length(rs) + linspace(0,1,length(rs)+1),'YTickLabel',rnaNames(rs));

rs =6:10;
subplot(1,4,2); 
imRNAc = Ncolor( IncreaseContrast( uint16(cat(3,imRNA{rs})),'high',.9995) );  
 imagesc(imRNAc); set(gca,'YDir','normal');
colormap(hsv(length(rs)))
cb2 = colorbar; set(cb2,'YTick',0.5/length(rs) + linspace(0,1,length(rs)+1),'YTickLabel',rnaNames(rs));


rs =[10:14];
subplot(1,4,3); 
imRNAc = Ncolor( IncreaseContrast( uint16(cat(3,imRNA{rs})),'high',.9995) );  
 imagesc(imRNAc); set(gca,'YDir','normal');
colormap(hsv(length(rs)))
cb2 = colorbar; set(cb2,'YTick',0.5/length(rs) + linspace(0,1,length(rs)+1),'YTickLabel',rnaNames(rs));

rs =[14:16,16,16];
subplot(1,4,4); 
imRNAc = Ncolor( IncreaseContrast( uint16(cat(3,imRNA{rs})),'high',.9995) );  
 imagesc(imRNAc); set(gca,'YDir','normal');
colormap(hsv(length(rs)))
cb2 = colorbar; set(cb2,'YTick',0.5/length(rs) + linspace(0,1,length(rs)+1),'YTickLabel',rnaNames(rs));

%% convert RNA binned image into cell x gene table
RNA_image = cat(3,imRNA{:});
[nX,nY,nG] = size(RNA_image);

isEmpty = sum(RNA_image,3)==0; % find empty cells
[Y,X] = meshgrid(1:nY,1:nX); % all bins (including ones with no RNA)
x = X(~isEmpty); % keep only the non-empty ones
y = Y(~isEmpty);

nCells = sum(~isEmpty(:));
RNA_cellXgene = nan(nCells,nG);
for g=1:nG % loop over genes
    curr_RNA_im = RNA_image(:,:,g);
    RNA_cellXgene(:,g) = curr_RNA_im(~isEmpty);
end


figure(1); clf; imagesc(RNA_cellXgene)



%% Hierarchical clustering of RNA counts 
sel = randi(nCells,1000,1); % select a small subset of genes for speed
% (warning a bit slow)
rna_norm = RNA_cellXgene(sel,:)'; % change sel to ':' to process all genes
[nGenes,nCells] = size(rna_norm)
rna_norm = (rna_norm - repmat(mean(rna_norm,2),1,nCells))./ repmat(std(rna_norm,[],2),1,nCells);
figure(1); clf; imagesc((rna_norm)); colorbar; clim([-3 3]);
set(gca,'YTick',1:nGenes,'YTickLabel',rnaNames);
cgo = clustergram(rna_norm);
cgo.Colormap = colormap('default');
cgo.RowLabels = rnaNames;
cgo.DisplayRange = 3;

%% UMAP
% 
umap_path = 'C:\Users\Alistair\Documents\code\external\mlab_umap\';
addpath(genpath(umap_path))

hoxN = rna_cnts;
U = run_umap(hoxN);
%% examine UMAP results
% plot some RNA labels on top of the UMAP
figure(3); clf; plot(U(:,1),U(:,2),'k.','MarkerSize',1);
cmap1 = hsv(5);
for r=1:5
    isOn = hoxN(:,r) > quantile(hoxN(:,r),.8);
    figure(3); hold on; plot(U(isOn,1),U(isOn,2),'o','color',cmap1(r,:),'MarkerSize',7-r);
end
title('umap')
legend('data',rnaNames{:});

% see how different RNA levels correlate with the UMAP vectors
figure(4); clf; 
subplot(1,2,1); PlotCorr(U(:,1),hoxN(:,4),'hex',true,'log',false); xlabel('U1'); ylabel(rnaNames{4})
subplot(1,2,2); PlotCorr(U(:,2),hoxN(:,9),'hex',true,'log',false);  xlabel('U1'); ylabel(rnaNames{9})

%% Using the manual draw tool to select regions based on the 3D map
%% Intialize region selector 
N = 3; % decide how many different regions you plan to draw
idx = {};
cmap = hsv(N);
k=0;
%% draw up to N polygons
% rerun this block up to N ten times.  Each time will increase the counter
% (k), and store a new set of indices (idx)
k=k+1;
ax = gca;
roi = drawpolygon(ax) ;  %  creates the ROI in the axes specified by ax.
polyXY = roi.Position;
idx{k} = inpolygon(xy(:,1),xy(:,2),polyXY(:,1),polyXY(:,2));
plot(xy(idx{k},1),xy(idx{k},2),'o','color',cmap(k,:));
% could easily extend this to grab RNA as well, and we could count the
% total and relative abundance of all the different RNA species in the area
% selected

%% compute DNA contact maps from these N polygons areas
nG = length(idx);
mapG = cell(length(idx),1);
for k=1:length(idx)
mapG{k} = maps(:,:,idx{k});
end

cMap =ContactFrac(maps,'threshold',150);

cmaps = cell(nG,1);
figure(1); clf; figure(2); clf; figure(3); clf;
for g=1:nG
    figure(1); subplot(1,nG,g); imagesc(nanmedian(mapG{g},3)); caxis([200,550]); 
    [cmaps{g},omap] = ContactFrac(mapG{g},'threshold',150);
    figure(3); subplot(1,nG,g); imagesc(cmaps{g}); caxis([0.05, .5]); 
end


figure(2); clf;
for g=1:nG
figure(2); subplot(1,nG,g); imagesc(cmaps{g} - cMap); caxis([-.5, .5]); 
end
GetColorMap('RedWhiteBlueK'); colorbar;

cMapAxis = nanmedian(cat(3,cmaps{[1:8,10]}),3);
for g=1:nG
figure(2); subplot(1,nG,g); imagesc(cmaps{g} - cMapAxis); caxis([-.5, .5]); 
end
GetColorMap('RedWhiteBlueK'); colorbar;
