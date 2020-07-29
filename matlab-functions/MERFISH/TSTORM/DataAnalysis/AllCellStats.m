function AllCellStats(savePath,varargin)


load([savePath,'cntData.mat'],'rnaCnts','perCnts','spotData','spotDaxData','libGenes','libExpect','libCodes','minPhotsPerStain','maxDtoCentroid','codebook','hybeBrightness');      

numGenes = size(spotData{1},1);

colordef white;

% Cell-Cell correlation Plots
cellcorr = figure(12); clf;
subplot(1,2,1); PlotCorr(rnaCnts(:,1),rnaCnts(:,2)); xlabel('cell 1'); ylabel('cell 2');
subplot(1,2,2); PlotCorr(rnaCnts(:,3),rnaCnts(:,4)); xlabel('cell 3'); ylabel('cell 4');
set(gcf, 'Units','Inches','Position',[0,0,14,5],'PaperUnits', 'Inches', 'PaperSize', [14 5]);
PresentationPlot('MarkerWidth',20);
export_fig(cellcorr,[savePath,'cellcell_correlation','.png']);

% FPKM correlation Plots
allCnts = sum(rnaCnts,2); 
allCntData = flipud(sortrows([libExpect,allCnts,(1:numGenes)']));
disp(allCntData(:,2))
fpkmFigAll = figure(13); clf; PlotCountsVsFPKM(libExpect,allCnts,libGenes);
ylabel('all counts');
export_fig(fpkmFigAll,[savePath,'fpkmVScounts_AllCells','.png']);

allPerCnts = sum(perCnts,2);
fpkmPerAll = figure(14); clf; PlotCountsVsFPKM(libExpect,allPerCnts,libGenes);
ylabel('perfect match counts');
export_fig(fpkmPerAll,[savePath,'fpkmVsPerCounts_AllCells','.png']);     

% Bar Graphs counts per gene
perCntsFig = figure(16); clf; BarGraphAllGenes(perCnts,libGenes) ;
title(['Perfect matches per gene cluster. ', '(',num2str(sum(sum(perCnts))),' genes)']); 
set(gcf, 'Units','Inches','Position',[0,0,14,5],'PaperUnits', 'Inches', 'PaperSize', [14 5]);
export_fig(perCntsFig,[savePath,'barPerCnts','.png']); 
saveas(perCntsFig,[savePath,'perCntsFig','.fig']); 

corrCntsFig = figure(15); clf; BarGraphAllGenes(rnaCnts,libGenes,'codebook',codebook) ;
set(gcf, 'Units','Inches','Position',[0,0,14,5],'PaperUnits', 'Inches', 'PaperSize', [14 5]);
export_fig(corrCntsFig,[savePath,'barAllCnts','.png']);
saveas(corrCntsFig,[savePath,'barAllCnts','.fig']); 




%% Per Bit Analysis

colordef white; 

% figure(1); clf; bar(sum(perCnts,2)./sum(rnaCnts,2));

[numGenes,numHybes] = size(libCodes); 

allGenes = 1:numGenes;
negIdx = [StringFind(libGenes,'blank');StringFind(libGenes,'notarget')];
realGenes = allGenes; 
realGenes(negIdx) = [];
cntrlGenes = allGenes(negIdx); 
numRealGenes = length(realGenes);
numControls = length(cntrlGenes); 



% Combine spot matrices for each cell so we have one per gene
spotMatrices = cell(numGenes,1);
spotDaxMatrices = cell(numGenes,1); 
for n=1:numGenes
    for i=1:length(spotData);
        try
        spotMatrices{n} = [spotData{i}{n} ,spotMatrices{n}];
        spotDaxMatrices{n} = [spotDaxData{i}{n} ,spotDaxMatrices{n}];
        catch
        end
    end
end

%
if length(minPhotsPerStain) == 1;
minPhotsPerStain = repmat(minPhotsPerStain,1,numHybes);
end

correctLocsPerHybe = NaN*zeros(numHybes,numRealGenes); 
allMissedBits = NaN*zeros(numHybes,numRealGenes); 
allGainedBits = NaN*zeros(numHybes,numRealGenes); 
expectedBits = zeros(numHybes,1); 
cntsInHybe = zeros(numHybes,1);
detectRate  = NaN*zeros(numHybes,numRealGenes); 
detectRateN  = NaN*zeros(numHybes,numRealGenes); 
totSpots = 0;
totLocs = cell(numRealGenes,1); 
i = 0;
for n=realGenes % i = 3
    i=i+1;
    totLocs{i} = libCodes(n,:)*spotMatrices{n};
    numSpots = size(spotMatrices{n},2);
    correctLocs = repmat(libCodes(n,:)',1,numSpots) .* spotMatrices{n};
    locsTrueBits = correctLocs(libCodes(n,:),:);

    expectedBits = expectedBits + libCodes(n,:)'*numSpots; 
    
   % Find missed bits from corrected words
    missedBits = sum(correctLocs < repmat(minPhotsPerStain',1,numSpots),2);
    missedBits( ~libCodes(n,:),:) = 0;
    allMissedBits(:,i) = missedBits; 
    
    % find gained bits from corrected words
    gainedBits = repmat(~libCodes(n,:)',1,numSpots) .* spotMatrices{n};
    allGainedBits(:,i) = sum(gainedBits >= repmat(minPhotsPerStain',1,numSpots),2);
    correctLocsPerHybe(:,i) = mean(correctLocs, 2);
    
    onBits = find(libCodes(i,:));
    for j=onBits % i = 3; j = 1
        detectRate(j,i) = sum(spotMatrices{i}(j,:)>= minPhotsPerStain(j))/length(spotMatrices{i}(j,:));
        detectRateN(j,i) = numSpots*detectRate(j,i); 
        cntsInHybe(j) = cntsInHybe(j) + numSpots;
    end
    totSpots = totSpots + numSpots;
end
missedRate = nansum(allMissedBits,2)./expectedBits;
gainedRate = nansum(allGainedBits,2)./expectedBits;
detectionRate = nansum(detectRateN,2)./cntsInHybe;

% Compute Average and Min localizations among called spots
 max(correctLocsPerHybe,[],2)
 correctLocsPerHybe(correctLocsPerHybe==0)=NaN;
 correctAndCalled = correctLocsPerHybe;
 correctAndCalled(correctAndCalled <=  repmat(minPhotsPerStain',1,numRealGenes)  ) = NaN;
 minCorrectAndCalled = min(correctAndCalled,[],2);
 

% Same analysis on the Negative Controls 
errorLocsPerHybe = NaN*zeros(numHybes,numControls);
minErrorLocsPerHybe = NaN*zeros(numHybes,numControls);
errorLocs = cell(numControls,1); 
allErrorLocs = [];
k = 0; 
for i=cntrlGenes % i = 3
    k = k+1;
    errorLocs{k} = libCodes(i,:)*spotMatrices{i};
    numSpots = size(spotMatrices{i},2);
    errLocs = repmat(libCodes(i,:)',1,numSpots) .* spotMatrices{i};
    locsTrueBits = errLocs(libCodes(i,:),:);
    errorLocsPerHybe(:,k) = mean(errLocs, 2);
    allErrorLocs = cat(2,allErrorLocs,errLocs);
end
max(errorLocsPerHybe,[],2);
errorLocsPerHybe(errorLocsPerHybe==0)=NaN;
min(errorLocsPerHybe,[],2);
nanmean(errorLocsPerHybe,2);

 errorLocsPerHybe(errorLocsPerHybe==0)=NaN;
 wrongAndCalled = errorLocsPerHybe;
 wrongAndCalled(wrongAndCalled <  repmat(minPhotsPerStain',1,numControls)  ) = NaN;
 minWrongAndCalled = min(wrongAndCalled,[],2);

 
 % Plotting
 
 % Plot missed and gained bits;
 gainLossBit = figure(2); clf; 
subplot(2,1,1); bar(sum(allMissedBits,2)); title('missing bit');
subplot(2,1,2); bar(sum(allGainedBits,2)); title('gained bit');
PresentationPlot();
export_fig(gainLossBit,[savePath,'gainLossBit.png']); 

 bitSummary = figure(2); clf; 
subplot(4,1,1); bar(hybeBrightness); title('median spot brightness');
 subplot(4,1,2); bar(sum(allMissedBits,2)); title('missing bit');
subplot(4,1,3); bar(sum(allGainedBits,2)); title('gained bit');
subplot(4,1,4); bar(expectedBits); title('expected counts per bit');
PresentationPlot();
set(gcf, 'Units','Inches','Position',[0,0,5,10],'PaperUnits', 'Inches', 'PaperSize', [5 10]);
export_fig(bitSummary,[savePath,'bitSummary.png']); 

 % Plot missed and gained bit rate;
 gainLossBitRate = figure(2); clf; 
subplot(2,1,1); bar(missedRate); title('missing bit rate');
subplot(2,1,2); bar(gainedRate); title('gained bit rate');
PresentationPlot();
export_fig(gainLossBitRate,[savePath,'gainLossBitRate.png']); 

figDetectionSummary = figure(6); clf; 
subplot(2,1,1); bar( nanmean(correctLocsPerHybe,2) );
aveLocs = nanmean( nanmean(correctLocsPerHybe,2) );
title(['brightness per hybe in mRNA clusters. Ave=',num2str(aveLocs,3)]); 
xlabel('hybe number'); ylabel('brightness');

subplot(2,1,2); bar(detectionRate);
xlabel('hybe number');
ylabel('detection efficiency'); 
title('detection efficiency');
PresentationPlot();
export_fig(figDetectionSummary,[savePath,'figDetectionSummary.png']);

rnaCntsFig = figure(7); clf; set(gcf,'color','w'); 
imagesc(log10(rnaCnts)); cbar = colorbar;  colormap(hot(256)); 
set(gca,'YTick',1:numGenes,'YTickLabel',libGenes);
xlabel('cell','FontSize',14); ylabel('gene','FontSize',14); 
ylabel(cbar,'log(10) gene counts','FontSize',14);
title('Error Corrected Counts','FontSize',15);
set(gcf, 'Units','Inches','Position',[0,0,10,10],'PaperUnits', 'Inches', 'PaperSize', [10 10]);
export_fig(rnaCntsFig,[savePath,'rnaCntsFig.png']);

perCntsFig = figure(8); clf; set(gcf,'color','w'); 
imagesc(log10(perCnts)); cbar = colorbar;  colormap(hot(256)); 
set(gca,'YTick',1:numGenes,'YTickLabel',libGenes);
xlabel('cell','FontSize',14); ylabel('gene','FontSize',14); 
ylabel(cbar,'log(10) gene counts','FontSize',14);
title('Perfect Match Counts','FontSize',15);
set(gcf, 'Units','Inches','Position',[0,0,10,10],'PaperUnits', 'Inches', 'PaperSize', [10 10]);
export_fig(perCntsFig,[savePath,'perCntsFig.png']);

% save([savePath,'brightnessData.mat'],'spotMatrices','libGenes');

%% Spot Matrices analysis


close all;
PlotSpotMatrices(savePath,spotMatrices,libGenes,'flag','dao');
PlotSpotMatrices(savePath,spotDaxMatrices,libGenes,'flag','dax');

