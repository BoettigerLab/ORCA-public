function [bestThetas,initRho,finalRho] = ScanThresholds(cellProperties,warps,thresholds,codebook,probeData,reverseBits,initThresh,newPath)

%% Threshold scanning!
% Parameters
parameters.verbose = true; 

maxDtoCentroid = 1;
minPhotsPerStain = 0; 


[~,libExpect,~,libGenes,libCodes] = ExpectedLocsPerHybe2(codebook,probeData); 



numIterations = 5; % max numIterations
numCells = length(cellProperties);
numThresholds = length(thresholds); 
[numGenes,numHybes] = size(libCodes);
numThreshCombos = numHybes*numThresholds*numIterations;

% cells to train on
numTestCells = 15;% numCells;% % 5; 
testCells =  randi(numCells,1,numTestCells); % 1:numCells; % 

if numHybes == 14
    isErrorCorrecting = false;
elseif numHybes == 12 || numHybes == 16
    isErrorCorrecting = true;
end


% precompute nearest neighbor codes for speed
[~,~,~,~,~,allCorCodes] = DecodeConvImage(ones(numHybes,1),libCodes);

% initialize data storage arrays;
perCnts = NaN(numGenes,numTestCells); 
allCnts = NaN(numGenes,numTestCells); 
rhoA = NaN(numHybes,numThresholds,numIterations); 
matchRate = NaN(numHybes,numThresholds,numIterations); 

thetas = initThresh*ones(1,numHybes); 
maxRho = 0;
maxMatch = 0;
t = 0; 
rhoL = NaN(numThreshCombos,1); 
hybeL  = NaN(numThreshCombos,1); 
threshPerms = initThresh*ones(numHybes,numThreshCombos);
for i=1:numIterations;  
    for h=1:numHybes
        for j=1:numThresholds
            t=t+1;
            thetas(h) = j;
            threshPerms(:,t) = thetas';

            spotMatrix = cell(numTestCells,1); 
            k=0; 
            noSpots = zeros(1,numHybes); 
            for c=testCells
                k = k+1; 
                try
                    idx = sub2ind([numHybes,numThresholds], (1:numHybes)',thetas'); 
                    imLists = cellProperties(c).imLists(idx);
                    [xf,yf,flists] = FilterConvData(imLists,{warps{c}.tform},'showPlots',false,'endframe',inf);

                    spotPositions = [cat(1,xf{:}),cat(1,yf{:})];
                    [wordsDetected,~,brightnessPerSpot] = FindWords(flists,spotPositions,...
                                                    'maxDtoCentroid',maxDtoCentroid,...
                                                    'minPhotsPerStain',minPhotsPerStain );   

                    if reverseBits
                        wordsDetected = fliplr(wordsDetected); 
                        brightnessPerSpot = fliplr(brightnessPerSpot);
                    end
                        [perCnts(:,k),allCnts(:,k),spotMatrix{k}] =...
                                    DecodeConvImage(wordsDetected,libCodes,...
                                    'correctErrors',isErrorCorrecting,...
                                    'allCorCodes',allCorCodes,...
                                    'brightnessPerSpot',brightnessPerSpot);

                    % compute correlation to FPKM    
                    x = libExpect;  y =  allCnts(:,k);
                    non0 = x>0 & y>0;
                    nonNaN = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y) );
                    rho = corr(log10(x(non0 & nonNaN)),log10(y(non0 & nonNaN)));
                    disp(['theta ',num2str(t),' of ',num2str(numThreshCombos),' cell ',num2str(c),' of ',num2str(numTestCells)]);
                catch er
                    warning(er.getReport);
                    allCnts = zeros(numGenes,numTestCells);
                end
                    
                noSpots = noSpots + cellfun(@isempty, xf);
            end
            if parameters.verbose
                disp([num2str( 100*t/numThreshCombos),'% complete...']);   
                disp(['max match rate = ',num2str(maxMatch)]);
                disp(['rho at max match rate = ',num2str(maxRho)]);
            end
            
            if isErrorCorrecting
             % Combine spot matrices for each cell so we have one per gene
                spotMat = cell(numGenes,1);
                spotMatNorm = cell(numGenes,1);
                for n=1:numGenes
                    for a=1:numTestCells;
                        try
                        spotMat{n} = [spotMatrix{a}{n} ,spotMat{n}];
                        catch
                        end
                    end
                    numGeneSpots = size(spotMat{n},2);
                    spotMatNorm{n} = spotMat{n}./repmat(thresholds(thetas)',1,numGeneSpots);
                end   
                % figure(1); clf; imagesc(spotMat{1})
                [missedRate,gainedRate,~] = GetBitFlipRates(libCodes,spotMat);
                matchRate(h,j,i) = 1 - mean(missedRate) - mean(gainedRate);
            else
                
                
                
                matchRate(h,j,i) = 1 - mean(missedRate) - mean(gainedRate);
            end
            if sum(noSpots  > numTestCells/2) >0 % if more than half the cells have no spots in 1 or more hybes
                matchRate(h,j,i) = NaN;
            end
            
             % compute correlation to FPKM    
            x = libExpect;  y =  sum(allCnts(:,:),2);
            non0 = x>0 & y>0;
            nonNaN = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y) );
            rhoA(h,j,i) = corr(log10(x(non0 & nonNaN)),log10(y(non0 & nonNaN)));
            if t==1
                initRho = rhoA(h,j,i);
            end
            rhoL(t) = rhoA(h,j,i);
            hybeL(t) = matchRate(h,j,i); 
        end 
        
        [maxMatch,maxId] = max(matchRate(h,:,i));
        maxRho = rhoA(h,maxId,i);
        disp(['max match rate = ',num2str(maxMatch),' theta ID=',num2str(maxId)]);
        disp(['init Rho = ',num2str(initRho)]);
        disp(['rho at max match rate = ',num2str(maxRho)]);
     
        % update thresholds to max matchRate
        thetas(h) = maxId;
        
        % plot Progress
        figure(2); clf; subplot(2,1,1);
        plot(rhoL,'r'); hold on; plot(hybeL,'b'); hold on;
        legend('fpkm correlation','hybe accuracy','Location','SouthEast'); 
        xlabel('iteration');  ylabel('FPKM rho (log_1_0)'); 
        title(['max rho = ',num2str( max(rhoA(:)) )]);
        subplot(2,1,2); imagesc(threshPerms); 
        xlabel('iteration'); ylabel('threshold per hybe');
        PresentationPlot(); pause(.1); 
    end
    
    % determine if thresholds have converged or not
    disp('finished iteration'); 
    threshKept = NaN(numHybes,1);  
    if i>1
        for h=1:numHybes
            threshKept(h) = threshPerms(h,t) == threshPerms(h,(i-1)*numHybes*numThresholds);
        end
        if sum(threshKept) == numHybes
            disp(['thresholds converged at iteration ',num2str(i)])
            break
        end
    end
    
end

bestThetas = thetas; 


figure(2); clf; 
rho_vs_matchCorr = figure(1); clf; plot(rhoA(:),matchRate(:),'b.');
xlabel('FPKM correlation'); ylabel('hybe accuracy'); 
PresentationPlot; 


iterationProgress = figure(2); clf; subplot(2,1,1);
plot(rhoA(:),'r'); hold on; plot(matchRate(:),'b'); hold on;
legend('fpkm correlation','hybe accuracy','Location','SouthEast'); 
ylabel('FPKM rho (log_1_0)'); xlabel('iteration'); 
title(['max rho = ',num2str( max(rhoA(:)) )]);
subplot(2,1,2); imagesc(threshPerms); 
xlabel('iteration');
ylabel('threshold per hybe');
PresentationPlot();

heatmapsRhoAndAccuracy = figure(3); clf; subplot(2,1,1);
imagesc(reshape(rhoA,numHybes,numThresholds*numIterations));
ylabel('FPKM rho (log_1_0)'); xlabel('threshold combination'); 
title(['max rho = ',num2str( max(rhoA(:)) )]);
subplot(2,1,2);%  imagesc(threshPerms'); 
imagesc(reshape(matchRate,numHybes,numThresholds*numIterations));
xlabel('threshold combination');
ylabel('hybe');
PresentationPlot();

% apply bestThetas to all cells
allCnts = NaN(numGenes,numCells);
for c=1:numCells
    idx = sub2ind([numHybes,numThresholds], (1:numHybes)',bestThetas'); 
    imLists = cellProperties(c).imLists(idx);
    [xf,yf,flists] = FilterConvData(imLists,{warps{c}.tform},'showPlots',false,'endframe',inf);

    spotPositions = [cat(1,xf{:}),cat(1,yf{:})];
    [wordsDetected,~,brightnessPerSpot] = FindWords(flists,spotPositions,...
                                    'maxDtoCentroid',maxDtoCentroid,...
                                    'minPhotsPerStain',minPhotsPerStain );   

    if reverseBits
        wordsDetected = fliplr(wordsDetected); 
        brightnessPerSpot = fliplr(brightnessPerSpot);
    end
    [~,allCnts(:,c)] =...
                    DecodeConvImage(wordsDetected,libCodes,...
                    'correctErrors',isErrorCorrecting,...
                    'allCorCodes',allCorCodes,...
                    'brightnessPerSpot',brightnessPerSpot);

end

% compute correlation to FPKM    
x = libExpect;  y =  sum(allCnts,2);
non0 = x>0 & y>0;
nonNaN = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y) );
finalRho = corr(log10(x(non0 & nonNaN)),log10(y(non0 & nonNaN)))

corrFig = figure(10); clf; PlotCorr(x,y,libGenes);

% Save DATA
% -------------------------------------------------------------
save([newPath,'thresholdScanData.mat'],'bestThetas','rhoA','matchRate','threshPerms','thresholds'); 

export_fig(corrFig,[newPath,'corrFig.png']);
export_fig(rho_vs_matchCorr,[newPath,'rho_vs_matchCorr.png']);
export_fig(iterationProgress,[newPath,'iterationProgress.png']);
export_fig(heatmapsRhoAndAccuracy,[newPath,'heatmapsRhoAndAccuracy.png']);
saveas(corrFig,[newPath,'corrFig.fig']);
saveas(rho_vs_matchCorr,[newPath,'rho_vs_matchCorr.fig']);
saveas(iterationProgress,[newPath,'iterationProgress.fig']);
saveas(heatmapsRhoAndAccuracy,[newPath,'heatmapsRhoAndAccuracy.fig']);
