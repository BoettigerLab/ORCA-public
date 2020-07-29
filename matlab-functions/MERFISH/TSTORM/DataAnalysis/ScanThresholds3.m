function [bestThetas,initRho,finalRho] = ScanThresholds3(cellProperties,thresholds,initThresh,newPath,varargin)
% [bestThetas,initRho,finalRho] = ScanThresholds2(cellProperties,thresholds,initThresh,newPath,varargin)
% 

%--------------------------------------------------------------------------
% Optional Input Parameters
%--------------------------------------------------------------------------

% general parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'scanThresholds','boolean',true};
defaults(end+1,:) = {'showFigs','string','on'};

% Parameters for spotfitting
defaults(end+1,:) = {'maxDtoCentroid', 'positive', 1}; 
defaults(end+1,:) = {'minPhotsPerStain', 'positive',0};
defaults(end+1,:) = {'bitOrder','array',[]};

% Parameters for Extract Code And FPKM
defaults(end+1,:) = {'FPKMData', 'struct', []}; % Delimiters for bin files
defaults(end+1,:) = {'codebookPath', 'string', ''}; 
defaults(end+1,:) = {'numHybs', 'positive',[]};

% parameters for estimating hybe accuracy
defaults(end+1,:) = {'method', {'minErrorRatePerSpot','minErrorRatePerGene'},'minErrorRatePerGene'};
defaults(end+1,:) = {'oneWeight','nonnegative',0};

parameters = ParseVariableArguments(varargin, defaults, mfilename);
if isempty(parameters.FPKMData)
    error('parameters.FPKMData is required');
end
if isempty(parameters.codebookPath)
    error('parameters.codebookPath is required');
end
if isempty(parameters.numHybs)
    error('parameters.numHybs is required');
end

%% Threshold scanning!

[libGenes,libCodes,libExpect] = ExtractCodeAndFPKM(parameters.FPKMData,parameters.codebookPath,parameters.numHybs);


if parameters.bitOrder(1) < parameters.bitOrder(end)
   reverseBits = false;
else
   reverseBits = true;
end

oneWeight = parameters.oneWeight;
numIterations = 5; % max numIterations
numCells = length(cellProperties);
numThresholds = length(thresholds); 
[numGenes,numHybes] = size(libCodes);
numThreshCombos = numHybes*numThresholds*numIterations;

% indexes for blanks and real genes
idxBlank = StringFind(libGenes,'blank');
idxNoTarget = StringFind(libGenes,'notarget');
idxBlank = [idxBlank; idxNoTarget];
idxGene = 1:numGenes;
idxGene(idxBlank) = [];

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
initRho = 0; 

if parameters.scanThresholds; 
    % initialize data storage arrays;
    rhoA = NaN(numHybes,numThresholds,numIterations); 
    matchRate = NaN(numHybes,numThresholds,numIterations); 
    
    perCnts = NaN(numGenes,numTestCells,numHybes,numThresholds,numIterations);
    recThetas = NaN(numHybes,numIterations,numHybes,numThresholds); 
    
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
                recThetas(:,i,h,j) = thetas'; 

                spotMatrix = cell(numTestCells,1); 
                k=0; 
                wordsDetected = cell(numTestCells,1);
                noSpots = zeros(1,numHybes); 
                 
                allCnts = NaN(numGenes,numTestCells); 
                for c=testCells
                    k = k+1; 
                    try
                        idx = sub2ind([numHybes,numThresholds], (1:numHybes)',thetas'); 
                        imLists = cellProperties(c).imLists(idx);
                        warps = cellProperties(c).warps; 
                        [xf,yf,flists] = FilterConvData(imLists,{warps.tform},'showPlots',false,'endframe',inf);

                        spotPositions = [cat(1,xf{:}),cat(1,yf{:})];
                        [wordsDetected{k},~,brightnessPerSpot] = FindWords(flists,spotPositions,...
                                                        'maxDtoCentroid',parameters.maxDtoCentroid,...
                                                        'minPhotsPerStain',parameters.minPhotsPerStain );   

                        if reverseBits
                            wordsDetected{k} = fliplr(wordsDetected{k}); 
                            brightnessPerSpot = fliplr(brightnessPerSpot);
                        end
                            [perCnts(:,k,h,j,i),allCnts(:,k),spotMatrix{k}] =...
                                        DecodeConvImage(wordsDetected{k},libCodes,...
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
                
                % 
                totBlanks = sum(sum(perCnts(idxBlank,:,h,j,i)));
                totGenes =  sum(sum(perCnts(idxGene,:,h,j,i)));
                sumGenes = sort(nansum(perCnts(idxGene,:,h,j,i),2));
                sumBlanks =  sum(perCnts(idxBlank,:,h,j,i),2);
                maxBlank =  max(sumBlanks);
                msBlank =  mean(sumBlanks) + std(sumBlanks);
                
                
                
                if parameters.verbose
                    disp([num2str( 100*t/numThreshCombos),'% complete...']);  
                    disp(['totBlanks/totGenes ',num2str(totBlanks),'/',num2str(totGenes)]);
                    disp(['genes > maxBlank: ', num2str( sum(sumGenes > maxBlank) )]);
                    disp(['genes > msBlank: ', num2str( sum(sumGenes > msBlank) )]);     
                    disp(['max match rate = ',num2str(maxMatch)]);
                    disp(['rho at max match rate = ',num2str(maxRho)]);
                end

                if isErrorCorrecting
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
                    
                    %% Num ones
                    allWordsDetected = cat(1,wordsDetected{:});
                    num1s = sum(sum(allWordsDetected,2) == 1);
                    oneScore = 1/(1+num1s);
                            
                    switch parameters.method
                        case 'minErrorRatePerGene'
                            [perHybeMissRate,perHybeGainRate] = GetBitFlipRates(libCodes,spotMat);
                             hybeAcc = 1 - mean(perHybeMissRate) - mean(perHybeGainRate);
                            matchRate(h,j,i) = hybeAcc + oneWeight*oneScore;
                        case 'minErrorRatePerSpot'
                            totSpotMat = cat(2,spotMat{:});
                            totSpotMat = logical(totSpotMat);
                            missedBitColumns = sum(totSpotMat,1) == 3;
                            gainedBitColumns = sum(totSpotMat,1) == 5;
                            % perfectBitColumns = sum(totSpotMat,1) == 4;
                            perHybeMissRate = sum(totSpotMat(:,missedBitColumns),2)/size(totSpotMat,2);
                            perHybeGainRate = sum(totSpotMat(:,gainedBitColumns),2)/size(totSpotMat,2);
                            hybeAcc =  1 - mean(perHybeMissRate) - mean(perHybeGainRate);
                            matchRate(h,j,i) = hybeAcc + oneWeight*oneScore;
                        case 'minVarErrorRate' % not functional
                            totSpotMat = cat(2,spotMat{:});
                            totSpotMat = logical(totSpotMat);
                            missedBitColumns = sum(totSpotMat,1) == 3;
                            gainedBitColumns = sum(totSpotMat,1) == 5;
                            % perfectBitColumns = sum(totSpotMat,1) == 4;
                            perHybeMissRate = sum(totSpotMat(:,missedBitColumns),2)/size(totSpotMat,2);
                            perHybeGainRate = sum(totSpotMat(:,gainedBitColumns),2)/size(totSpotMat,2);
                            std(perHybeMissRate) + std(perHybeGainRate);
                    end
                else % for non error correcting codes
                    matchRate(h,j,i) = totGenes/totBlanks;
                    
%                     allWordsDetected = cat(1,spotMatrix{:}); % overloaded spotMatrix with cell array of wordsDetected  
%                     num3s = sum(allWordsDetected,2) == 3;
%                     num4s = sum(allWordsDetected,2) == 4;
%                     num5s = sum(allWordsDetected,2) == 5;
%                     matchRate(h,j,i) = sum(num4s)/( sum(num3s)+sum(num4s)+sum(num5s) );
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
            try
                if strcmp(parameters.showFigs,'on');
                    figure(2); clf; subplot(2,1,1);
                    plot(rhoL,'r'); hold on; plot(hybeL,'b'); hold on;
                    legend('fpkm correlation','hybe accuracy','Location','SouthEast'); 
                    xlabel('iteration');  ylabel('FPKM rho (log_1_0)'); 
                    title(['max rho = ',num2str( max(rhoA(:)) )]);
                    subplot(2,1,2); imagesc(threshPerms); 
                    xlabel('iteration'); ylabel('threshold per hybe');
                    PresentationPlot(); pause(.1); 
                end
            
            catch er
                disp(er.getReport);
            end
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
else
    scanData = dir([newPath,'thresholdScanData*.mat']);
    [~,newest] = max( [scanData.datenum] );
    load([newPath,scanData(newest).name]);
end

try
rho_vs_matchCorr = figure('Name','rho_vs_matchCorr','visible',parameters.showFigs); clf; 
plot(rhoA(:),matchRate(:),'b.');
xlabel('FPKM correlation'); ylabel('hybe accuracy'); 
PresentationPlot; 
catch er
    warning(er.getReport);
end

try
iterationProgress = figure('Name','iterationProgress','visible',parameters.showFigs);
clf; subplot(2,1,1);
plot(rhoA(:),'r'); hold on; plot(matchRate(:),'b'); hold on;
legend('fpkm correlation','hybe accuracy','Location','SouthEast'); 
ylabel('FPKM rho (log_1_0)'); xlabel('iteration'); 
title(['max rho = ',num2str( max(rhoA(:)) )]);
subplot(2,1,2); imagesc(threshPerms); 
xlabel('iteration');
ylabel('threshold per hybe');
PresentationPlot();
catch er
    warning(er.getReport);
end
    
try
heatmapsRhoAndAccuracy = figure('Name','heatmapsRhoAndAccuracy','visible',parameters.showFigs); 
clf; subplot(2,1,1);
imagesc(reshape(rhoA,numHybes,numThresholds*numIterations));
ylabel('FPKM rho (log_1_0)'); xlabel('threshold combination'); 
title(['max rho = ',num2str( max(rhoA(:)) )]);
subplot(2,1,2);%  imagesc(threshPerms'); 
imagesc(reshape(matchRate,numHybes,numThresholds*numIterations));
xlabel('threshold combination');
ylabel('hybe');
PresentationPlot();
catch er
    warning(er.getReport)
end
    

% apply bestThetas to all cells
allCnts = NaN(numGenes,numCells);
for c=1:numCells
    idx = sub2ind([numHybes,numThresholds], (1:numHybes)',bestThetas'); 
    imLists = cellProperties(c).imLists(idx);
    warps = cellProperties(c).warps;
    [xf,yf,flists] = FilterConvData(imLists,{warps.tform},'showPlots',false,'endframe',inf);

    spotPositions = [cat(1,xf{:}),cat(1,yf{:})];
    [wordsDetected,~,brightnessPerSpot] = FindWords(flists,spotPositions,...
                                    'maxDtoCentroid',parameters.maxDtoCentroid,...
                                    'minPhotsPerStain',parameters.minPhotsPerStain );   

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
finalRho = corr(log10(x(non0 & nonNaN)),log10(y(non0 & nonNaN)));

corrFig = figure('Name','corrFig'); clf; PlotCorr(x,y,libGenes);

% Save DATA
% -------------------------------------------------------------
savename = IncrementSaveName([newPath,'thresholdScanData.mat']);
save(savename,'bestThetas','rhoA','matchRate','threshPerms','thresholds','testCells'); 

try
    SetFigureSavePath(newPath);
    overwrite = false;
    SaveFigure(corrFig, 'overwrite', overwrite, ...
                        'formats', parameters.figFormats);
    savename = IncrementSaveName([newPath,'rho_vs_matchCorr.png']);
    export_fig(rho_vs_matchCorr, savename);
    savename = IncrementSaveName([newPath,'rho_vs_matchCorr.fig']);
    saveas(rho_vs_matchCorr, savename);
    
    savename = IncrementSaveName([newPath,'iterationProgress.png']);
    export_fig(iterationProgress, savename);
    savename = IncrementSaveName([newPath,'iterationProgress.fig']);
    saveas(iterationProgress, savename);
    
    savename = IncrementSaveName([newPath,'heatmapsRhoAndAccuracy.png']);
    export_fig(heatmapsRhoAndAccuracy, savename);
    savename = IncrementSaveName([newPath,'heatmapsRhoAndAccuracy.fig']);
    saveas(heatmapsRhoAndAccuracy, savename);
    
    close corrFig rho_vs_matchCorr iterationProgress heatmapsRhoAndAccuracy;               
catch er
    warning(er.getReport);
    warning('error saving figures'); 
end

% export_fig(corrFig,[newPath,'corrFig.png']);
% export_fig(rho_vs_matchCorr,[newPath,'rho_vs_matchCorr.png']);
% export_fig(iterationProgress,[newPath,'iterationProgress.png']);
% export_fig(heatmapsRhoAndAccuracy,[newPath,'heatmapsRhoAndAccuracy.png']);
% saveas(corrFig,[newPath,'corrFig.fig']);
% saveas(rho_vs_matchCorr,[newPath,'rho_vs_matchCorr.fig']);
% saveas(iterationProgress,[newPath,'iterationProgress.fig']);
% saveas(heatmapsRhoAndAccuracy,[newPath,'heatmapsRhoAndAccuracy.fig']);
