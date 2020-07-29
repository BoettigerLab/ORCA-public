function [wordMap,sepMat] = MatchPixelsToWords(alignedIm,libCodes,c,varargin)
% Child function of Decode Pixels (step 2)

global figureSavePath;

troubleshoot = false;

defaults = cell(0,3);
defaults(end+1,:) = {'minDistFromZeros', 'nonnegative', .001}; % .004 Distance from all zero bit (speeds things up to make this cut before brightness cut0.
defaults(end+1,:) = {'minSeparationFromNext', 'nonnegative', .0001}; % Distance from next nearest bit
defaults(end+1,:) = {'quantileBlankBitRatioCut', 'fraction', .9}; % 
defaults(end+1,:) = {'quantileGeneBrightnessCut', 'fractions', .3}; % 
defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; % 
defaults(end+1,:) = {'moreVerbose', 'boolean', true}; % 
defaults(end+1,:) = {'overwrite', 'boolean', true}; % 
defaults(end+1,:) = {'savePath', 'string', figureSavePath}; % 
defaults(end+1,:) = {'idxBlanks', 'integer', []}; % 
defaults(end+1,:) = {'isGoodGene', 'integer', []}; % 
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% parameters = ParseVariableArguments([], defaults);

% required for blank cuts;
idxBlanks = parameters.idxBlanks;
isGoodGene = parameters.isGoodGene;


%% Main Function
 colordef white;
 testIm = double(alignedIm);
 if troubleshoot
    upsample = 3;
    % boxY = round((801:845)/upsample); 
    % boxX = round((1116:1160)/upsample);
    boxY = 261:290;     boxX = 301:330;
    testIm = double(alignedIm(boxY,boxX,:));
 end

testIm = testIm - min(testIm(:));
testIm = testIm./max(testIm(:));
[yDim,xDim,numBits] = size(testIm);
[numGenes,~] = size(libCodes);
if troubleshoot 
    figure(4); clf; Ncolor(testIm,hsv(numBits));
    colordef black;
    figure(5); clf; TileImage(uint16(2^16*testIm),'gain',2);
    [h,w] = size(testIm);
    set(gca,'XTick',1:5:4*w,'YTick',1:5:4*h);
    grid on; 
end

imWords = reshape(testIm,xDim*yDim,numBits); % dimenstions checked :)
% figure(6); clf; imagesc(imWords); colormap(jet(256)); caxis([0,1]); colorbar;
[dists,idx] = pdist2([double(libCodes)],imWords,'euclidean','Smallest',2);
distsJunk = pdist2([zeros(1,numBits); ones(1,numBits)],imWords,'euclidean');

idx(1,distsJunk(1,:)< parameters.minDistFromZeros) = numGenes+1;

sep = (dists(2,:)./dists(1,:)) -1; % separation from nearest alternative match
sepMat = reshape(sep,yDim,xDim);
% figure(6); clf; imagesc(sepMat); colormap hot; colorbar; set(gcf,'color','w'); 
reject = sep < parameters.minSeparationFromNext;
idx(1,reject) = numGenes+1;


bitRatio = cell(numGenes,1);
bitBrightness = cell(numGenes,1);
bitIDs = cell(numGenes,1);
for i=1:numGenes
    geneWords = imWords(idx(1,:)==i,:);
     bitRatio{i} = nanmedian(geneWords(:,libCodes(i,:)),2)./(nanmedian(geneWords(:,~libCodes(i,:)),2)+.0001); % avoiding inf errors
     bitBrightness{i} = nanmedian(geneWords(:,libCodes(i,:)),2);
     bitIDs{i} = find(idx(1,:)==i)';
end

blankBitRatio = cat(1,bitRatio{idxBlanks});
geneBitRatio = cat(1,bitRatio{isGoodGene});
blankBrightness = cat(1,bitBrightness{idxBlanks});
geneBrightness = cat(1,bitBrightness{isGoodGene});

minBitRatio = quantile(blankBitRatio,parameters.quantileBlankBitRatioCut);
minBrightness = quantile(geneBrightness,parameters.quantileGeneBrightnessCut);
rejectPixels = false(xDim*yDim,1); 
for i=1:numGenes
    idxReject = bitRatio{i} < minBitRatio;
    pixelReject = bitIDs{i}(idxReject);
    rejectPixels(pixelReject) = true;
    
    idxReject = bitBrightness{i} < minBrightness;
    pixelReject = bitIDs{i}(idxReject);
    rejectPixels(pixelReject) = true;
end
% sum(rejectPixels)/length(rejectPixels);
geneIdx = idx(1,:); 
geneIdx(rejectPixels) = numGenes + 1; 

wordMap = reshape(geneIdx,yDim,xDim); 



% just plotting
if parameters.saveFigures || strcmp(parameters.figVis,'on');
    % posGenes = cat(1,libGenes,'bkd');  % just for troubleshooting, requires libGenes 
    bitRatioImageBlanks = zeros(yDim,xDim);
    bitRatioImageBlanks( cat(1,bitIDs{idxBlanks}) ) = blankBitRatio;
    bitRatioImageGenes = zeros(yDim,xDim);
    bitRatioImageGenes( cat(1,bitIDs{isGoodGene}) ) = geneBitRatio;
    bitBrigtnessImageBlanks = zeros(yDim,xDim);
    bitBrigtnessImageBlanks( cat(1,bitIDs{idxBlanks}) ) = blankBrightness;
    bitBrigtnessImageGenes = zeros(yDim,xDim);
    bitBrigtnessImageGenes( cat(1,bitIDs{isGoodGene}) ) = geneBrightness;
    maxDisplayRatio = max(cat(1,bitRatio{:}))/4;
    maxDisplayBrightness = max(cat(1,bitBrightness{:}))/4;
    imageBlankComp = figure('Name','imageBlankComp','visible',parameters.figVis); clf; 
    subplot(2,2,1);  imagesc(bitRatioImageGenes); colorbar; caxis([0,maxDisplayRatio]); title('bit ratio genes'); % title(['bitRatio ',libGenes{i}]);
    subplot(2,2,2);  imagesc(bitRatioImageBlanks); colorbar; caxis([0,maxDisplayRatio]); title('bit ratio blanks'); %  title(['bitRatio ',libGenes{i}]);
    subplot(2,2,3);  imagesc(bitBrigtnessImageGenes); colorbar; caxis([0,maxDisplayBrightness]); title('brightness genes'); % title(['bitBrightness ',libGenes{i}]);
    subplot(2,2,4);  imagesc(bitBrigtnessImageBlanks); colorbar; caxis([0,maxDisplayBrightness]); title('brightness blanks'); % title(['bitBrightness ',libGenes{i}]);
    colormap(hot(256)); 
    set(gcf,'Units','Inches','Position',[.5 .5 10 10]);
    SaveFigure(imageBlankComp,'name',['imageBlankComp_',num2str(c,'%02d')],...
        'formats',{'png'},'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
        'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);

    figBlankComp = figure('Name','figBlankComp','visible',parameters.figVis); clf; 
    subplot(2,2,1); hist(blankBitRatio,linspace(0,maxDisplayRatio,100)); xlim([0,maxDisplayRatio]); 
    title(['bit ratio, blanks.  95%=',num2str(quantile(blankBitRatio,.95),3) ]);
    subplot(2,2,3); hist(geneBitRatio,linspace(0,maxDisplayRatio,100)); xlim([0,maxDisplayRatio]);
    title(['bit ratio, good genes. 50%=',num2str(quantile(geneBitRatio,.50),3) ]);
    subplot(2,2,2); hist(blankBrightness,linspace(0,maxDisplayBrightness,100)); xlim([0,maxDisplayBrightness]); 
    title(['bit brightenss, blanks.  95%=',num2str(quantile(blankBrightness,.95),3) ]);
    subplot(2,2,4); hist(geneBrightness,linspace(0,maxDisplayBrightness,100)); xlim([0,maxDisplayBrightness]); 
    title(['bit brightness, good genes. 50%=',num2str(quantile(geneBrightness,.50),3) ]);
    set(gcf,'Units','Inches','Position',[5,3,9,7]);
    SaveFigure(figBlankComp,'name',['figBlankComp_',num2str(c,'%02d')],...
        'formats',{'png'},'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
        'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);

    % just plotting
    wordMapFig = figure('Name','wordMapFig','visible',parameters.figVis); clf; 
    imagesc(wordMap);
    [X,Y] = meshgrid(1:xDim,1:yDim);    
    colormap([jet(256);zeros(1,3)]); hold on; colorbar;
    set(gcf,'Units','Inches','Position',[.5 .5 10 10]);
    if length(X(:)) < 1E3 && troubleshoot;  
        pixNames = posGenes(geneIdx);  text(X(:),Y(:),pixNames,'FontSize',5,'color','w');  
    end%  
    SaveFigure(wordMapFig,'name',['wordMapFig_',num2str(c,'%02d')],...
        'formats',{'png','fig'},'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
        'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);
end