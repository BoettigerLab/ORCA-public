function DecodePixels(dataToAnalyze,varargin)

% globals
global figureSavePath;

% parameters
defaults = cell(0,3);
defaults(end+1,:) = {'minDistFromZeros', 'nonnegative', .001}; % .004 Distance from all zero bit (speeds things up to make this cut before brightness cut0.
defaults(end+1,:) = {'minSeparationFromNext', 'nonnegative', .0001}; % Distance from next nearest bit
defaults(end+1,:) = {'quantileBlankBitRatioCut', 'fraction', .9}; % 
defaults(end+1,:) = {'quantileGeneBrightnessCut', 'fractions', .3}; % 
defaults(end+1,:) = {'borderSize', 'nonnegative', 45}; %  % down scale this if we drop upsampling.
defaults(end+1,:) = {'minPixels', 'nonnegative', 8}; %  % down scale this if we drop upsampling.
defaults(end+1,:) = {'calcSpotMat', 'boolean', false}; %  
defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; % 
defaults(end+1,:) = {'moreVerbose', 'boolean', true}; % 
defaults(end+1,:) = {'overwrite', 'boolean', true}; % 
defaults(end+1,:) = {'resumeRun', 'boolean', false}; % 
defaults(end+1,:) = {'startCell', 'integer', 1}; % 
defaults(end+1,:) = {'overwriteCellData', 'boolean', false}; % 
defaults(end+1,:) = {'maxCells', 'nonnegative', inf}; % 
defaults(end+1,:) = {'useParPool', 'boolean', false}; % 
defaults(end+1,:) = {'idxBits','array',[]};
defaults(end+1,:) = {'bitNames','array',{}};
defaults(end+1,:) = {'FPKMData','array',[]};
defaults(end+1,:) = {'finalDaxName','string','*Balanced*'};

parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main function

% 
if parameters.useParPool
    currPool = gcp;
    if ~isempty(currPool);  
        delete(currPool);  
    end
    parpool(20);
end

% maintaining backwards compatability with previous L7/L8 pipeline.  This can get cleaned up at some future timepoint.  
idxBits = parameters.idxBits;
bitNames = parameters.bitNames;
if isempty(idxBits)  % just for old L7 L8 pipepline
    idxBits = 1:dataToAnalyze.numHybs;
    if strcmp(dataToAnalyze.hybOrder,'rev');
        idxBits = fliplr(idxBits);
    end
    bitNames = cellstr(num2str(idxBits')); % will be used later to record bit and color channel 
end

if isempty(dataToAnalyze.FPKMData)
   dataToAnalyze.FPKMData = LoadByteStream(dataToAnalyze.FPKMDataPath);
end

%  Load codebook & compute expected counts from FPKM data
codebook =  dataToAnalyze.codebookPath; % [MERFISHdrive,'Libraries\Library11\V1\']; %'\\cajal\TSTORMdata\Library9\Lib9_140912\'; % '\\cajal\TSTORMdata\Library8\Lib8_140831\';
[libGenes,libExpect,libCodes,idxBlanks,isGene] = GetLibExpect(codebook,dataToAnalyze.FPKMData); % returns in descending FPKM order 
isGoodGene = 1:50; % use top 50-most expressed genes as valid count benchmark  
if length(idxBlanks) > 20
    idxBlanks = idxBlanks(1:20);
end
numGenes = length(libGenes);

% L7 data gets spatial image plots
exportSpatialImage = dataToAnalyze.libraryNum == 7;

%% Main function

newDataPath =dataToAnalyze.processedDataFolder;
temp = dir([newDataPath,parameters.finalDaxName,'.dax']);
daxNormName = {temp.name}';
numCells = min(length(daxNormName),parameters.maxCells);

startCell = parameters.startCell;
if parameters.resumeRun
    try 
        cellData = LoadByteStream([figureSavePath,'cellData.matb']);
        startCell = max([cellData.cellNum])+1;
        if startCell < numCells && startCell > 1
            disp(['Resuming analysis from cell ',num2str(startCell)]);
        end
    catch
    end
end

% delete(gcp);
% parpool(16);

if startCell < numCells    
    cellData(numCells).processedDataPath = newDataPath;  % structure to store cell data; 
    verbose = parameters.verbose; 
    geneDatas = cell(numCells,1);
    sepMats = cell(numCells,1);
    
    %% Loop over cells
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %% parfor loop
    if parameters.useParPool
        savePath = figureSavePath;
        parfor c = startCell:numCells;  % parfor 
            if  verbose;
                disp(['Analyzing Cell ',num2str(c)]);
            end

            % Load data
            cellName = [newDataPath,daxNormName{c}];
            alignedIm = LoadAndSortCellImage(cellName,idxBits,c,'parameters',parameters,'savePath',savePath);

            % Decode pixels
            [wordMap,sepMat] = MatchPixelsToWords(alignedIm,libCodes,c,'parameters',parameters,'savePath',savePath,...
                'isGoodGene',isGoodGene,'idxBlanks',idxBlanks);  % required if using comparison to blanks

            % Cluster pixels
            geneData = ClusterPixels(wordMap,alignedIm,libGenes,c,'parameters',parameters,'savePath',savePath,...
                'exportSpatialImage',exportSpatialImage) ;
            geneDatas{c} = geneData;
            sepMats{c} = sepMat;
        end % end loop over cells
        
        % sort things back into structure arrays
        for c=startCell:numCells;
            cellData(c).daxName = daxNormName{c};
            cellData(c).processedDataPath = newDataPath;  %
            cellData(c).cellNum = c;
            cellData(c).geneData =geneDatas{c};
            cellData(c).sepMat =sepMats{c};
        end
        disp('finished loop over cells'); 
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    else
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %% regular for  loop
        for c = startCell:numCells;  % parfor 
            if  verbose;
                disp(['Analyzing Cell ',num2str(c)]);
            end

            % Load data
            cellName = [newDataPath,daxNormName{c}];
            alignedIm = LoadAndSortCellImage(cellName,idxBits,c,'parameters',parameters);

            % Decode pixels
            [wordMap,sepMat] = MatchPixelsToWords(alignedIm,libCodes,c,'parameters',parameters,...
                'isGoodGene',isGoodGene,'idxBlanks',idxBlanks);  % required if using comparison to blanks

            % Cluster pixels
            geneData = ClusterPixels(wordMap,alignedIm,libGenes,c,'parameters',parameters,...
                'exportSpatialImage',exportSpatialImage) ;
            geneDatas{c} = geneData;
            sepMats{c} = sepMat;
        end % end loop over cells
        
        % sort things back into structure arrays
        for c=startCell:numCells;
            cellData(c).daxName = daxNormName{c};
            cellData(c).processedDataPath = newDataPath;  %
            cellData(c).cellNum = c;
            cellData(c).geneData =geneDatas{c};
            cellData(c).sepMat =sepMats{c};
        end
        disp('finished loop over cells'); 
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end


    %% combine all cell data for report
    allCellPerGeneArea = zeros(numGenes,numCells);
    allCellPerGeneCnt = zeros(numGenes,numCells);
    try
        for c=1:numCells
            allCellPerGeneArea(:,c) = cellfun(@sum,{cellData(c).geneData.area})  ;
            allCellPerGeneCnt(:,c) = cellfun(@length,{cellData(c).geneData.area})  ;
        end
        catch
    end
    
    saveName = [figureSavePath,'cellData.matb'];
    SaveAsByteStream(saveName,cellData,'overwrite',parameters.overwriteCellData);
    
     expPerCell = figure('Name','expPerCell','visible',parameters.figVis); clf; 
     set(gcf,'Units','Inches','Position',[.5 .5 3 12]);
     totCnts = sum(allCellPerGeneCnt,2);
     totAreas = sum(allCellPerGeneArea,2);
     subplot(1,2,1); barh(totCnts); set(gca,'YTick',1:numGenes,'YTickLabel',libGenes,'FontSize',5); title('spot counts','FontSize',12);
     subplot(1,2,2); barh(totAreas); set(gca,'YTick',1:numGenes,'YTickLabel',libGenes,'FontSize',5); title('total area','FontSize',12);
    SaveFigure(expPerCell,'name',['countsAllCells'],'formats',{'png','fig'},...
        'overwrite',parameters.overwrite,'saveData',1,...
            'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose);

    corrAreaFig = figure(50); clf; PlotCorr(libExpect',totAreas,libGenes); xlabel('FPKM'); ylabel('area');
    corrCntFig = figure(51); clf; PlotCorr(libExpect',totCnts,libGenes); xlabel('FPKM'); ylabel('count');
    SaveFigure(corrAreaFig,'name',['corrAreaFig'],'formats',{'png','eps','fig'},...
        'overwrite',parameters.overwrite,'saveData',1,...
            'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose);
    SaveFigure(corrCntFig,'name',['corrCntFig'],'formats',{'png','eps','fig'},...
        'overwrite',parameters.overwrite,'saveData',1,...
            'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose);

    if parameters.verbose
        meanPlusStdBlank = mean(totCnts(idxBlanks)) + 0*std(totCnts(idxBlanks));
        numGenesLessThanMeanPlusStd = sum(totCnts(isGene) < meanPlusStdBlank);
        disp(['# genes less than mean +/- std blanks (',num2str(meanPlusStdBlank),') = ',num2str(numGenesLessThanMeanPlusStd)]);

        maxBlank = max(totCnts(idxBlanks));
        numGenesLessThanMaxBlank = sum(totCnts(isGene) < maxBlank);
        disp(['# genes less than max blank (',num2str(maxBlank),') = ',num2str(numGenesLessThanMaxBlank)]);
    end


else
     disp('Data already analyzed! Skipping DecodePixels...');
end
