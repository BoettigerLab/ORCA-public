function [movieNW,infoNW,parameters] = CorrectDrift(movieN,infoPerCell,movieFid,newDataPath,varargin)
% Correct drift and chormatic aberration
% Uses either correlation based correction or fiducial mapping
% 

% parameters
defaults = cell(0,3);
defaults(end+1,:) = {'upsampleWarp', 'positive', 3}; % for correlation based drift alignment 
defaults(end+1,:) = {'minBeadContrast', 'fraction', .5}; % for correlation based drift alignment 
defaults(end+1,:) = {'numChns', 'positive', 1}; % for parsing multiple frames
defaults(end+1,:) = {'chromaticWarpFile', 'string', ''}; % chromewarps.mat file containing data from fidicial warp  
defaults(end+1,:) = {'chnNames', 'cell', {}}; % color channel names, match names in chromewarps.mat for the fidicial warp  s
defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
defaults(end+1,:) = {'showImages', 'boolean', false}; % 
defaults(end+1,:) = {'overwriteDax', 'boolean', false}; % 
defaults(end+1,:) = {'corrAlignImages', 'boolean', true}; % 
defaults(end+1,:) = {'fastAlign', 'boolean', true}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; 
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function


[nRows,nCols,numCells,numBits] = size(movieN);
% use just the center of the field of view for the correlation alignment
if parameters.fastAlign
    cropRows = round(.25*nRows):round(.75*nRows);
    cropCols = round(.25*nCols):round(.75*nCols);
else
    cropRows = 1:nRows;
    cropCols = 1:nCols;
end
movieFidwarped = movieFid;
movieNW = movieN; % just wapred; 

daxNoNormName = cell(numCells,1);
fidName = cell(numCells,1);

infoNW = infoPerCell;
infoFid= infoPerCell;

for c=1:numCells 
    daxNoNormName{c} = ['ImagesCell_',sprintf('%06d',c),'.dax'];
    fidName{c} = ['ImageFiducialsCell_',sprintf('%06d',c),'.dax'];
    daxExists = ~isempty( dir([newDataPath,daxNoNormName{c}]) );
    fidExists = ~isempty( dir([newDataPath,fidName{c}]) );
    if parameters.verbose 
        disp(['warping and saving data for cell ',num2str(c)]);
    end
    testData = movieN(:,:,c,1);
    notEmpty = sum(testData(:)) > 0;
    
    if (~(daxExists && fidExists) || parameters.overwriteDax) && notEmpty    
        h1 = imadjust(movieFid(cropRows,cropCols,c,1),stretchlim(movieFid(cropRows,cropCols,c,1),[parameters.minBeadContrast,1])); % 
        hcurr = 0;
        for b=1:numBits
            if parameters.corrAlignImages
                h = 1+ floor((b-.1)/parameters.numChns); % assumes each channel of same hybe uses same fiducial for drift / stage alignment 
                if h ~= hcurr
                    hH = imadjust(movieFid(cropRows,cropCols,c,h),stretchlim(movieFid(cropRows,cropCols,c,h),[parameters.minBeadContrast,1]));
                    figAlign = figure('Name','figAlign','visible',parameters.figVis); clf; 
                    [xshift,yshift] = CorrAlign(h1,hH,'region',200,...
                        'showplot',parameters.saveFigures,'upsample',parameters.upsampleWarp);
                    set(figAlign,'Units','Inches','Position',[2,.1,18,5],'color','w');
                    title(['xshift=',num2str(xshift,3),' yshift=',num2str(yshift,3)]);
                    SaveFigure(figAlign,'name',['Hybe',sprintf('%02d',h),'_Cell',sprintf('%02d',c),'_Fid'],...
                        'formats',{'png'},'overwrite',true,'verbose',false,...
                        'closeFig',parameters.closeOnSave,'saveData',parameters.saveFigures);
                    hcurr = h; 
                end
                movieNW(:,:,c,b) = TranslateImage(movieN(:,:,c,b),xshift,yshift,'upsample',parameters.upsampleWarp); 
                movieFidwarped(:,:,c,h) = TranslateImage(movieFid(:,:,c,h),xshift,yshift,'upsample',parameters.upsampleWarp);
            elseif parameters.fiducialAlign
                
                   Warp2BestPair(fedPos{1},fedPos{i}, ...
                    'parameters', parameters, ...
                    'showPlots', false); 
            else
                movieNW(:,:,c,b) = movieN(:,:,c,b); 
            end
            
            if ~isempty(parameters.chromaticWarpFile)
                 movieNW(:,:,c,b) = WarpImage(convBk,parameters.chnNames{b},...
                                         parameters.chromaticWarpFile,...
                                         'verbose',parameters.verbose);
            end

            if b==1  % update info file for this cell
                infoN = infoPerCell(c);
                infoN.localPath = newDataPath; 
                infoN.number_of_frames = numBits;
            end

            if parameters.showImages
                figure(1); clf; 
                subplot(1,2,1); imagesc(movieNW(:,:,c,b)); colorbar; colormap(gray(256));
                subplot(1,2,2); imagesc(temp); colorbar; colormap(gray(256));
                set(gcf,'Units','Inches','Position',[6,5.5,15,5],'color','w');
            end   
        end
        
        % Write dax files
        infoNW(c) = infoN; % info file for raw data file
        infoNW(c).localName = regexprep(daxNoNormName{c},'.dax','.inf'); 
        daxOutNoNorm = squeeze(movieNW(:,:,c,:)); 
        WriteDAXFiles(daxOutNoNorm,infoNW(c),'verbose',false); 
        
        infoFid(c) = infoN; % info file for fiducial dax
        infoFid(c).number_of_frames = 1;
        infoFid(c).localName = regexprep(fidName{c},'.dax','.inf'); 
        fiducialDax = squeeze(movieFidwarped(:,:,c,:));
        WriteDAXFiles(fiducialDax,infoFid(c),'verbose',false);              
    else
        if parameters.verbose
            disp(['Skipping warp for cell ',num2str(c),' Using existing file...']);
        end
    end
end

