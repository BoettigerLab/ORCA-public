function [mosaic,pars,imData,imBk,stageXYc] = DaxToMosaic(daxFiles,varargin)
% Take a stack of filepaths to daxfiles, load the dax files and their info
% files, and arrange them in a mosaic tile.  


%% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'pix_to_mm','positive',6.55};  % 6.55 = scope 1. 6.45 = scope 2
defaults(end+1,:) = {'trimBorder','nonnegative',0};
defaults(end+1,:) = {'selectFrame','integer',2};
defaults(end+1,:) = {'projectAll','boolean',false};
defaults(end+1,:) = {'transpose','boolean',true}; % see the mosaic parameters in the Hal parameter file used to collect the movie 
defaults(end+1,:) = {'fliplr','boolean',true};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'flipud','boolean',false};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'troubleshoot','boolean',false};  
defaults(end+1,:) = {'rescale','positive',1}; %
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'showProgress','boolean',false};
defaults(end+1,:) = {'saveProjections','boolean',true};  % save the max projected images for quick loading 
defaults(end+1,:) = {'loadProjections','boolean',true};  % load previousl save the max projected images for quick loading 
% parameters to flatten background
defaults(end+1,:) = {'flatten','boolean',false};
defaults(end+1,:) = {'flatField','freeType',[]}; % pass array
defaults(end+1,:) = {'outOfFocus','integer',1};
defaults(end+1,:) = {'flattenBlur','positive',0};
% parameters for correlation based registration
defaults(end+1,:) = {'corrAlign','boolean',false};
defaults(end+1,:) = {'corrAlignHigh','boolean',.9999};
defaults(end+1,:) = {'corrAlignLow','boolean',.9};
defaults(end+1,:) = {'shiftsXY','freeType',[]};
defaults(end+1,:) = {'showplots','boolean',false};
defaults(end+1,:) = {'minOverlap','fraction',.05};
defaults(end+1,:) = {'maxShift','nonnegative',500};
defaults(end+1,:) = {'saveCorrFolder','string',''};
defaults(end+1,:) = {'minGrad', 'float', 6000}; % -inf  min gradient to be called a good alignment match 
% parameters for blendEdges
defaults(end+1,:) = {'blendEdges','boolean',false};
defaults(end+1,:) = {'blendRadius','integer',100};
defaults(end+1,:) = {'bufferEdge','integer',50};
% parameters useful for iterative mosaic overlays
defaults(end+1,:) = {'offset','float',0}; % this is useful to pass from previous mosaics to have multiple files in same register
defaults(end+1,:) = {'nRows','freeType',[]};
defaults(end+1,:) = {'nCols','freeType',[]};
defaults(end+1,:) = {'buffer','nonnegative',1000};
defaults(end+1,:) = {'sortIndex','freeType',[]};
% parameters for mosaic tile if data is already loaded
defaults(end+1,:) = {'imData','cell',{}}; % this doesn't work at all right now for reasons unclear 
defaults(end+1,:) = {'stageXY','array',[]}; % probably having to do with the shift here. 
defaults(end+1,:) = {'imBk','array',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% main function

if ~isempty(pars.saveCorrFolder)
    SetFigureSavePath(pars.saveCorrFolder);
    saveFigure = true;
else
    saveFigure = false;
end

imBk = [];

%% load data
if isempty(pars.imData)
    if pars.verbose
        disp('loading mosaic...');
    end

    % load data
    sc = pars.pix_to_mm; % 6.5 5
    numIms = length(daxFiles);
    imData = cell(numIms,1);
    imBkd = cell(numIms,1);
    stagePos = zeros(numIms,2);
    for f = 1:numIms
        if pars.projectAll
            [folder,filename] = fileparts(daxFiles{f});
            maxFile = [folder,filesep,'max_',filename,'.dax'];
            if exist(maxFile,'file')~=0 && pars.loadProjections  % use previously written
                [dax,infoFile] = ReadDax(maxFile,'verbose',false);
                if pars.flatten
                    imBkd{f} = dax;
                end
            else
                [dax,infoFile] = ReadDax(daxFiles{f},'verbose',false);
                 % select out of focus frames for flattening image background
                 if pars.flatten
                    imBkd{f} = median(dax(:,:,pars.outOfFocus),3);
                 end
                dax = max(dax,[],3); 
                if pars.saveProjections % save max projections
                    infoOut = infoFile;
                    infoOut.number_of_frames= 1;
                    infoOut.localName = ['max_',infoOut.localName];
                    WriteDAXFiles(dax,infoOut,'verbose',false);
                end
            end
            
        else
            [dax,infoFile] = ReadDax(daxFiles{f},'startFrame',pars.selectFrame,'endFrame',pars.selectFrame,'verbose',false);

            % select out of focus frames for flattening image background
            if pars.flatten
                temp = ReadDax(daxFiles{f},'verbose',false,'startFrame',min(pars.outOfFocus),'endFrame',max(pars.outOfFocus),'verbose',false);
                imBkd{f}  = median(temp,3); 
            end
        end
        imData{f} = dax;
        stagePos(f,:) = [infoFile.Stage_X,infoFile.Stage_Y]*sc;    
    end
    [h,w,~] = size(dax); % assumes all dax are the same size;
    pars.imHeight = h - 2*pars.trimBorder;
    pars.imWidth = w - 2*pars.trimBorder;

    if pars.flatten && isempty(pars.imBk)
        if isempty(pars.flatField)
            imBk = cat(3, imBkd{:}); 
            imBk = nanmedian( imBk,3);
        else 
            imBk = pars.flatField;
        end
        % figure(100); clf; imagesc(IncreaseContrast(imBk,'high',.9999));
    end

    % use optional offset
    if pars.offset == 0
        %  This is helpful to align multiple mosaics to the same stage coordinates
        pars.offset = - min(stagePos) + [h/2 w/2] + pars.buffer*[1 1];
    end
    
    % update stageXY
    stageXY = stagePos + pars.offset;
    
else % data is already loaded and was passed. 
    if pars.verbose
        disp('loading skipped. using passed data');
    end
    imData = pars.imData;
    stageXY = pars.stageXY;
    imBk = pars.imBk;
    numIms = length(imData); 
    [h,w] = size(imData{1});
end
%-------------------------------------------------------------------------%

%=========================================================================%
% array data in matrices
%=========================================================================%


if isempty(pars.nRows)
    %  This is helpful to align multiple mosaics to the same stage coordinates
    pars.nRows = ceil( max(stageXY(:,2)) + h/2+ pars.buffer*1 ); % 9/21
end
if isempty(pars.nCols)
    %  This is helpful to align multiple mosaics to the same stage coordinates
    pars.nCols = ceil( max(stageXY(:,1)) + w/2+ pars.buffer*1 ); % 9/21
end


% imData = imData(sortIndex); 
% shiftsXY = zeros(numIms,2);
% pars.sortIndex = sortIndex;

if isempty(pars.sortIndex)
    [sortIndex,~] = SortFOVsByOverlap([pars.nRows,pars.nCols], [pars.imWidth, pars.imHeight], stageXY);
    pars.sortIndex = sortIndex;
    stageXYc = stageXY(sortIndex,:);
    shiftsXY = zeros(numIms,2);
else
    sortIndex = pars.sortIndex;
    stageXYc = pars.stageXY(sortIndex,:);
    shiftsXY = pars.shiftsXY(sortIndex,:);
end

imData = imData(sortIndex); 



if ~pars.troubleshoot
    % mosaic = zeros(pars.nRows*pars.rescale,pars.nCols*pars.rescale,'uint16');
    mosaic  = nan(pars.nRows,pars.nCols,'double');
    mosaic2 = nan(pars.nRows,pars.nCols,'double');
    for f=1:numIms   
        dax = double(imData{f});        
        if pars.flatten
            dax = 2^8*double(dax)./double(imBk);
        end        
        if pars.transpose
            dax = dax';
        end
        if pars.fliplr
            dax = fliplr(dax);
        end
        if pars.flipud
            dax = flipud(dax);
        end
        if pars.trimBorder > 0
           dax(1:pars.trimBorder,:) = []; % nan
           dax(end-pars.trimBorder+1:end,:) = [];
           dax(:,1:pars.trimBorder) = [];
           dax(:,end-pars.trimBorder+1:end) = [];
        end
        rs = round(stageXYc(f,2) - h/2 +pars.trimBorder : stageXYc(f,2) + h/2 - 1 -pars.trimBorder);
        cs = round(stageXYc(f,1) - w/2 +pars.trimBorder : stageXYc(f,1) + w/2 - 1 -pars.trimBorder);    
        h1 = h - 2*pars.trimBorder;
        w1 = w - 2*pars.trimBorder;

        
       %  if pars.corrAlign || pars.blendEdges   % moved 
            % use overlap to fine tune position. 
                % expand the 'current image' just for the purposes of improved alignment
            % rs, cs - current image position
            pars.alignBuffer = 200; % number of pixels to pad current region for aligment
            c1 = max(cs(1)-pars.alignBuffer ,1);
            c2 = min(cs(end)+pars.alignBuffer, pars.nCols);
            r1 = max(rs(1)-pars.alignBuffer, 1);
            r2 = min(rs(end)+pars.alignBuffer, pars.nRows);
            
            rs2 = r1:r2;
            cs2 = c1:c2;
            currImage = mosaic(rs2,cs2); % mosaic(rs,cs);  
       
            % deal with edge effects
            [daxH,daxW] = size(dax);
            [curH,curW] = size(currImage);
            vPad = (curH-daxH)/2;
            hPad = (curW-daxW)/2;
            daxPad = padarray(dax,[vPad,hPad],nan);
            % extract just the non-nan portion of current frame
            %    (there has to be a more elegant way)
        if pars.corrAlign || pars.blendEdges     
            idData = ~isnan(currImage) & ~(currImage==0);
            [ys,xs] = ind2sub([length(rs2),length(cs2)],find(idData));
            y1 = max(min(ys)-pars.alignBuffer,1);
            y2 = min(max(ys)+pars.alignBuffer,length(rs2));
            x1 = max(min(xs)-pars.alignBuffer,1); % first row with data in current mosaic
            x2 = min(max(xs)+pars.alignBuffer,length(cs2)); % last row with data in current mosaic
            hasData = sum( idData(:) ) > 0; % min # overlapping pixels
            hasOverlap = sum( idData(:) ) > pars.minOverlap*h1*w1; % min # overlapping pixels 
        end
        if pars.corrAlign && hasOverlap  && isempty(pars.shiftsXY)
             im1 = IncreaseContrast(uint16(currImage(y1:y2,x1:x2)),'high',pars.corrAlignHigh,'low',pars.corrAlignLow);
             im2 = IncreaseContrast(uint16(daxPad(y1:y2,x1:x2)),'high',pars.corrAlignHigh,'low',pars.corrAlignLow);                
             if pars.showplots
                figure(10); clf;
                      im3=cat(3,im1,im2);
                      Ncolor(im3);   
                figure(11); clf;
                    subplot(1,2,1); imagesc(im1);
                    subplot(1,2,2); imagesc(im2); colormap('parula');
                
             end
             
            if ~isempty(pars.saveCorrFolder); figOut = figure(11); clf; end
            [xshift,yshift] = CorrAlign(im1,im2,'maxShift',pars.maxShift,...
                'minGrad',pars.minGrad,...
                'showplot', pars.showplots | saveFigure);  % the functional line
            if saveFigure
                SaveFigure(figOut,'name',['corrFig_',num2str(f,'%03d')],'formats',{'png','fig'},'overwrite',true);
            end
            
            shiftsXY(f,:) = [xshift,yshift];          
        elseif pars.corrAlign && ~isempty(pars.shiftsXY)
            if pars.veryverbose
                disp('using previously calculated shifts');
            end
            xshift = shiftsXY(f,1);
            yshift = shiftsXY(f,2);
        else
            xshift = 0;
            yshift = 0;
            if pars.veryverbose
                disp(['Pos. ',num2str(f),' insufficient overlap detected (or corrAlign is off)']);
                disp('troubleshoot');
                               
                im1 = IncreaseContrast(uint16(currImage),'high',pars.corrAlignHigh,'low',pars.corrAlignLow);
                im2 = IncreaseContrast(uint16(daxPad),'high',pars.corrAlignHigh,'low',pars.corrAlignLow);                
                figure(10); clf;
                      im3=cat(3,im1,im2);
                      Ncolor(im3);   
                figure(11); clf;
                    subplot(1,2,1); imagesc(im1);
                    subplot(1,2,2); imagesc(im2); colormap('parula');
                     
            end
            
        end
        if pars.corrAlign
           %-- (old version, Translate dax) --%
           % dax = TranslateImage(dax,xshift,yshift,'padValue',0,'trim',true);          
           % update stage position on mosaic rather than translate dax
           cs = cs + xshift;
           rs = rs + yshift; 
           cs(cs<1) = 1;
           rs(rs<1) = 1;            
           % stageXYc(f,:) = stageXYc(f,:) + [yshift,xshift];   % validated "+"        
           stageXYc(f,:) = stageXYc(f,:) + [xshift,yshift];   % validated "+"        
           if pars.veryverbose
               disp(['correcting drift ',num2str([xshift,yshift])]);
           end
        end 
        
        if pars.blendEdges && hasData
            % the right extreme of data in the mosaic image is: x2
            % the bottom extreme of data in the mosaic image is: y2
            mosDat = mosaic(rs,cs); 
            blend1 = GetBlendMask(mosDat); 
            blend2= 1-blend1;           
            im3 = nansum( cat(3, mosDat.*blend1, dax.*blend2 ), 3);
            
            if pars.showProgress
            figure(12); clf;
            subplot(2,2,1); imagesc(mosDat.*blend1);  colorbar;
            subplot(2,2,2); imagesc(dax.*blend2);  colorbar;
            subplot(2,2,3); imagesc(blend1); colorbar;
            subplot(2,2,4); imagesc(blend2);  colorbar;
            
            figure(13); clf; 
            subplot(1,3,1); imagesc(mosDat);
            subplot(1,3,2); imagesc(dax);
            subplot(1,3,3); imagesc(IncreaseContrast(uint16(im3),'high',.999));
            end
            
             mosaic(rs,cs) = im3;
        else
            currImage =  nanmedian( cat(3, mosaic(rs,cs), dax), 3);
            mosaic(rs,cs) = currImage;
        end
    
        if pars.showProgress
            figure(100); clf; 
             imagesc(IncreaseContrast(uint16(mosaic2),'high',.9999));
            hold on; text(cs(1),rs(1),num2str(f),'color','y'); title('mosaic2 (old)')
            
            figure(101); clf;
            imagesc(IncreaseContrast(uint16(mosaic),'high',.9999));
            hold on; text(cs(1),rs(1),num2str(f),'color','y'); title('mosaic (new)');
            
            mosaic2 = mosaic;
            pause(1);
        end
        if pars.verbose
            disp(f/numIms);
        end
        
    end
    mosaic = uint16(mosaic);
    
else
    % To troubleshoot, randomly assign different colors (R/G/B) to each
    % tile.  This can help check the quality of the pix_to_mm chosen. 
    mosaic = zeros(pars.nRows,pars.nCols,3,'uint16');
    for f=1:numIms   
        k=randi(3);
        rs = round(stageXYc(f,2) - h/2 : stageXYc(f,2) + h/2 - 1) ;
        cs = round(stageXYc(f,1) - w/2 : stageXYc(f,1) + w/2 - 1) ;
        dax = imData{f};
        if pars.transpose
            dax = dax';
        end
        if pars.fliplr
            dax = fliplr(dax);
        end
        if pars.flipud
            dax = flipud(dax);
        end
        mosaic(rs,cs,k) = dax;    
    end
    figure(); clf; imagesc(10*mosaic); colormap gray; title(sc);
end
% stageXY; % save original stage positions, listed relative to ORIGINAL image order
% shiftsXY; % save shifts, relative to the SORTED image order 

% resort back to original stacking order
[~,sortBack] = sort(sortIndex);
imData = imData(sortBack);
shiftsXY = shiftsXY(sortBack,:);
stageXYc = stageXYc(sortBack,:);

% everything exported will now be sorted relative to the original order
pars.stageXY = stageXY; % save original stage positions, listed relative to ORIGINAL image order
pars.shiftsXY = shiftsXY; % save shifts, 
pars.stageXY00 = stageXYc - repmat([w/2,h/2],numIms,1); % move position to upper-left of image. Listed relative to SORTED image order   


%%
% figure(2); clf; imagesc(IncreaseContrast(mosaic,'high',.9999)); colormap gray;

