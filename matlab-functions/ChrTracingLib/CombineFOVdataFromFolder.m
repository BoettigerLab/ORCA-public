function dataOut = CombineFOVdataFromFolder(dataFolder,varargin)
% dataOut = CombineFOVdataFromFolder(dataFolder);s


% default optional parameters
defaults = cell(0,3);
% internal pars
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'selectReads','integer',0};
defaults(end+1,:) = {'filterReads','fraction',0};
defaults(end+1,:) = {'filterFids','fraction',0};
defaults(end+1,:) = {'datTables','freeType',{}};
defaults(end+1,:) = {'mosaicPars','freeType',[]};
defaults(end+1,:) = {'selectFOV','freeType',[]};
defaults(end+1,:) = {'subFolders','boolean',true};
% TableToDistMap pars
defaults(end+1,:) = {'bins','integer',0}; % fixed matrix size rather than auto determine 
defaults(end+1,:) = {'impute','boolean',false}; % don't use
defaults(end+1,:) = {'removeEmpty','boolean',true};
defaults(end+1,:) = {'byHybe','boolean',true};
defaults(end+1,:) = {'dims','freeType',3};
defaults(end+1,:) = {'tableFormat',{'new','old','chrom'},'new'};
defaults(end+1,:) = {'chromCorrect','boolean',true};
defaults(end+1,:) = {'sortMethod',{'byPanel','byRead'},'byPanel'};
defaults(end+1,:) = {'dataBasedDriftCorrect','boolean',false};
% RemoveBadPts pars
defaults(end+1,:) = {'maxJump','nonnegative',250};
defaults(end+1,:) = {'maxDist','nonnegative',2000};
% for mosaic
defaults(end+1,:) = {'transpose','boolean',true}; % see the mosaic parameters in the Hal parameter file used to collect the movie 
defaults(end+1,:) = {'fliplr','boolean',true};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'flipud','boolean',false};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'trimBorder','nonnegative',0};
defaults(end+1,:) = {'reZero','array',[0,0]};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% find files 
if isempty(pars.datTables)
    if pars.subFolders
        datTables = FindFileInSubFolders(dataFolder,'*_AllFits.csv');
    else
        datTables =strcat(dataFolder, cellstr(ls( [dataFolder,'*_AllFits.csv']) ));
    end
    if isempty(datTables)
        datTables =strcat(dataFolder, cellstr(ls( [dataFolder,'*_AllFits.csv']) ));
    end
    if isempty(datTables)
       error(['no AllFits.csv found in ',dataFolder]); 
    end
else
    datTables = pars.datTables;
end
    
nFOV = length(datTables);
imTables = cell(nFOV,1);
distMaps = cell(nFOV,1);
distMapFilts = cell(nFOV,1);
polymers = cell(nFOV,1);
polymerFilts = cell(nFOV,1);
polymerFill = cell(nFOV,1);
distMapFill = cell(nFOV,1);
spotProps = cell(nFOV,1);
if pars.selectReads 
    nReads = length(pars.selectReads);
    rs = pars.selectReads;
else
    nReads = 0;
end

if isempty(pars.selectFOV)
    selectFOV = 1:nFOV;
else
    selectFOV = pars.selectFOV;
end
goodTables = true(1,nFOV);

for f=selectFOV % f=5
    if ischar(datTables{f})
        imTables{f} = readtable(datTables{f});
        % Fix buggy datasets that don't have the FOV recorded
        [~,tabName] = fileparts(datTables{f});
        fovNum = str2double(tabName(4:6));
        tabFOV = imTables{f}.fov(1);
        if tabFOV~=fovNum
            newTab = imTables{f};
            newTab.fov = fovNum*ones(height(newTab),1);
            newTab.fs = CantorPair(newTab.fov,newTab.s);
            newTab.fsr = CantorPair(newTab.fs,newTab.readout);
            imTables{f} = newTab;
            writetable(newTab,datTables{f});
            display(['updated fov number in table ',datTables{f}]);
        end
        
    else
        imTables{f} = datTables{f};
    end
    if isempty(imTables{f})
       goodTables(f) = false; 
    else
        
        % Filter the data based on brightness before further analysis
        if pars.filterFids % remove x% dimest fiducials.
         cut = quantile(imTables{f}.fid_h,pars.filterFids);
         sel =  imTables{f}.fid_h < cut;
         imTables{f}(sel,:) = [];
        end
        if pars.filterReads % remove x% dimest reads among all spots that had that read
            for r=1:nReads
                sel = imTables{f}.hybe==r;
                cut = quantile(imTables{f}.h(sel),pars.filterReads);
                if isnan(cut) 
                    cut = 0;
                end
                sel = imTables{f}.hybe==r & imTables{f}.h < cut;
                imTables{f}(sel,:) = [];
            end
        end
        
        % Convert table to distance map and polymer structure
        [distMaps{f},polymers{f},spotProps{f}] = TableToDistMap(imTables{f},...
              'parameters',pars);
%             'byHybe', pars.byHybe,'removeEmpty',pars.removeEmpty,...
%             'dims',pars.dims);      
        
        
        % Pad or Trim Distance Maps to make select reads if requested
        nHybes = size(distMaps{f},1);
        if nReads < nHybes && pars.selectReads(end)
            distMapFilt = distMaps{f}(rs,rs,:);
            polymerFilt = polymers{f}(rs,:,:);
        elseif nReads > nHybes && pars.selectReads(end)
            distMapFilt = padarray( distMaps{f}, [nReads-nHybes,nReads-nHybes],nan,'post');
            polymerFilt =padarray( polymers{f}, [nReads-nHybes,0],nan,'post');
        else
            distMapFilt = distMaps{f};
            polymerFilt = polymers{f};
        end
        % remove bad points from distance map
        [nReads,~,nCells] = size(distMapFilt);
        [badPts,strays] = RemoveBadPts(distMapFilt,...
            'maxJump',pars.maxJump,...
            'maxDist',pars.maxDist,...  
            'verbose',pars.verbose); % 300  
        distMapFilt(badPts) = NaN;
        distMapFilts{f} = distMapFilt;

        % remove bad points from polymers                 
        stray = permute(repmat(strays,1,1,4),[1,3,2]);
        polymerFilt(stray) = NaN;
        polymerFilts{f} = polymerFilt;
        
        % interpolate data for smoother plots
        polymerFill{f} = zeros(nReads,3,nCells);
        distMapFill{f} = nan(nReads,nReads,nCells);
        for n=1:nCells
            polyTemp = fillmissing(polymerFilt(1:nReads,1:3,n),'linear');
            polymerFill{f}(:,:,n) =polyTemp;
            distMapFill{f}(:,:,n) = squareform(pdist( polyTemp ) );
        end   
        
        % add unique spotIDs into the spot table
        imTables{f}.uniqueID = imTables{f}.s + imTables{f}.fov*1E5;

        % start a new table of mosaic positions
        if ~isempty(pars.mosaicPars)
            locusX = spotProps{f}.locusX;
            locusY = spotProps{f}.locusY;
            locusX(locusX==0) = NaN;
            locusY(locusY==0) = NaN;
            if pars.trimBorder > 0
                locusX = locusX - pars.trimBorder;
                locusY = locusY - pars.trimBorder;
            end
            if pars.transpose
                tempX = locusX;
                tempY = locusY;
                locusX = tempY;
                locusY = tempX;
            end
            if pars.fliplr
                locusX = pars.mosaicPars.imWidth - locusX+1;  
            end
            if pars.flipud
                locusY = pars.mosaicPars.imHeight - locusY+1;  % 
            end
            spotProps{f}.mosaicX = locusX + pars.mosaicPars.stageXY00(f,1) - pars.reZero(1);
            spotProps{f}.mosaicY = locusY + pars.mosaicPars.stageXY00(f,2) - pars.reZero(2);
%            figure(3); clf; plot(spotProps{f}.mosaicX,spotProps{f}.mosaicY,'k.');
%            spotProps{f}.mosaicX = (pars.mosaicPars.imHeight - spotProps{f}.locusY) + pars.mosaicPars.stageXY00(f,1) ;
%            spotProps{f}.mosaicY =  spotProps{f}.locusX + pars.mosaicPars.stageXY00(f,2) ;
%            figure(3); clf; plot(spotProps{f}.mosaicX,spotProps{f}.mosaicY,'k.');
%            spotProps{f}.fov = f*ones(height(spotProps{f}),1);
%            spotProps{f}.s = floor(1:.5:height(spotProps{f})/2+.5)';
        end
    end
end

try
    % flatten tables
    dataOut.distMapRaw = distMaps; 
    dataOut.polymerRaw = polymers;
    nReads = cellfun(@(x) size(x,1), distMapFilts);
    allReads = nReads == max(nReads);
    if any(allReads==0) && pars.verbose
        warning('Variable number of reads found. Some data may be exculded from the filtered distance map and polymer matrices');
    end
    dataOut.distMapFilt = cat(3,distMapFilts{allReads});
    dataOut.distMapFill = cat(3,distMapFill{allReads});
    dataOut.polymerFilt = cat(3,polymerFilts{allReads});
    dataOut.polymerFill = cat(3,polymerFill{allReads});
    dataOut.datTable = cat(1,imTables{goodTables});
    spotStats = cat(1,spotProps{goodTables});
    dataOut.spotProps  = spotStats;
catch er
    disp(er.getReport);
    stopNow = true;
    disp(stopNow);
    error(er.getReport);
end

if pars.verbose
    disp(['finished loading all data from folder ',dataFolder]);
end