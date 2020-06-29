function dataOut = CombineFOVdataFromFolder(dataFolder,varargin)
% dataOut = CombineFOVdataFromFolder(dataFolder);


% default optional parameters
defaults = cell(0,3);
defaults(end+1,:) = {'maxJump','nonnegative',250};
defaults(end+1,:) = {'maxDist','nonnegative',2000};
defaults(end+1,:) = {'byHybe','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'selectReads','integer',0};
defaults(end+1,:) = {'filterReads','fraction',0};
defaults(end+1,:) = {'filterFids','fraction',0};
defaults(end+1,:) = {'datTables','freeType',[]};
defaults(end+1,:) = {'removeEmpty','boolean',true};
defaults(end+1,:) = {'mosaicPars','freeType',[]};
defaults(end+1,:) = {'selectFOV','freeType',[]};
defaults(end+1,:) = {'dims','integer',3};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% find files 
if isempty(pars.datTables)
    datTables = FindFileInSubFolders(dataFolder,'*_AllFits.csv');
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
totalSpots = 0;
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

for f=selectFOV % f=5
    imTables{f} = readtable(datTables{f});
    if ~isempty(imTables{f})
        
        % Filter the data based on brightness before further analysis
        if pars.filterFids % remove x% dimest fiducials. 
         sel = imTables{f}.hybe==1;
         cut = quantile(imTables{f}.fid_h(sel),pars.filterFids);
         sel = imTables{f}.hybe==1 & imTables{f}.fid_h < cut;
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
            'byHybe', pars.byHybe,'removeEmpty',pars.removeEmpty,...
            'dims',pars.dims);      
        
        
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
        [badPts,strays] = RemoveBadPts(distMapFilt,'maxJump',pars.maxJump,'maxDist',pars.maxDist,'verbose',pars.verbose); % 300  
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
        imTables{f}.uniqueID = imTables{f}.s + totalSpots;
        totalSpots = totalSpots + 2*max(imTables{f}.s);

        % start a new table of mosaic positions
        if ~isempty(pars.mosaicPars)
            spotProps{f}.mosaicX = (pars.mosaicPars.imHeight - spotProps{f}.locusY) + pars.mosaicPars.stageXY00(f,1) ;
            spotProps{f}.mosaicY =  spotProps{f}.locusX + pars.mosaicPars.stageXY00(f,2) ;
            spotProps{f}.fov = f*ones(height(spotProps{f}),1);
            spotProps{f}.s = floor(1:.5:height(spotProps{f})/2+.5)';
        end
    end
end

try
    % flatten tables
    dataOut.distMapRaw = distMaps; 
    dataOut.polymerRaw = polymers;
    dataOut.distMapFilt = cat(3,distMapFilts{:});
    dataOut.distMapFill = cat(3,distMapFill{:});
    dataOut.polymerFilt = cat(3,polymerFilts{:});
    dataOut.polymerFill = cat(3,polymerFill{:});
    dataOut.datTable = cat(1,imTables{:});
    spotStats = cat(1,spotProps{:});
    spotStats.uniqueID = (1:height(spotStats))';
    dataOut.spotProps  = spotStats;
catch er
    error(er.getReport);
end

if pars.verbose
    disp(['finished loading all data from folder ',dataFolder]);
end