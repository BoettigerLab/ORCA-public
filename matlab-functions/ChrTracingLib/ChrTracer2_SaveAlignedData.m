function [fidMapAlign,datMapAlign,fidFlatFrames,datFlatFrames] = ChrTracer2_SaveAlignedData(eTable,regData,varargin)

defaults = cell(0,3);
% general FOV parameters
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'saveData', 'boolean', false}; 
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'dataFolder','string',''};
defaults(end+1,:) = {'refHybe','integer',1};
% general pars
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'filename', 'string', 'ConvZscan'}; 
defaults(end+1,:) = {'overwrite', 'boolean', 'false'}; 
defaults(end+1,:) = {'memMapData','boolean',true};
defaults(end+1,:) = {'loadFlatData','boolean',true};
defaults(end+1,:) = {'loadFlatFid','boolean',true};
defaults(end+1,:) = {'useTableName','boolean',false}; % to be made true in future release  
defaults(end+1,:) = {'saveMaxProject','boolean',true}; 
defaults(end+1,:) = {'warpMaps','cell',{}}; 
defaults(end+1,:) = {'warpChns','integer',0}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

tic;
% ------------- Extract hyb folders ------------------
hybFolders = eTable.FolderName;
numHybes = length(hybFolders); 

% ------------- Extract or create DataFolders Array --------------
% check if DataFolder was specified in addition to the standard HybeFolder
% (this allows multiple parent directories to be included in the same
% experiment table). 
dataFolder = pars.dataFolder;
if sum(strcmp(eTable.Properties.VariableNames,'DataFolder'))
   dataFolders = eTable.DataFolder;
else
   dataFolders = cellstr(repmat(dataFolder,numHybes,1));
end


% ------ Identify fiducial channel and parameter channel ------
% Determine number of channels based on first hybe 
h = 1; % (could be made more general later) 
try
    channels = strsplit(regexprep(eTable.channels{h},'[^A-Za-z0-9,]',''),',');
catch
    error('channel should be formated as a string, add braces or a leading quote');
end
if ~iscell(eTable.fiducialChannel(h))
    fidChannel = num2str(eTable.fiducialChannel(h));
else
    fidChannel = regexprep(eTable.fiducialChannel{h},'[\[,\]]','');
end 
numChns = length(channels);
[~,fidChn] = intersect(channels,fidChannel);
dataChns = 1:numChns; 
dataChns(fidChn) = [];
numDataChns = length(dataChns); 
if isempty(fidChn)
    error(['could not find fid channel ',fidChannel, ' in channel list ',channels{:}]);
end

%------------------- File names ---------------------------------------
if pars.useTableName
    nameField = find(contains(eTable.Properties.VariableNames,'name','IgnoreCase',true));
    if length(nameField)>1
        nameField = nameField(end);
    end
    hybeNames = regexprep(eTable{:,nameField},{'_',' '},{'',''});
else
    hybeNames = cellstr(strcat('h',num2str( (1:numHybes)','%03d')));
end

% ------------------- Test if data is already written ------------------
fidNames = cellstr(ls([pars.saveFolder,'fov',num2str(pars.fov,'%03d'),'*fid.dax']));
datNames = cell(numHybes,numDataChns);
for n=1:numDataChns
    foundDat = cellstr(ls([pars.saveFolder,'fov',num2str(pars.fov,'%03d'),'_*data',num2str(n),'.dax']));
    if isempty(foundDat{1}) % old naming convention, will remove in future release 
        foundDat = cellstr(ls([pars.saveFolder,'fov',num2str(pars.fov,'%03d'),'_*dat',num2str(n),'.dax']));
    end
    datNames(1:length(foundDat),n) = foundDat;
end
foundAllFid = length(fidNames) == numHybes;
foundAllDat = sum(cellfun(@isempty,datNames(:))) == 0;
if foundAllFid && foundAllDat && ~pars.overwrite
    skip = true;
    cprintf([1 .25 0],'found existing aligned dax files, skipping alignment');
elseif  foundAllFid && foundAllDat && pars.overwrite
    skip = false;
    cprintf([1 .25 0],'overwriting existing aligned dax files');
else
    skip = false;
end
fidNames = strcat(pars.saveFolder,fidNames);
for n=1:numDataChns
    datNames(:,n) = strcat(pars.saveFolder,datNames(:,n));
end

% load chromatic warp maps if requested
if ~isempty(pars.warpMaps)
    nWarps = length(pars.warpMaps);
    warpMaps = cell(nWarps,1);
   for n=1:nWarps
       load(pars,warpMaps,'warpMap');
       warpMaps{n} = warpMap;
   end
else
    warpMaps = [];
end


fidFlatFrames = cell(numHybes,1);
datFlatFrames = cell(numHybes,numDataChns); 
%------------------ Save new aligned dax file ---------------------------% 
if ~skip
    datNames = cell(numHybes,numDataChns);
    fidNames = cell(numHybes,1);
    % loop through all hybes, but do refHybe first
    allHybes = 1:numHybes;
    allHybes(pars.refHybe) = []; 
    allHybes = [pars.refHybe,allHybes];
    for h=allHybes
        if pars.verbose
            disp(['loading data from hybe ',num2str(h), '  folder: ',hybFolders{h}]);
        end
        % check if file already exists
        fidInfName = ['fov',num2str(pars.fov,'%03d'),'_',hybeNames{h},'_fid','.inf'];
        if ( exist([pars.saveFolder,fidInfName],'file')==0 || pars.overwrite ) || h==pars.refHybe

            currFolder = [dataFolders{h},filesep,hybFolders{h},filesep];
            daxFiles =  cellstr(ls([currFolder,pars.filename,'*.dax']));
            % Read in data for this hybe    
            [dax,info] = ReadDax([currFolder,daxFiles{pars.fov}],'verbose',false);
            info.localPath = pars.saveFolder;
            % Stage_X/Stage_Y varies between hybes for unknown reasons
            %   since we have corrected positions to subpixel accuracy we want
            %   to record the corrected positions. 
            if h==pars.refHybe
               refStageX = info.Stage_X;
               refStageY = info.Stage_Y;
            else
                info.notes = ['Orig Stage_XY=(',num2str(info.Stage_X),',',num2str(info.Stage_Y),')'];
                info.Stage_X = refStageX;
                info.Stage_Y = refStageY;       
            end

            % split out data frames and fiducial frames 
            keepFrames = str2num(regexprep(eTable.framesToKeep{h},'[\[,\]]','')); %#ok<ST2NM>
            fidFrames = (1+fidChn):numChns:max(keepFrames);
            fiducialFrames = dax(:,:,intersect(fidFrames,keepFrames));
            dataFrames = cell(numDataChns,1);
            for n=1:numDataChns
                dataKeepFrames = (1+dataChns(n)):numChns:max(keepFrames);
                dataFrames{n} = dax(:,:,intersect(dataKeepFrames,keepFrames)); % fix me
            end

            % Align and save with correct infoFile and updated stage positions
            fidReg = ApplyReg(fiducialFrames,regData(h));
            fidInfo = info;
            fidInfo.number_of_frames = size(fidReg,3); 
            fidInfo.localName = fidInfName;
            WriteDAXFiles(fidReg,fidInfo,'verbose',pars.verbose);
            fidNames{h} = [pars.saveFolder,fidInfo.localName(1:end-4),'.dax'];
            fidFlatFrames{h} = max(fidReg,[],3); % record max project
            for n = 1:numDataChns
                datReg = ApplyReg(dataFrames{n},regData(h));
                datInfo = info;
                % Apply chromatic correction if requested
                if ~isempty(warpMaps)
                    datReg = imwarp(datReg,warpMaps{n},'OutputView',imref3d(size(datReg))); 
                    datInfo.notes = ['chromatic warp: ' pars.warpMaps{n}];
                end
                datInfo.number_of_frames = size(fidReg,3);
                datInfo.localName = ['fov',num2str(pars.fov,'%03d'),'_',hybeNames{h},'_data',num2str(n),'.inf']; 
                WriteDAXFiles(datReg,datInfo,'verbose',pars.verbose);
                datNames{h,n} = [pars.saveFolder,datInfo.localName(1:end-4),'.dax'];
                datFlatFrames{h,n} = max(datReg,[],3);
            end  
        else
            fidNames{h} = [pars.saveFolder,fidInfName(1:end-4),'.dax'];
            for n = 1:numDataChns
                localName = ['fov',num2str(pars.fov,'%03d'),'_',hybeNames{h},'_data',num2str(n),'.inf']; 
                datNames{h,n} = [pars.saveFolder,localName(1:end-4),'.dax'];
            end  
            if pars.verbose
                disp(['found existing ', fidInfName ,'.  skipping to next image']);
            end
        end
    end
end

%------------------------------------------------------------------------% 


% ------------ Create new memory maps of saved AlignedDax -----------------

fidMapAlign = cell(numHybes,1);
datMapAlign = cell(numHybes,numDataChns);
if pars.memMapData
    if pars.verbose
        disp('memory mapping data');
    end
    for h=1:numHybes
        [dax,fidMapAlign{h}] = ReadData(fidNames{h},'loadFile',pars.loadFlatFid);
        fidFlatFrames{h} = max(dax,[],3);
        infoOut = fidMapAlign{h}.infoFile;
        maxName = ['max_',infoOut.localName];
        if pars.saveMaxProject && ~exist([pars.saveFolder,filesep,maxName],'file') % write max projection for speed.
            infoOut.number_of_frames= 1;
            infoOut.localName = maxName;
            WriteDAXFiles(fidFlatFrames{h},infoOut,'verbose',pars.verbose);
        end
        for n=1:numDataChns
            [dax,datMapAlign{h,n}] = ReadData(datNames{h,n},'loadFile',pars.loadFlatData);
            datFlatFrames{h,n} = max(dax,[],3);
            infoOut = datMapAlign{h}.infoFile;
            maxName = ['max_',infoOut.localName];
            if pars.saveMaxProject && ~exist([pars.saveFolder,filesep,maxName],'file') % write max projection for speed.
                infoOut.number_of_frames= 1;
                infoOut.localName = maxName;
                WriteDAXFiles(datFlatFrames{h},infoOut,'verbose',pars.verbose);
            end
        end
    end
    fidFlatFrames = cat(3,fidFlatFrames{:});
    datFlatFrames = cat(3,datFlatFrames{:});
end
% % data returned by function
% dataOut.fidMapData = fidMapAlign;
% dataOut.datMapData  = datMapAlign;


t = toc;
disp(['data loaded in ',num2str(t/60),' min']);

%%
function [movie,mapData] = ReadData(fileName,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'startFrame','integer',1}; % 0 = first frame
defaults(end+1,:) = {'endFrame','integer',0}; % 0 = last frames
defaults(end+1,:) = {'selectFrames','integer',0}; % 0 = all frames
defaults(end+1,:) = {'roi','integer',0}; % 0 = full field
defaults(end+1,:) = {'mapName','string',''};
defaults(end+1,:) = {'loadFile','boolean',true};
pars = ParseVariableArguments(varargin,defaults,'ReadData');
% pars = ParseVariableArguments([],defaults,'ReadData');

% parse info file 
infoFile = ReadInfoFile(fileName, 'verbose', pars.verbose);
framesInDax = infoFile.number_of_frames;
frameDim = [infoFile.frame_dimensions(1)/infoFile.binning(1),...
            infoFile.frame_dimensions(2)/infoFile.binning(2)];
frameDim = uint32(frameDim);
framesInDax = uint32(framesInDax);
if ~isempty( strfind(infoFile.data_type,'little endian') )
    binaryFormat = 'l';
else
    binaryFormat = 'b';
end

if pars.roi == 0
   xi = 1; xe = frameDim(1);
   yi = 1; ye = frameDim(2);
else
    xi = pars.roi(1);
    xe = pars.roi(2);
    yi = pars.roi(3);
    ye = pars.roi(4);
end
if pars.startFrame == 0
    startFrame = 1;
else
    startFrame = pars.startFrame;
end
if pars.endFrame == 0
    endFrame = framesInDax;
else
    endFrame = pars.endFrame;
end
if sum(pars.selectFrames) == 0
    selectFrames = uint32(startFrame):uint32(endFrame);
else
    selectFrames = uint32(pars.selectFrames);
end
framesToLoad = length(selectFrames);


% main function
memoryMap = memmapfile(fileName, ...
            'Format', 'uint16', ...
            'Writable', false, ...
            'Offset', 0); %  (startFrame-1)*frameSize*16/8);  

if pars.loadFile
    [ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);
    inds = sub2indFast([frameDim(2),frameDim(1),framesInDax],ri(:),ci(:),zi(:));
    movie = memoryMap.Data(inds);
    if strcmp(binaryFormat,'b')
        movie = swapbytes(movie);
    end
    xs = xe-xi+uint32(1);
    ys = ye-yi+uint32(1);
    movie = reshape(movie,[ys,xs,framesToLoad]); % figure(7); clf; imagesc(max(movie,[],3)); colormap gray;
else
    movie = [];
end 
mapData.memMap = memoryMap;
mapData.frameDim =frameDim;
mapData.framesInDax = framesInDax;
mapData.selectFrames = selectFrames;
mapData.name = pars.mapName;   
mapData.binaryFormat = binaryFormat;
mapData.infoFile = infoFile;
