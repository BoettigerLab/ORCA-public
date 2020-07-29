function [dataOut,pars] = ChrTracer_LoadData_light(expTableXLS,varargin)
% Load an expTableXLS table and use it to organize hybes into a fiducial
% stack and data stack.  Correct both stacks using image correlation of
% fiducial frames. 


% dataFolder = 'E:\Nasa\2017-03-23_K562-chr21_combined';

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
global scratchPath; 

defaults = cell(0,3);
% key parameters
% general FOV parameters
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryVerbose', 'boolean', false}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'showExtraPlots', 'boolean', false}; 
defaults(end+1,:) = {'saveData', 'boolean', false}; 
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'saveFolder','string',''};
% load data pars
defaults(end+1,:) = {'filename', 'string', 'ConvZscan'}; 
defaults(end+1,:) = {'dataTypes','cell',{'Hyb','Repeat','Toehold','H','T','R','any'}};
defaults(end+1,:) = {'overwrite','boolean',false}; % realign data even if existing files
% registration parameters
defaults(end+1,:) = {'rotation','boolean',false}; % also compute and correct rotation angle
defaults(end+1,:) = {'maxD','positive',15}; % distance in pixels to match objects prior to computing rotation
defaults(end+1,:) = {'alignmentBoxWidth', 'positive', inf}; % pixels.  size of region to use for frame alignment
defaults(end+1,:) = {'alignUpsample', 'positive', 1}; % upsample data for coarse alignment (slow, should be unnecessary);
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .999}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'previousAlignFrames','array',[]}; % previous align frames
defaults(end+1,:) = {'targetHybes','integer',[]};
defaults(end+1,:) = {'goodHybes','boolean',[]};
defaults(end+1,:) = {'corrAngles','freeType',[]}; % rotate angles to check for better alignment

pars = ParseVariableArguments(varargin, defaults, mfilename);
% pars = ParseVariableArguments([], defaults, mfilename);

%%
disp('Loading data...');
tic

% Step 1: Match file names
% load experiment table to know what probes are in which folders and
% whether they contain hybe data, repeat-hybe data, or toehold data
dataFolder = fileparts(expTableXLS); % get folder name;
eTable = readtable(expTableXLS);
hybFolderID = false(height(eTable),1);
if sum(  strcmp(pars.dataTypes,'any') ) == 0
    for i=1:length(pars.dataTypes)
        hybFolderID = hybFolderID | strcmp(eTable.DataType,pars.dataTypes{i});
    end
else
    hybFolderID = true(height(eTable),1);
end
numHybes = sum(hybFolderID);
eTable = eTable(hybFolderID,:);
hybFolders = eTable.FolderName;

% create a save folder
if isempty(pars.saveFolder) || strcmp(pars.saveFolder,scratchPath)
    saveFolder = SetFigureSavePath( [dataFolder,'_CT2out\'],'makeDir',true);
else
    saveFolder = SetFigureSavePath(pars.saveFolder);
end

% setup dataFolder
if sum(strcmp(eTable.Properties.VariableNames,'DataFolder'))
   dataFolder = eTable.DataFolder;
else
   dataFolder = cellstr(repmat(dataFolder,numHybes,1));
end

% Estimate number of channels based on first hybe (could be made more general later) 
% Identify fiducial channel and parameter channel
h = 1;
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
dataChns = 1:numChns; dataChns(fidChn) = [];
numDataChns = length(dataChns); 
if isempty(fidChn)
    error(['could not find fid channel ',fidChannel, ' in channel list ',channels{:}]);
end

% we need
% an (x,y)-shift for each hybe (computed from fiducial)
% memory map for each point 

fidNames = cellstr(ls([saveFolder,'fov',num2str(pars.fov,'%03d'),'_h*fid.dax']));
datNames = cell(numHybes,numDataChns);
for n=1:numDataChns
    foundDat = cellstr(ls([saveFolder,'fov',num2str(pars.fov,'%03d'),'_h*dat',num2str(n),'.dax']));
    datNames(1:length(foundDat),n) = foundDat;
end
foundAllFid = length(fidNames) == numHybes;
foundAllDat = sum(cellfun(@isempty,datNames(:))) == 0;

if foundAllFid && foundAllDat && ~pars.overwrite
    skip = true;
    cprintf([1 .25 0],'found existing aligned dax files, skipping alignment');
    regData = [];
    fiducialAlignFrames = [];
elseif  foundAllFid && foundAllDat && pars.overwrite
    skip = false;
    cprintf([1 .25 0],'overwriting existing aligned dax files');
else
    skip = false;
end
    


% Load Data
fidMapData = cell(numHybes,1);
datMapData = cell(numHybes,numDataChns);
fiducialFrames = cell(numHybes,1); 
goodHybes = true(1,numHybes);


if ~skip
    
    
for h=1:numHybes % load data
    if pars.verbose
        disp(['loading data from hybe ',num2str(h), '  folder: ',hybFolders{h}]);
    end
    currFolder = [dataFolder{h},'\',hybFolders{h},'\'];
    daxFiles =  cellstr(ls([currFolder,pars.filename,'*.dax']));
    try
        fidName = [currFolder,regexprep(daxFiles{pars.fov},'.dax',''),'_fid'];
        keepFrames = str2num(regexprep(eTable.framesToKeep{h},'[\[,\]]','')); %#ok<ST2NM>
        fidFrames = (1+fidChn):numChns:max(keepFrames);
        keepFidFrames = intersect(fidFrames,keepFrames);
        [dax,fidMapData{h}] = ReadData([currFolder,daxFiles{pars.fov}],...
            'verbose',false,'selectFrames',keepFidFrames,'mapName',fidName);
        
        fiducialFrames{h} = max(dax,[],3); % save some memory
        for n=1:numDataChns
            dataKeepFrames = (1+dataChns(n)):numChns:max(keepFrames);
            datFrames = intersect(dataKeepFrames,keepFrames);
            datName = [currFolder,regexprep(daxFiles{pars.fov},'.dax',''),'_data',num2str(n)];
            [~,datMapData{h,n}] = ReadData([currFolder,daxFiles{pars.fov}],...
                'verbose',false,'selectFrames',datFrames,'mapName',datName);
        end
    catch er
        goodHybes(h) = false;
        
            warning(er.message);
        if pars.stopOnError
            warning(er.message);
            error(['failed to load data from hybe ',num2str(h)]); 
        end
        if pars.verbose
            disp(['failed to load data from hybe ',num2str(h)]); 
        end
    end
end

% Align Data
[fiducialAlignFrames,regData,goodHybes] = RegisterImages_light(fiducialFrames,'parameters',pars,'goodHybes',goodHybes);

% Save AlignedDax
fidNames = SaveAlignedDax(fidMapData,regData,'tag','fid','saveFolder',saveFolder,'fov',pars.fov);
datNames = cell(numHybes,numDataChns);
for n=1:numDataChns
    datNames(:,n) = SaveAlignedDax(datMapData(:,n),regData,'tag',['dat',num2str(n)],'saveFolder',saveFolder,'fov',pars.fov);
end

end

% Map AlignedDax
fidMapAlign = cell(numHybes,1);
datMapAlign = cell(numHybes,numDataChns);
for h=1:numHybes
    [~,fidMapAlign{h}] = ReadData([saveFolder,fidNames{h}],'loadFile',false);
    for n=1:numDataChns
        [~,datMapAlign{h,n}] = ReadData([saveFolder,datNames{h,n}],'loadFile',false);
    end
end
if skip
    temp = ReadFromMemMap(fidMapAlign{1});
    fiducialAlignFrames = {max(temp,[],3)};
end
% allFrames = cat(3,fiducialFrames{:});
% im = max(allFrames,[],3);
% figure(1); clf; imagesc(30*im);

% data returned by function
dataOut.fidMapData = fidMapAlign;
dataOut.datMapData  = datMapAlign;
dataOut.regData = regData;
dataOut.goodHybes = goodHybes;
dataOut.dataFolder = dataFolder;
dataOut.saveFolder = saveFolder;
dataOut.eTable = eTable;
dataOut.fiducialFrames = fiducialAlignFrames;  % (drop this?) 

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
