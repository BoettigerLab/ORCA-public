function [dataOut,pars] = ChrTracer2_LoadData(expTableXLS,varargin)
% Load an expTableXLS table and use it to organize hybes into a fiducial
% stack and data stack.  Correct both stacks using image correlation of
% fiducial frames. 
% For speed only a single frame from each stack is loaded and aligned.


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
defaults(end+1,:) = {'regData','freeType',[]}; % optionally, pass a precomputed set of shifts to use. 
% other parameters
defaults(end+1,:) = {'defaultSaveDir','string','I:\Analyzed\'};
defaults(end+1,:) = {'midFrame','integer',50}; % if data exists, just load this frame  

pars = ParseVariableArguments(varargin, defaults, mfilename);
% pars = ParseVariableArguments([], defaults, mfilename);

%%
disp('Loading data...');

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
    folderName = strsplit(dataFolder,filesep);
    folderName = cat(2,folderName{3:end});
    saveFolder = SetFigureSavePath( [pars.defaultSaveDir,folderName,'_CT2out\'],'makeDir',true);
else
    saveFolder = SetFigureSavePath(pars.saveFolder);
end

% check if DataFolder is specified in addition to the standard Hybe Folder
% (this allows multiple parent directories to be included in the same
% experiment table). 
if sum(strcmp(eTable.Properties.VariableNames,'DataFolder'))
   dataFolders = eTable.DataFolder;
else
   dataFolders = cellstr(repmat(dataFolder,numHybes,1));
end

% Estimate number of channels based on first hybe (could be made more general later) 
% Identify fiducial channel and parameter channel
h = 1;
[fidChn,channels] = GetFidChn(eTable,h);
numChns = length(channels);
dataChns = 1:numChns; dataChns(fidChn) = [];
numDataChns = length(dataChns); 

    
% Test if data is already written
fidNames = cellstr(ls([saveFolder,'fov',num2str(pars.fov,'%03d'),'_h*fid.dax']));
datNames = cell(numHybes,numDataChns);
for n=1:numDataChns
    foundDat = cellstr(ls([saveFolder,'fov',num2str(pars.fov,'%03d'),'_h*dat',num2str(n),'.dax']));
    if isempty(foundDat{1})  % a little historical inconsistency problem 
        foundDat = cellstr(ls([saveFolder,'fov',num2str(pars.fov,'%03d'),'_h*data',num2str(n),'.dax']));
    end
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

fiducialFrames = cell(numHybes,1); 
goodHybes = true(1,numHybes);


if ~skip  
    % Step 1: Load data
    %   load in reference image for x,y pixel registration
    for h=1:numHybes 
        currFolder = [dataFolders{h},'\',hybFolders{h},'\'];
        daxFiles =  cellstr(ls([currFolder,pars.filename,'*.dax']));
        try
            fidChn = GetFidChn(eTable,h);
            fidFrame = (1+fidChn); 
            fiducialFrames{h}  = ReadDax([currFolder,daxFiles{pars.fov}],...
                'verbose',false,'startFrame',fidFrame,'endFrame',fidFrame);
        catch er
            goodHybes(h) = false;      
            if pars.stopOnError
                warning(er.message);
                error(['failed to load data from hybe ',num2str(h)]); 
            end
            if pars.verbose
                warning(er.message);
                disp(['failed to load data from hybe ',num2str(h)]); 
            end
        end
    end

    % Step 2: Compute registration
    % export the registration map and the result of registration. 
    [fiducialAlignFrames,regData,goodHybes] = RegisterImages_light(fiducialFrames,'parameters',pars,'goodHybes',goodHybes);
end

% If registered data already exists, just load that
%   we don't save the first in focus frame, just the z-stack, so now we
%   need to grab a frame half-way through the stack, 'midFrame', which
%   should be in focus.  
if skip
    fiducialAlignFrames = cell(numHybes,1); 
    for h=1:numHybes
        fiducialAlignFrames{h}  = ReadDax([saveFolder,fidNames{h}],...
                'verbose',false,'startFrame',pars.midFrame,'endFrame',pars.midFrame);
    end
end

dataOut.regData = regData;
dataOut.goodHybes = goodHybes;
dataOut.dataFolder = dataFolders;
dataOut.saveFolder = saveFolder;
dataOut.eTable = eTable;
dataOut.fiducialFrames = fiducialAlignFrames;  



function [fidChn,channels] = GetFidChn(eTable,h)
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
    [~,fidChn] = intersect(channels,fidChannel);
    if isempty(fidChn)
        error(['could not find fid channel ',fidChannel, ' in channel list ',channels{:}]);
    end
