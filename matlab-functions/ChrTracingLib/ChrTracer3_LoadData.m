function [dataOut,pars] = ChrTracer3_LoadData(expTableXLS,varargin)
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
defaults(end+1,:) = {'alignMaxProject','boolean',false}; % use maxprojection instead of mid-frame 
% registration parameters
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .999}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'regData','freeType',[]}; % optionally, pass a precomputed set of shifts to use. 
% CorrAlignFast Parameters
defaults(end+1,:) = {'maxSize', 'positive', 400}; % rescale all images to this size for alignment
defaults(end+1,:) = {'fineBox', 'freeType', []};  % perform fine scale alignment using a box of this size around the brightest point.
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'minGrad', 'float', -inf};
defaults(end+1,:) = {'angles','float',0}; % -10:1:10
defaults(end+1,:) = {'scales','float',1}; % -10:1:10
defaults(end+1,:) = {'fineMaxShift', 'nonnegative', 10};
defaults(end+1,:) = {'fineAngles','float',0}; % -1:.1:1
defaults(end+1,:) = {'fineScales','float',1}; % 0.95:0.01:.1.05
defaults(end+1,:) = {'fineCenter','array',[0,0]};
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'showplot', 'boolean', true};
defaults(end+1,:) = {'fastDisplay', 'boolean', true};
defaults(end+1,:) = {'displayWidth', 'integer', 500};
defaults(end+1,:) = {'showExtraPlot', 'boolean', false};
defaults(end+1,:) = {'minFineImprovement', 'float', .5}; 
defaults(end+1,:) = {'showCorrAlign', 'boolean', false};

% other parameters
defaults(end+1,:) = {'defaultSaveDir','string',''};

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
    if ~isempty(pars.defaultSaveDir)
        folderName = strsplit(dataFolder,filesep);
        folderName = cat(2,folderName{3:end});
        saveFolder = SetFigureSavePath( [pars.defaultSaveDir,folderName,'_CT2out\'],'makeDir',true);
    else
        saveFolder = SetFigureSavePath( [dataFolder,'_CT2out\'],'makeDir',true);
    end
else
    saveFolder = SetFigureSavePath(pars.saveFolder);
end

% check if DataFolder is specified in addition to the standard Hybe Folder
% (this allows multiple parent directories to be included in the same
% experiment table). 
dataInTableFolder = false;
if sum(strcmp(eTable.Properties.VariableNames,'DataFolder'))
    if ~isempty(eTable.DataFolder(1))
        dataFolders = eTable.DataFolder;
    else
        dataInTableFolder = true;
    end
else
    dataInTableFolder = true;
end
if dataInTableFolder
   dataFolders = cellstr(repmat(dataFolder,numHybes,1));
end


% convert annotations of all frames
[isFidChn,frameChannels,channels,midFrame] = GetFidChnFromTable(eTable);
numChns = length(channels);
numDataChns = numChns -1; 

    
% Test if data is already written
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

fiducialFrames = cell(numHybes,1); 
if ~skip  
    % Step 1: Load data
    %   load in single frame reference image for x,y pixel registration
    for h=1:numHybes 
        currFolder = [dataFolders{h},'\',hybFolders{h},'\'];
        daxFiles =  cellstr(ls([currFolder,pars.filename,'*.dax']));
        try
            % isFidChn = GetFidChnFromTable(eTable); % Currently don't allow this to vary as a function of hybe.  
            if ~pars.alignMaxProject
                fidFrame = find(isFidChn,1,'first');
                fiducialFrames{h}  = ReadDax([currFolder,daxFiles{pars.fov}],...
                    'verbose',false,'startFrame',fidFrame,'endFrame',fidFrame);
            else
                dax = ReadDax([currFolder,daxFiles{pars.fov}],'verbose',false);
                fiducialFrames{h} = max(dax(:,:,isFidChn),[],3);
                % save max project
                
            end
        catch er   
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
    [fiducialAlignFrames,regData] = RegisterImagesFast(fiducialFrames,'parameters',pars);
end

% If registered data already exists, just load that
%   we don't save the first in focus frame, just the z-stack, so now we
%   need to grab a frame half-way through the stack, 'midFrame', which
%   should be in focus.  
if skip
    fiducialAlignFrames = cell(numHybes,1); 
    for h=1:numHybes
        fiducialAlignFrames{h}  = ReadDax([saveFolder,fidNames{h}],...
                'verbose',false,'startFrame',midFrame,'endFrame',midFrame);
    end
end

dataOut.regData = regData;
dataOut.dataFolder = dataFolders;
dataOut.saveFolder = saveFolder;
dataOut.eTable = eTable;
dataOut.fiducialFrames = fiducialAlignFrames;  
dataOut.numDataChns = numDataChns;
dataOut.frameChannels = frameChannels; % currently assuming constant 

