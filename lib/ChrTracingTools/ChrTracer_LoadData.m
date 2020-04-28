function [dataOut,pars] = ChrTracer_LoadData(expTableXLS,varargin)
% Load an expTableXLS table and use it to organize hybes into a fiducial
% stack and data stack.  Correct both stacks using image correlation of
% fiducial frames. 


% dataFolder = 'E:\Nasa\2017-03-23_K562-chr21_combined';

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
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
% load data pars
defaults(end+1,:) = {'filename', 'string', 'ConvZscan'}; 
defaults(end+1,:) = {'dataTypes','cell',{'Hyb','Repeat','Toehold','H','T','R','any'}};
% registration parameters
defaults(end+1,:) = {'rotation','boolean',false}; % also compute and correct rotation angle
defaults(end+1,:) = {'maxD','positive',15}; % distance in pixels to match objects prior to computing rotation
defaults(end+1,:) = {'alignmentBoxWidth', 'positive', inf}; % pixels.  size of region to use for frame alignment
defaults(end+1,:) = {'alignUpsample', 'positive', 1}; % upsample data for coarse alignment (slow, should be unnecessary);
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .999}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignToFirst', 'boolean', true}; % align to first image? if false will align to previous non-empty image data
defaults(end+1,:) = {'minFracObj','fraction',.65}; % min fraction of objects from previous hybe found in target hybe to be acceptable 
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'previousAlignFrames','array',[]}; % previous align frames
defaults(end+1,:) = {'targetHybes','integer',[]};
defaults(end+1,:) = {'goodHybes','boolean',[]};
defaults(end+1,:) = {'flattenBackground','nonnegative',0};
defaults(end+1,:) = {'corrAngles','freeType',[]}; % rotate angles to check for better alignment

pars = ParseVariableArguments(varargin, defaults, mfilename);
% pars = ParseVariableArguments([], defaults, mfilename);

%  expTableXLS = 'E:\Nasa\2017-05-23_IMR90-chr21\ExperimentLayout.xlsx'
% expTableXLS = [dataFolder,'\ExperimentLayout.xlsx'];
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
saveFolder = SetFigureSavePath( [dataFolder,'_Analyzed\'],'makeDir',true);

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



% Load Data
fiducialNames = cell(numHybes,1);
dataNames = cell(numHybes,1);
fiducialFrames = cell(numHybes,1);
dataFrames = cell(numHybes,numDataChns);
%fiducialFrames = zeros(yDim,xDim,numFrames,numHybes,'uint16'); % it should be faster to have these as matrices rather than cells
%dataFrames = zeros(yDim,xDim,numFrames,numHybes*numDataChns,'uint16'); % it should be faster to have these as matrices rather than cells
goodHybes = true(1,numHybes);
for h=1:numHybes % load data
    if pars.verbose
        disp(['loading data from hybe ',num2str(h), '  folder: ',hybFolders{h}]);
    end
    currFolder = [dataFolder{h},'\',hybFolders{h},'\'];
    daxFiles =  cellstr(ls([currFolder,pars.filename,'*.dax']));
    try
        fiducialNames{h} = [currFolder,regexprep(daxFiles{pars.fov},'.dax',''),'_fid'];
        dataNames{h} = [currFolder,regexprep(daxFiles{pars.fov},'.dax',''),'_data'];
        dax = ReadDax([currFolder,daxFiles{pars.fov}],'verbose',false);
        keepFrames = str2num(regexprep(eTable.framesToKeep{h},'[\[,\]]','')); %#ok<ST2NM>
        fidFrames = (1+fidChn):numChns:max(keepFrames);
        fiducialFrames{h} = dax(:,:,intersect(fidFrames,keepFrames));
        for n=1:numDataChns
            dataKeepFrames = (1+dataChns(n)):numChns:max(keepFrames);
            dataFrames{h,n} = dax(:,:,intersect(dataKeepFrames,keepFrames)); % fix me
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
[fiducialAlignFrames,dataAlignFrames,goodHybes] = RegisterImages(fiducialFrames,dataFrames,'parameters',pars,'goodHybes',goodHybes);

% data returned by function
dataOut.rawFiducialFrames = fiducialFrames;
dataOut.rawDataFrames = dataAlignFrames;
dataOut.fiducialFrames = fiducialAlignFrames; % 
dataOut.dataFrames = dataAlignFrames; % 
dataOut.goodHybes = goodHybes; % 
dataOut.dataFolder = dataFolder;
dataOut.saveFolder = saveFolder;
dataOut.eTable = eTable;

