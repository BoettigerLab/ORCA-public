function [fiducialFrames,eTable,foundFiles] = ChrTracer3_LoadExpTable(expTableXLS,saveFolder,varargin)
% Load experiment table and max-project data. Return max-project (non-drift
% corrected) fiduical data for registration
% 
% To fix: 
% should probably skip all the way to spot-fitting if the aligned data is
% detected.

defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryverbose', 'boolean', false}; 
defaults(end+1,:) = {'dataFolder', 'string', ''}; 
defaults(end+1,:) = {'dataTypes', 'cell', {'any'}}; 
defaults(end+1,:) = {'daxRootDefault', 'string', 'ConvZscan*.dax'}; 
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
defaults(end+1,:) = {'maxProject', 'boolean', true}; 
defaults(end+1,:) = {'saveProject', 'boolean', true}; 
defaults(end+1,:) = {'selectFOVs', 'freeType', inf}; 
defaults(end+1,:) = {'stopOnError', 'boolean', true}; 
%  pars = ParseVariableArguments([],defaults,'ChrTracer3_LoadExpTable');
pars = ParseVariableArguments(varargin,defaults,mfilename);

%%
% Load the experiment table
eTable = readtable(expTableXLS);

if isempty(pars.dataFolder)
    dataFolder = fileparts(expTableXLS); % get folder name;
else
    dataFolder = pars.dataFolder;
end

%% Check OPTIONAL Inputs
% ---- OP1: load only a subset of data  ----- 
% load experiment table to know what probes are in which folders and
% whether they contain hybe data, repeat-hybe data, or toehold data
hybFolderID = false(height(eTable),1);
if ~any(strcmp(pars.dataTypes,'any'))  % if specific datatypes are requested
    for i=1:length(pars.dataTypes)
        hybFolderID = hybFolderID | strcmp(eTable.DataType,pars.dataTypes{i});
    end
else
    hybFolderID = true(height(eTable),1);
end
numHybes = sum(hybFolderID);
eTable = eTable(hybFolderID,:);
hybFolders = eTable.FolderName;

% ----- OP2: specify full datapath in Excel file ---
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

% ---- OP3: specify daxfile names ----
daxNameCol = contains(eTable.Properties.VariableNames,'DaxRoot');
if any(daxNameCol)
    daxRoots = eTable{:,daxNameCol};
else
    daxRoots = repmat({pars.daxRootDefault},numHybes,1);
end

% find total number of FOVs by looking in hyb folder 1
if ~any(pars.selectFOVs) || any(isinf(pars.selectFOVs))
    h = 1;
    currFolder = [dataFolders{h},'\',hybFolders{h},'\'];
    daxFiles =  cellstr(ls([currFolder,daxRoots{h}]));
    numFOVs = length(daxFiles);
    pars.selectFOVs  = 1:numFOVs;
else
    numFOVs = length(pars.selectFOVs);
end

%% Check if data is already written
% convert annotations of all frames
[isFidChn,~,channels] = GetFidChnFromTable(eTable);
numChns = length(channels);
numDataChns = numChns -1; 

[foundFiles,fovToSkip] = ChrTracer3_FindAlignedFiles(saveFolder,pars.selectFOVs,numHybes,numDataChns);
disp('Experiment table loaded');

% % update numFOVs
selectFOVs = pars.selectFOVs;
% selectFOVs(fovToSkip) = [];
% numFOVs = length(selectFOVs); 

skip = numFOVs < 1;

%% load data (default is load xy max project)
fiducialFrames = cell(numHybes,numFOVs); 
if ~skip  
    disp('extracting fiducial data for global drift correction');
    disp(['processing ',num2str(numHybes),' hybes containing ',num2str(numFOVs),' total FOVs']);
    if numFOVs > 5 && numHybes > 5
        disp('this may take some time...');
    end
    for h=1:numHybes 
        currFolder = [dataFolders{h},'\',hybFolders{h},'\'];
        daxFiles =  cellstr(ls([currFolder,daxRoots{h}]));
        for f=selectFOVs
            currDax = [currFolder,daxFiles{f}];
            try
                if pars.maxProject
                    maxName = ['fidMax_',daxFiles{f}];
                    maxFile = [currFolder,maxName];
                    if exist(maxFile,'file')
                        dax = ReadDax(maxFile,'verbose',pars.veryverbose);
                    else
                        [dax,infoFile] = ReadDax(currDax,'verbose',pars.veryverbose);
                        dax = max(dax(:,:,isFidChn),[],3);
                        if pars.saveProject % save max project  
                            infoOut = infoFile;
                            infoOut.number_of_frames= 1;
                            infoOut.localName = regexprep(maxName,'.dax','.inf');
                            WriteDAXFiles(dax,infoOut,'verbose',pars.veryverbose);
                        end
                    end  
                else
                    fidFrame = find(isFidChn,1,'first');
                    dax = ReadDax(currDax,'verbose',false,'startFrame',fidFrame,'endFrame',fidFrame);              
                end
                fiducialFrames{h,f} = dax; 
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
        if pars.verbose
           disp([num2str(h/numHybes*100,3),'% complete']); 
        end
    end
end

