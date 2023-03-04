function [daxInHybFov,daxInfoHybFov,rawNameInHybFov,maxNameInHybFov] = LoadDaxFromEtable(eTableXLS,varargin)
% load a dax file using the eTable to determine 
%  hybNumber, fov, dataType, and dataChannel
% by default, the max-projection is returned.  Set 'maxProject',false, to
% return the full 3D data from the requested channel.  Max projections are
% saved in a protected name structure when computed, so they will be loaded
% faster in the future. They are saved in the same folder with the data.
% If multiple hybs, multiple FOVs, or multiple data channels are requested,
% the fxn returns a data structure. 
% 
% 
% daxInHybFov = nHybs x nFOVs x datachannels
% dataType 'all'  daxInHybFov  datachannels are fid, dataRed, dataBlue
% 
% To add:
%  CheckRegData -- this should be checked upstream, not here.
%
%

defaults = cell(0,3);
defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'all'}; % changed from 'data'
defaults(end+1,:) = {'dataFolder','string',''};
defaults(end+1,:) = {'saveFolder','string',''}; % place to save the registered data files 
defaults(end+1,:) = {'hybNumber', 'integer', 1};  
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'dataChannel', 'integer', inf}; % 1
defaults(end+1,:) = {'hybFolder', 'string', ''}; 
defaults(end+1,:) = {'maxProject', 'boolean', true}; % return max projection
defaults(end+1,:) = {'overwrite', 'boolean', false}; % return max projection
defaults(end+1,:) = {'selectFrames', 'integer', 0}; % 0 for all frames
defaults(end+1,:) = {'daxRootDefault', 'string', 'ConvZscan*.dax'}; 
defaults(end+1,:) = {'saveProject', 'boolean', true};  % write a maxProjection if one does not exist
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryverbose', 'boolean', false}; 
defaults(end+1,:) = {'checkRegData', 'boolean', true};  % see if regData exists, if so, just return names 
defaults(end+1,:) = {'readDax', 'boolean', true}; % make false to just return names 
defaults(end+1,:) = {'fixDrift', 'boolean', false}; % optionally, try to find and fovNNN_regData.csv and correct drift
defaults(end+1,:) = {'driftFolder', 'string', ''}; % folder in which to find and fovNNN_regData.csv 
defaults(end+1,:) = {'hybType',{'all','B','H','R','T'},'all'}; % 
defaults(end+1,:) = {'simplifyOutput', 'boolean', true}; % convert cell-array to image matrix if only 1 movie is requested
pars = ParseVariableArguments(varargin,defaults,mfilename);

%%
% eTableXLS = 'T:\2019-04-02_DemoCT3_dat\ExperimentLayout.xlsx'
% Required aspects of data layout:
%   eTableXLS sits in dataFolder
%   dataFolder contains hybFolders, whose names are listed in the hybFolder
%   column of the eTable.

% allow either a pre-loaded table or filepath to table to be loaded
%   this saves time if this function is called frequently;
if isempty(eTableXLS)
    [eFile,folder] = uigetfile('*.xlsx');
   eTableXLS = [folder,eFile]; 
end
if ischar(eTableXLS)
    eTable = readtable(eTableXLS);
    dataFolder = [fileparts(eTableXLS),filesep];
elseif istable(eTableXLS)
    eTable = eTableXLS;
    dataFolder = pars.dataFolder;
end

[~,frameChannels,channels,~,dataChns] = GetFidChnFromTable(eTable);  %
fidChnName = setdiff(channels,dataChns);
hybFolders = eTable.FolderName;

if isempty(pars.saveFolder)
    saveFolder = dataFolder;
else
    saveFolder = pars.saveFolder;
end

% Get requested hyb set
if ~isempty(pars.hybFolder)
    hybs = find(strcmp(pars.hybFolder,hybFolders));
    if isempty(hybs)
       error(['could not find folder ',pars.hybFolder,' in ',dataFolder]); 
    end
else
    hybs = pars.hybNumber;
end
if isinf(hybs)
    hybs = 1:length(hybFolders);
end
if strcmp(pars.hybType,'B') % return only barcode hybs
    barcodeHybes = strcmp(eTable.DataType,'B');
    hybs = Row(intersect(hybs,find(barcodeHybes)));
end

% get requested FOV
currFolder = [dataFolder,hybFolders{1},filesep];
rawDaxFiles = cellstr(ls([currFolder,pars.daxRootDefault]));
if isinf(pars.fov)
    pars.fov = 1:length(rawDaxFiles);
end


% Check if regData.csv files already exist
%   if so, offer to skip
if pars.checkRegData && ~isempty(pars.saveFolder) && pars.readDax
    regDatas = ls([saveFolder,'fov*regData.csv']);
    if ~isempty(regDatas)
        fovsWithReg = str2num(regDatas(:,4:6))'; %#ok<ST2NM>
        allDone = isempty(setdiff(fovsWithReg,pars.fov));
        if allDone 
            answer = questdlg(['Found regData.csv for all ',num2str(length(pars.fov)),...
            ' requested FOVs. Skip loading max-projection?'], ... % question
                        'Prompt ', ...  % pop-up label
                        'Yes','No','No'); % op1 op2 default
            if strcmp(answer,'Yes')
                pars.readDax = false;
            end
        end
    end
end

% determine requested frames and assign requested flags
if strcmp(pars.dataType,'all')
    numDatas     = length(channels);
    reqFrameList = cell(numDatas,1);
    dataFlagList = cell(numDatas,1);
    reqFrameList{1} = strcmp(fidChnName,frameChannels);
    dataFlagList{1} = 'fidMax';
    for d=1:numDatas-1
        reqFrameList{d+1} = strcmp(dataChns{d},frameChannels);
        dataFlagList{d+1} = ['dat',dataChns{d},'Max'];
    end
elseif strcmp(pars.dataType,'data')
    numDatas     = length(dataChns);
    reqFrameList = cell(numDatas,1);
    dataFlagList = cell(numDatas,1);
    for d=1:numDatas
        reqFrameList{d} = strcmp(dataChns{d},frameChannels);
        dataFlagList{d} = ['dat',dataChns{d},'Max'];
    end
elseif strcmp(pars.dataType,'fiducial')
    numDatas        = 1;
    reqFrameList{1} = strcmp(fidChnName,frameChannels);
    dataFlagList{1} = 'fidMax';
else
    cprintf([1 0 0],'valid data types: "all", "fiducial", "data"');
    error(['found unrecognized dataType ',pars.dataType]);
end

daxInHybFov     = cell(max(hybs),max(pars.fov),numDatas);
daxInfoHybFov   = cell(max(hybs),max(pars.fov),numDatas);
rawNameInHybFov = cell(max(hybs),max(pars.fov),numDatas);
maxNameInHybFov = cell(max(hybs),max(pars.fov),numDatas);
% Loop over FOVs (default just FOV 1)
for f = pars.fov
    if length(pars.fov) > 5 && pars.verbose
        disp(['Loading data ',num2str(f/max(pars.fov)*100),'% complete']); 
    end
    % Loop over hybs (default, just hyb 1)
    for h = hybs
        % Loop over data types  (default just data)
        for d=1:numDatas
            reqFrames   = reqFrameList{d};
            maxFlag     = [dataFlagList{d},'_'];           
            currFolder  = [dataFolder,hybFolders{h},filesep];
            rawDaxFiles = FindFiles([currFolder,pars.daxRootDefault],'fullPath',false); % cellstr(ls([currFolder,pars.daxRootDefault]));
            currDax     = [currFolder,rawDaxFiles{f}];
            maxName     = [maxFlag,rawDaxFiles{f}];
            maxFile     = [currFolder,maxName];
            try
                if pars.readDax
                    if pars.maxProject
                        if exist(maxFile,'file') && ~pars.overwrite
                            [dax,infoOut] = ReadDax(maxFile,'verbose',pars.veryverbose);
                        else
                            % This is inefficient for the first pass --
                            %   The same dax file gets read in 3 seperate
                            %   times (though matlab speeds rpt reading up
                            %   a bit) in order to separately save the
                            %   fiducial, dataChn1, dataChn2.  
                            % proposed soln:
                            % Maybe, given that we have the eTable, if we
                            % are asked for only a max project and we don't
                            % have one, we create all 3, even if we only
                            % return 1. 
                            [dax,infoFile] = ReadDax(currDax,'verbose',pars.veryverbose);
                            infoOut = infoFile;
                            dax = max(dax(:,:,reqFrames),[],3);
                            if pars.saveProject % save max project                                  
                                infoOut.number_of_frames= 1;
                                infoOut.localName = regexprep(maxName,'.dax','.inf');
                                WriteDAXFiles(dax,infoOut,'verbose',pars.veryverbose);
                            end
                        end
                    elseif pars.selectFrames == 0
                        [dax,infoOut] = ReadDax(currDax,'verbose',false);              
                        dax = dax(:,:,reqFrames);
                    elseif pars.selectFrames ~= 0
                        dax = ReadDax(currDax,'verbose',false,'startFrame',min(pars.selectFrames),'endFrame',max(pars.selectFrames));              
                    end                
                    daxInHybFov{h,f,d} = dax;
                    daxInfoHybFov{h,f,d} = infoOut;
                end
                rawNameInHybFov{h,f,d} = currDax;
                maxNameInHybFov{h,f,d} = maxFile;
            catch er   
                if pars.verbose
                    warning(er.message);
                    disp(['failed to load ', currDax]); 
                end
            end
        end
    end
end

% --- CorrectDrift
if pars.fixDrift && pars.readDax
    if isempty(pars.driftFolder)
        pars.driftFolder = pars.saveFolder;
    end
    for f=pars.fov
        try
            regTable = readtable([pars.driftFolder,'fov',num2str(f,'%03d'),'_regData.csv']);
            regStruct = table2struct(regTable);
            for h= hybs
                for d=1:numDatas
                    daxInHybFov{h,f,d} = ApplyReg(daxInHybFov{h,f,d},regStruct(h));
                end
            end
        catch er
            warning(['Requested drift fix for FOV ', num2str(f), ' unable to fix drift']); 
            warning(er.message);
        end
    end
end

%--- Correct Chromatic Aberration 


% if only a single movie is requested, return the dax file. 
if length(pars.fov) == 1 && length(hybs) == 1 &&  numDatas == 1 && pars.simplifyOutput
   daxInHybFov = cat(3,daxInHybFov{hybs,pars.fov,:});
end
if pars.verbose
   disp('data loaded'); 
end
