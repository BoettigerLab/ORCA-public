function [fileNames,fidFlatFrames,datFlatFrames] = ChrTracer3_SaveRegData(eTable,regData,varargin)
% this replaces ChrTracer3_SaveAlignedData
% this only saves the registered DaxFiles
% MemoryMapping will be handled in a separate step. 
%   Here the data is processed 1 hybe at a time.  Memory maps will allow
%   all the data from all hybes to read in at once without overloading the
%   system RAM.
% 

defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryverbose', 'boolean', false}; 
defaults(end+1,:) = {'filename', 'string', 'ConvZscan'};  % eTable value should be able to override this.  See LoadExpTable.  Needs fix 
defaults(end+1,:) = {'useTableName','boolean',true}; % if "names" are given in the table file, use these as file names.
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
defaults(end+1,:) = {'selectFOVs','freeType',inf};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'dataFolder','string',''};
defaults(end+1,:) = {'refHybe','integer',1};
defaults(end+1,:) = {'loadFlatData','boolean',true};
defaults(end+1,:) = {'loadFlatFid','boolean',true};
defaults(end+1,:) = {'saveMaxProject','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% To Check:
% maybe some of these next steps should be passed on from LoadExpTable

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
[isFidChannel,frameChannels,~,~,dataChns] = GetFidChnFromTable(eTable);
nBufFrames = eTable.bufferFrames(1);
totFrames  = eTable.totalFrames(1);
keepFrames = nBufFrames+1:totFrames-nBufFrames;
numDataChns = length(dataChns);

%------------------- File names ---------------------------------------
% if hyb names are supplied (e.g. gene names, use them)
% otherwise, just record hyb numbers
try
    nameField = find(contains(eTable.Properties.VariableNames,'name','IgnoreCase',true));
    if length(nameField)>1
        nameField = nameField(end);
    end
    tableNames = regexprep(eTable{:,nameField},{'_',' '},{'',''});
catch
    if pars.verbose
       disp('error reading hyb names from table. Just recording numbers'); 
       pars.useTableName = false;
       tableNames = {};
    end
end
if pars.useTableName
    hybeNames = tableNames;
else
    hybeNames = cellstr(strcat('h',num2str( (1:numHybes)','%03d')));
end



if ~any(pars.selectFOVs) || any(isinf(pars.selectFOVs))
    numFOVs = length(regData);
    pars.selectFOVs = 1:numFOVs;
else
    numFOVs = length(pars.selectFOVs);
end

% initialize output datastructures
fileNames(numFOVs).fidNames = cell(numHybes,1);
fileNames(numFOVs).datNames = cell(numHybes,numDataChns);
fidFlatFrames = cell(numHybes,1,numFOVs);
datFlatFrames = cell(numHybes,numDataChns,numFOVs); 

% select Hybes 
% (ordered with ref-hyb first, entered still in to referenced position)
allHybes = 1:numHybes;
allHybes(pars.refHybe) = []; 
allHybes = [pars.refHybe,allHybes]; 

% check if files are already written
selectFOVs = pars.selectFOVs;
[foundFiles,fovToSkip] = ChrTracer3_FindAlignedFiles(pars.saveFolder,selectFOVs,numHybes,numDataChns);
if ~pars.overwrite
    selectFOVs(fovToSkip) = []; % remove previously analyzed fovs from list of to drift correct   
end
fileNames(fovToSkip) = foundFiles(fovToSkip); % add file names from previously analyzed data to list of file names
idToSkip = find(fovToSkip);

noData = cellfun(@isempty, regData);
noData(fovToSkip) = false;
if any(noData)
    skipMe = Column(find(noData))';
    cprintf([1 .25 0],['Warning: did not find registration data for FOV ',num2str(skipMe),' . skipping this FOV']);
    idToSkip = [idToSkip,skipMe];
end


            
if pars.loadFlatFid
    if length(idToSkip)>1 && pars.verbose
        disp('loading previously aligned fiducial projections...');
    end
    for f = idToSkip
        for h=pars.refHybe
            try
                [~,maxName]=fileparts(fileNames(f).fidNames{h});
                maxName = ['max_',maxName,'.inf']; %#ok<AGROW>
                hasMax = exist([pars.saveFolder,filesep,maxName],'file')==2;
                if hasMax 
                    fidFlatFrames{h,1,f} = ReadDax([pars.saveFolder,maxName],'verbose',pars.veryverbose);
                else
                    [dax,infoFile] = ReadDax(foundFiles(f).fidNames{h,1},'verbose',pars.veryverbose);
                    fidFlatFrames{h,1,f} = max(dax,[],3);
                    if pars.saveMaxProject
                        infoOut = infoFile;
                        infoOut.number_of_frames= 1;
                        infoOut.localName = maxName;
                        WriteDAXFiles(fidFlatFrames{h,1,f},infoOut,'verbose',pars.verbose);
                    end
                end
            catch er
                warning(er.getReport); 
                disp(['something went wrong on FOV ',num2str(f)]);
                idToSkip(idToSkip==f) = [];
            end
        end
    end
end
if pars.loadFlatData 
    if length(idToSkip)>1 && pars.verbose
        disp('loading previously aligned data projections...');
    end
    for f = idToSkip
        for h=allHybes
            for n=1:numDataChns
                try
                    [~,maxName]=fileparts(fileNames(f).datNames{h,n});
                    maxName = ['max_',maxName,'.inf']; %#ok<AGROW>
                    hasMax = exist([pars.saveFolder,filesep,maxName],'file')==2;
                    if hasMax 
                        datFlatFrames{h,n,f} = ReadDax([pars.saveFolder,maxName],'verbose',pars.veryverbose);
                    else
                        [dax,infoFile] = ReadDax(foundFiles(f).datNames{h,n},'verbose',pars.veryverbose);
                        datFlatFrames{h,n,f} = max(dax,[],3);
                        if pars.saveMaxProject
                            infoOut = infoFile;
                            infoOut.number_of_frames = 1;
                            infoOut.localName = maxName;
                            WriteDAXFiles(datFlatFrames{h,n,f},infoOut,'verbose',pars.verbose);
                        end
                    end
                catch er
                    warning(er.getReport); 
                    disp(['something went wrong on FOV ',num2str(f)]);
                    idToSkip(idToSkip==f) = [];
                end
            end
        end
    end
end
    
for f=selectFOVs
    if isempty(regData)
        warning('No registration data found. Run "Fix Global Drift" and try again');
        break;
    end
    
    try
        tic   
        %------------------ Save new aligned dax file ---------------------------%
        % loop through all hybes, but do refHybe first
        for h=allHybes
            if pars.verbose
                disp(['loading data from fov ',num2str(f),' hybe ',num2str(h), '  folder: ',hybFolders{h}]);
            end
            % check if file already exists
            fidInfName = ['fov',num2str(f,'%03d'),'_',hybeNames{h},'_fid','.inf'];
            if ( exist([pars.saveFolder,fidInfName],'file')==0 || pars.overwrite ) || h==pars.refHybe

                % Read in data for this hybe 
                currFolder = [dataFolders{h},filesep,hybFolders{h},filesep];
                daxFiles =  cellstr(ls([currFolder,pars.filename,'*.dax']));
                [dax,info] = ReadDax([currFolder,daxFiles{f}],'verbose',false,...
                    'startFrame',nBufFrames+1,'endFrame',totFrames-nBufFrames);

                %-----Update the info file for these dax movies------------
                info.localPath = pars.saveFolder;
                % Stage_X/Stage_Y varies between hybes for unknown reasons
                %   since we have corrected positions to subpixel accuracy we want
                %   to record the corrected positions, which are those of
                %   the first "reference" hyb. The originals will be added
                %   to the notes.
                if h==pars.refHybe
                   refStageX = info.Stage_X;
                   refStageY = info.Stage_Y;
                else
                    info.notes = [info.notes,' Orig Stage_XY=(',num2str(info.Stage_X),',',num2str(info.Stage_Y),')'];
                    info.Stage_X = refStageX;
                    info.Stage_Y = refStageY;       
                end

                % split out data frames and fiducial frames 
                fiducialFrames = dax(:,:,isFidChannel(keepFrames));
                dataFrames = cell(numDataChns,1);
                for n=1:numDataChns 
                    isCurrDat = StringFind(frameChannels,dataChns{n},'boolean',true);
                    dataFrames{n} = dax(:,:,isCurrDat(keepFrames));
                end

                % --- Align and save with correct infoFile and updated stage positions
                fidReg = ApplyReg(fiducialFrames,regData{f}(h));
                fidInfo = info;
                fidInfo.number_of_frames = size(fidReg,3); 
                fidInfo.localName = fidInfName;
                WriteDAXFiles(fidReg,fidInfo,'verbose',pars.veryverbose);
                fileNames(f).fidNames{h} = [pars.saveFolder,fidInfo.localName(1:end-4),'.dax'];
                if pars.loadFlatFid
                    maxName = ['max_',fidInfo.localName];
                    fidFlatFrames{h,f} = max(fidReg,[],3); % record max project
                    if pars.saveMaxProject && (~exist([pars.saveFolder,filesep,maxName],'file') || pars.overwrite)
                        fidInfo.number_of_frames= 1;
                        fidInfo.localName = maxName;
                        WriteDAXFiles(fidFlatFrames{h,f},fidInfo,'verbose',pars.veryverbose);
                    end
                end

                % Fix data-channel
                for n = 1:numDataChns
                    % fix XY-drift (maybe we should do something about z-drift
                    % in these steps too?)
                    datReg = ApplyReg(dataFrames{n},regData{f}(h));

                    % setup info file and file name
                    datInfo = info;
                    datInfo.number_of_frames = size(fidReg,3);
                    datInfo.localName = ['fov',num2str(f,'%03d'),'_',hybeNames{h},'_data',num2str(n),'.inf']; 
                    datInfo.notes = dataChns{n};

                    % write the dax file and record the filename/location.
                    WriteDAXFiles(datReg,datInfo,'verbose',pars.veryverbose);
                    fileNames(f).datNames{h,n} = [pars.saveFolder,datInfo.localName(1:end-4),'.dax'];

                    % save max project file
                    if pars.loadFlatData
                        maxName = ['max_',datInfo.localName];
                        datFlatFrames{h,n,f} = max(datReg,[],3);
                        if pars.saveMaxProject && (~exist([pars.saveFolder,filesep,maxName],'file') || pars.overwrite) % write max projection for speed.
                            datInfo.number_of_frames= 1;
                            datInfo.localName = maxName;
                            WriteDAXFiles(datFlatFrames{h,n,f},datInfo,'verbose',pars.veryverbose);
                        end
                    end
                end  
            else
                fileNames(f).fidNames{h} = [pars.saveFolder,fidInfName(1:end-4),'.dax'];
                for n = 1:numDataChns
                    localName = ['fov',num2str(f,'%03d'),'_',hybeNames{h},'_data',num2str(n),'.inf']; 
                    fileNames(f).datNames{h,n} = [pars.saveFolder,localName(1:end-4),'.dax'];
                end  
                if pars.verbose
                    disp(['found existing ', fidInfName ,'.  skipping to next image']);
                end
            end
        end
    toc
    catch er
        warning(er.getReport);
        disp(['something went wrong on FOV ',num2str(f)]);
    end

end
%------------------------------------------------------------------------% 

