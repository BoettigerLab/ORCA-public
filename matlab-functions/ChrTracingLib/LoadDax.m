function [chnStk,imProps,imData] = LoadDax(daxFileName,varargin)
% this function uses the xml data to guess the channel organization
% following curent lab conventions:
%   Z-stack
%   optional buffer frames
%   alternate channels at each setp of the z-stack
%   movie saved as a single z-stack
% 
%  it also uses the xml file mosaic data to figure out the xy-pixel size
%  and the xml file z-stack data to figure out the z-step size
%   (Though in some cases this was miscalibrated in Hal and will need
%   updating). 
% 
% -------------------------------------------------------------------------
% 
%% Assumption for autopopulation
% 1. number of channels and channel names are extracted from the shutters
%      requires taht shutters file reads shutters_chn1name_chn2name_series
% 2. hyb number is extracted from folder name 

defaults = cell(0,3);
% key parameters
% defaults(end+1,:) = {'imProps', 'struct',[]};
defaults(end+1,:) = {'verbose', 'boolean',true};
defaults(end+1,:) = {'overwrite', 'boolean',false};
defaults(end+1,:) = {'maxProject', 'boolean',false};
defaults(end+1,:) = {'channel', 'integer',0}; % 0 = load all.  otherwise [1,3] loads first and third
defaults(end+1,:) = {'skipFirst', 'integer',0};
defaults(end+1,:) = {'writeData', 'boolean',true};  
defaults(end+1,:) = {'justImProps', 'boolean',false};  
defaults(end+1,:) = {'driftFolder', 'string',''};  % if empty no drift correct, otherwise, search for an fov_alignTable.csv
defaults(end+1,:) = {'hyb', 'integer',0}; % extract from folder name by default
pars = ParseVariableArguments(varargin,defaults,mfilename); 

daxFileName = regexprep(daxFileName,{'.xml','.inf'},{'.dax','.dax'}); % allow xml of inf input
xmlFile = regexprep(daxFileName,'.dax','.xml');
xmlFile = regexprep(xmlFile,'_C1',''); % remove camera flags;
xmlFile = regexprep(xmlFile,'_C2',''); % remove camera flags;
if exist(xmlFile,'file')
    xmlStruct = ReadXML(xmlFile);
else
    warning('no xml file detected, defaulting to ReadDax');
    chnStk = ReadDax(daxFileName);
    imData = [];
    imProps = [];
    return
end

[imProps.dataFolder,imProps.daxName] = fileparts(daxFileName);
imProps.nFrames  =       xmlStruct.settings.acquisition.number_frames;
imProps.imSize =         [xmlStruct.settings.camera1.x_pixels,xmlStruct.settings.camera1.y_pixels];
imProps.camera.flipHorizontal = xmlStruct.settings.camera1.flip_horizontal;
imProps.camera.flipVertical = xmlStruct.settings.camera1.flip_vertical;
imProps.camera.transpose = xmlStruct.settings.camera1.transpose;
imProps.mosaic.flipHorizontal = xmlStruct.settings.mosaic.flip_horizontal;
imProps.mosaic.flipVertical = xmlStruct.settings.mosaic.flip_vertical;
imProps.mosaic.transpose = xmlStruct.settings.mosaic.transpose;

try
    imProps.stageXY  =       str2num(xmlStruct.settings.acquisition.stage_position); %#ok<ST2NM>
catch
    warning('stage position missing')
    imProps.stageXY  = [NaN NaN];
end

% imProps.nChns =[];
% imProps.preStartFrames = [];
% imProps.zStepSize = [];
% imProps.zSteps  = [];
% imProps.postImageFrames = [];

% software_z_scan
if isfield(xmlStruct.settings.focuslock,'software_z_scan')
    % these parameters have changed names in different versions of hal
    if isfield(xmlStruct.settings.focuslock.software_z_scan,'frames_to_pause')
        imProps.nChns = xmlStruct.settings.focuslock.software_z_scan.frames_to_pause; % UPDATE me for Old Hal
    elseif isfield(xmlStruct.settings.focuslock.software_z_scan,'frames_per_step')
        imProps.nChns =  xmlStruct.settings.focuslock.software_z_scan.frames_per_step; % UPDATE me for Old Hal
    else
        warning('unable to determine nChns');
        imProps.nChns = NaN;
    end
    if isfield(xmlStruct.settings.focuslock.software_z_scan,'deadtime')  % Update me (sometimes this gets different names)
        imProps.preStartFrames = xmlStruct.settings.focuslock.software_z_scan.deadtime;
    else
        warning('unable to determine deadtime, assuming 10')
        imProps.preStartFrames = 10;
    end
    % these software_z_scan properties have not changed names
    imProps.zStepSize =      xmlStruct.settings.focuslock.software_z_scan.step_size;  % step size in nm, as recorded by hal
    imProps.zSteps =         2*xmlStruct.settings.focuslock.software_z_scan.range/imProps.zStepSize;  % total scan is +/- scan range  
else
    % python2 hal doesn't have 'software_z_scan'
    imProps.nChns =          xmlStruct.settings.focuslock.zscan_frames_to_pause; % UPDATE me for Old Hal  
    warning('unable to determine deadtime, assuming 10')
    imProps.preStartFrames = 10; % UPDATE me for Old Hal
    imProps.zStepSize =      xmlStruct.settings.focuslock.zscan_step*1e3;  % step size in nm, as recorded by hal
    imProps.zSteps =         2*xmlStruct.settings.focuslock.zscan_stop*1e3/imProps.zStepSize;  % total scan is +/- scan range  
end

if isstr(imProps.nFrames)
    imProps.nFrames = str2double(imProps.nFrames);
end

    
imProps.postImageFrames= imProps.nFrames -imProps.preStartFrames - imProps.zSteps*imProps.nChns;
% auto-detect pixel size
objFields = fields(xmlStruct.settings.mosaic);
objUsed = strcmp(objFields,xmlStruct.settings.mosaic.objective);
magData = strsplit(xmlStruct.settings.mosaic.(objFields{objUsed}),',');
imProps.xy2um = str2double(magData(2));  % xy-pixel size in um
imProps.z2um = imProps.zStepSize/1e3;  % z-pixel size in um 

% auto-detect FOV
nameParts = strsplit(imProps.daxName,'_');
imProps.fov = str2double(nameParts{end});
% auto-detect Hyb (might give NaN if 
nameParts = strsplit(imProps.dataFolder,filesep);
if isempty(nameParts{end})
    nameParts = nameParts{end-1};
else
    nameParts = nameParts{end};
end
nameParts = strsplit(nameParts,'_');
imProps.Hyb = str2double(nameParts{end}); 


% autodetect channel names from shutters
shutterName =           xmlStruct.settings.illumination.shutters;
[~,shutter] = fileparts(shutterName);
shutter = regexprep(shutter,{'shutters_','_series','_parallel'},{'','',''});
shutter = regexprep(shutter,'confocal_',''); % confocal flag is not a channel
imProps.chnNames = strsplit(shutter,'_');
imProps.nChns = length(imProps.chnNames);
imProps.projection = false;

if pars.justImProps
    chnStk = imProps;
else
    %% Load image data for FOV 1
    % write max projections for each channel now. Future functions that require
    % only max projections will check if these exist and load them in place of
    % loading the full data and re-projecting it, reducing computational time

    if pars.channel == 0
        chnsToLoad = 1:imProps.nChns; % load all channels
    elseif isinf(pars.channel)
        chnsToLoad = imProps.nChns; % load last channel (typically fiduical)
    else
        chnsToLoad = Row(pars.channel); % load select channels
    end

    imProps.chnsLoaded = length(chnsToLoad); 

    % a folder to save Max Projections
    maxProjectionFolder = [imProps.dataFolder,'\maxProjections\'];
    if ~exist(maxProjectionFolder,'dir')
        mkdir(maxProjectionFolder);   
    end

    % An enforced naming convention for Max Projections
    nameRoot = strsplit(imProps.daxName,'_'); % hal appends '_';
    nameRoot = nameRoot{1};
    imMaxNames = cell(1,imProps.nChns);
    projectionExist = false(1,imProps.nChns);
    projectionExist(chnsToLoad) = false;
    for c=1:imProps.nChns  % always check all files
        maxName = regexprep(imProps.daxName,nameRoot,['Max',imProps.chnNames{c}]);
        imMaxNames{c} = [maxProjectionFolder,maxName,'.dax'];
        if exist(imMaxNames{c},'file')
            projectionExist(c) = true;
        end
    end


    % load max projections if requested, and they exist and 
    if pars.maxProject && all(projectionExist) && ~pars.overwrite
        chnStk = zeros(imProps.imSize(2),imProps.imSize(1),length(chnsToLoad),'uint16');
        k=0;
       for c=chnsToLoad  % only load requested files for speed
           k=k+1;
           chnStk(:,:,k) = ReadDax(imMaxNames{c},'verbose',pars.verbose);
       end
       imProps.chnsLoaded = size(chnStk,3);
       imProps.projection = true;
    else
        % else load data.  
        [dax,infoFile] = ReadDax(daxFileName,'verbose',pars.verbose);
        % remove buffer frames
        dax = dax(:,:,imProps.preStartFrames+1:end-imProps.postImageFrames+pars.skipFirst);
        % update start frame based on buffer frames.  Begin at the first
        % non-buffer frame that also has 
        chnShift = rem(imProps.preStartFrames,imProps.nChns);
        if chnShift~=0
            imMaxNames = circshift(imMaxNames,-chnShift);
            chnsAll = 1:imProps.nChns;
            chnsAll = circshift(chnsAll,chnShift);
        end
        % parse channels
        chnStk = zeros(imProps.imSize(2),imProps.imSize(1),imProps.zSteps,imProps.nChns,'uint16');
        for c=1:imProps.nChns
            chnStk(:,:,:,c) = dax(:,:,pars.skipFirst+c:imProps.nChns:end);   
            % Create and save max-projections for later
            %       whether or not they are requested is requested yet.
            chnMax = squeeze(max(chnStk(:,:,:,c),[],3));
            if (~projectionExist(c) || pars.overwrite)  &&  pars.writeData
                WriteDax(chnMax,'infoFile',infoFile,'saveFullName',imMaxNames{c},'confirmOverwrite',false);
            end
        end
        % return only the requested channels
        if pars.maxProject
            chnStk = squeeze(max(chnStk,[],3));
            if chnShift~=0
                chnStk = chnStk(:,:,chnsAll);
                imProps.chnsLoaded = size(chnStk,3);
            end
            chnStk = chnStk(:,:,chnsToLoad); 
            imProps.projection = true;
        else
            if chnShift~=0
                chnStk = chnStk(:,:,:,chnsAll);
                imProps.chnsLoaded = size(chnStk,4);
            end
            chnStk = chnStk(:,:,:,chnsToLoad); 
        end
    end

    % drift correct if requested
    if pars.hyb==0 % use guess
    h = imProps.Hyb;
    else
        h=pars.hyb; % use passed value
    end
    if ~isnan(h) && ~isempty(pars.driftFolder)
        driftFileOld = [pars.driftFolder,'fov',num2str(imProps.fov+1,'%03d'),'_regData.csv'];
        driftFileNew = [pars.driftFolder,'alignTable_fov',num2str(imProps.fov+1,'%03d'),'.csv'];
        if exist(driftFileNew,'file') 
            alignTable = readtable(driftFileNew);
            fidAlign = table2struct(alignTable(alignTable.hyb==h,:));
        elseif exist(driftFileOld,'file')
            alignTable = readtable(driftFileOld);
            fidAlign = table2struct(alignTable(h,:));
        else 
            fidAlign = [];
        end
        if ~isempty(fidAlign)
            if length(size(chnStk)) == 4
                for d4=1:size(chnStk,4)                
                    for d3=1:size(chnStk,3)
                        chnStk(:,:,d3,d4) = ApplyReg(chnStk(:,:,d3,d4),fidAlign);
                    end
                end
            elseif length(size(chnStk)) == 3
                for d3=1:size(chnStk,3)
                    chnStk(:,:,d3) = ApplyReg(chnStk(:,:,d3),fidAlign);
                end
            else
                chnStk = ApplyReg(chnStk,fidAlign);
            end    
        end
    end

    % figure(43); clf;
    % for c=1:size(chnStk,3)
    %     subplot(1,size(chnStk,3),c); imagesc(chnStk(:,:,c));
    % end

    % optionally, export all data if requested
    if nargout ==3
        imData.name = daxFileName;
        imData.xml  = xmlStruct;
        imData.infoFile = infoFile;
    end
end