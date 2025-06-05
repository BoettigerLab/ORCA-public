function [imH,imMax,imOut] = ShowDaxSeries(daxFiles,varargin)
%% Inputs:
% takes a list of daxfiles, load them as a multi-color cell array
% 
%% Notes:
%   defaults to 2 channels, display only 1st channel 
%    as this is our most common data type (data + fiducial). 
% 
%  Updated 5/24/24 to allow multiple data channels
% 
%%

defaults = cell(0,3);
defaults(end+1,:) = {'daxType',{'maxProject','Zstack'},'Zstack'};
defaults(end+1,:) = {'analysisFolder','string',''}; % leave blank to save in [dataFolder,'Analysis\']; 
defaults(end+1,:) = {'numChannels','integer',2};
defaults(end+1,:) = {'driftPars','struct',[]};
defaults(end+1,:) = {'driftRefHyb','integer',1};
defaults(end+1,:) = {'selectChannels','integer',1};
defaults(end+1,:) = {'selectHybes','integer',0};
defaults(end+1,:) = {'framesPerTime','integer',0}; % load all
defaults(end+1,:) = {'channelOrder',{'alternate','sequential'},'alternate'};
defaults(end+1,:) = {'figHandle','handle',200};
defaults(end+1,:) = {'childHandle','handle',200};
defaults(end+1,:) = {'autoContrastMin','fraction',.3};
defaults(end+1,:) = {'autoContrastMax','fraction',.9999};
defaults(end+1,:) = {'alignContrastHigh','fraction',1};
defaults(end+1,:) = {'alignContrastLow','fraction',0};
defaults(end+1,:) = {'currTime','integer',0};  % current time to show
defaults(end+1,:) = {'showOverlay','boolean',true};
defaults(end+1,:) = {'xylim','array',[]};
defaults(end+1,:) = {'channelNameXls','string',''};
defaults(end+1,:) = {'driftCorrectUsingChn','integer',0}; % 0 is don't run
defaults(end+1,:) = {'overwriteDrift','boolean',false}; 
defaults(end+1,:) = {'offsetFrame','integer',0}; 
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename); 
   
%% Main function

% only use a subset of autodetected Hybs 
%   (maybe we should require this is done 1 level up, before calling)
if pars.selectHybes ~=0
    daxFiles = daxFiles(pars.selectHybes);
end

% if using drift correction we need a place to save files 
if pars.driftCorrectUsingChn
    % navigate up 1 directory from the hyb folders, and create an analysis
    % folder (if it does not exist).  This is where we will save the
    % registration information. 
    [dataFolder,fname] = fileparts(daxFiles{1});
    % seps = strfind(dataFolder,'\');  % 
    % hybNames = cellfun(@(x) x(seps(end)+1:length(dataFolder)),daxFiles,'UniformOutput',false);
    %   % the above only works if all the folder names are the same length
    nameParts = cellfun(@(x) strsplit(x,filesep), daxFiles,'UniformOutput',false);
    nameParts = cat(1,nameParts{:}); % this assumes all the daxFiles have the same degree of nesting.  This should rarely be a problem because they are often all in the same folder. Or in parallel RNA vs. DNA folders 
    hybNames = nameParts(:,end-1);

%     dataFolder = dataFolder(1:seps(end));
    if isempty(pars.analysisFolder)
        analysisFolder = [dataFolder,filesep,'Analysis\'];
    else
        analysisFolder = pars.analysisFolder;
    end
    if ~exist(analysisFolder,'dir')
        mkdir(analysisFolder)
    end
    
    alignTableOld = false;
    % we also get the FOV number from the file name 
    % and create a alignment file in this folder
    fnum = strsplit(fname,'_');
    f = str2double(fnum{end});
    alignTableFile = [analysisFolder,'alignTable_fov',num2str(f,'%03d'),'.csv'];
    alignTableFileOld = [analysisFolder,'fov',num2str(f+1,'%03d'),'_regData.csv'];
    if exist(alignTableFileOld,'file') && ~exist(alignTableFile,'file')
        if pars.overwriteDrift
            delete(alignTableFileOld);
        else
            alignDriftTable = readtable(alignTableFileOld);
            alignTableOld = true;
        end
    end

    if exist(alignTableFile,'file')
        if pars.overwriteDrift
            delete(alignTableFile);
        else
            alignDriftTable = readtable(alignTableFile);
        end
    end
end


try
%% now we acutally load the data
numHybes = length(daxFiles);
numSelectChns = length(pars.selectChannels);
imH = cell(numHybes,numSelectChns);

hs = 1:numHybes;
hs(pars.driftRefHyb) = [];
hs = [pars.driftRefHyb,hs]; % make the drift reference the first hyb  processed.   (1 by default) 

% not needed
% origOrder = [1+hs(2:pars.driftRefHyb),1,hs(pars.driftRefHyb+1:numHybes)];
% origOrder(hs)

hh = 0; % counter of the hybes we've processed
for h=hs
    hh=hh+1;
    if strcmp(pars.daxType,'Zstack')
        im = cell(pars.numChannels,1);
        dax = ReadDax(daxFiles{h},'verbose',false);
        if  strcmp(pars.channelOrder,'alternate')
            for c=1:pars.numChannels
                im{c} = dax(:,:,c+pars.offsetFrame:pars.numChannels:end);
            end
        end
    end
    
    if pars.selectChannels ~= 0 
        if pars.driftCorrectUsingChn
            imF = im{pars.driftCorrectUsingChn};
            if h==1  % align to hyb1
                refHybFid = max(imF,[],3);
            end
        end
        % im = im(pars.selectChannels);        
    end

    for c = pars.selectChannels
        imH{h,c} = im{c}; 
    end

 % as currently written will only load 1 select channel as overlay
        %  I think this is generally convienent for addressing color-specific
        %  background and color-specific aberrations

    if pars.driftCorrectUsingChn > 0
        alignData = [];
        if exist('alignDriftTable','var')
            try % just because the variable exists, we might not have processed this hyb yet. if not, we just reload. 
                if alignTableOld
                    if hh==1
                        alignData.xshift = 0;
                        alignData.yshift = 0;
                        alignData.theta = 0;
                        alignData.xshift2 = 0;
                        alignData.yshift2 = 0;
                        alignData.theta2 = 0;
                    else
                        alignData = alignDriftTable(h,:);
                    end
                else                
                    alignData = alignDriftTable(strcmp(alignDriftTable.hyb,hybNames{h}),:);  % we record the hyb name we aligned too
                    if height(alignData)>1
                        warning(['found multiple entries for ',hybNames{h}]);
                        alignData = alignData(1,:);
                    end
                end
            catch
            end
        end
        if isempty(alignData)
            if hh==1 
                refHybFid = max(imF,[],3);
                curHybFid = max(imF,[],3);
            else
                curHybFid = max(imF,[],3);
            end
            if pars.alignContrastHigh ~= 1
                refHybFid = IncreaseContrast(refHybFid,'high',pars.alignContrastHigh,'low',pars.alignContrastLow);
                curHybFid = IncreaseContrast(curHybFid,'high',pars.alignContrastHigh,'low',pars.alignContrastLow);
            end
            
            f1 = figure(30); clf; 
            if isempty(pars.driftPars)
                fidAlign = CorrAlignFast(refHybFid,curHybFid);  
            else
                fidAlign = CorrAlignFast(refHybFid,curHybFid,'parameters',pars.driftPars);  
            end
            fidAlign.fov = f;
            fidAlign.hyb = hybNames{h};
            alignTable = struct2table(fidAlign);
            writetable(alignTable,alignTableFile,'WriteMode','Append'); 
            if h==hs(1) && pars.verbose
                disp(['started writing ',alignTableFile])
            end
            if h==max(hs) && pars.verbose
                disp(['completed writing ',alignTableFile])
            end
            SetFigureSavePath([analysisFolder,'CorrAlign',filesep],'makeDir',true);
            imName = ['fov',num2str(f,'%03d'),'_H',num2str(h,'%03d')];
            SaveFigure(f1,'name',imName,'formats',{'png'},'overwrite',true);
            alignData = alignTable(strcmp(alignTable.hyb,hybNames{h}),:);
        end
        for c = pars.selectChannels
            imH{h,c} = ApplyReg(imH{h,c},alignData); 
        end
    end
end

catch er
    warning(er.getReport);
    disp('debug here.')
end


%%
% blank or incomplete dax files can create problems 
%  we introduce some error handling for this. 
try
    imMax_c = cell(numSelectChns,1);
    imOut_c = cell(numSelectChns,1);
    for c=pars.selectChannels
        try
        [imMax_c{c},imOut_c{c}] =  ShowXY(imH(:,c),'parameters',pars);
        catch er
            warning(er.getReport);
            disp('debug here.')
        end
    end
    hasData = ~cellfun(@isempty,imMax_c);
    imMax = cat(4,imMax_c{hasData});
    imOut = cat(4,imOut_c{hasData});
catch er
    warning(er.getReport);
    disp('debug here.')
end
