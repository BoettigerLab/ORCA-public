function [im,imMax,imOut,imMaxZ,imOutZ,daxData] = ShowDax(daxFile,varargin)
% a daxReader for my daxViewer
% 
% Can load simple frames, coupled to slider (see DaxApp)
% Can load multi-color images as overlays from alternating-channel data
% Can load 2-cam data as registered overlays
% Can load Hyb series as registered overlays 
% All overlays are shown as axes-linked subpanels as well.
% Can load 3D data as 3D stacks and show the XY and XZ
% Can load continuous time z-scans
% 


defaults = cell(0,3);
defaults(end+1,:) = {'daxType',{'Zstack','2cam'},'Zstack'};  %   ,'Movie'
defaults(end+1,:) = {'numChannels','integer',1};
defaults(end+1,:) = {'selectChannels','integer',0};
defaults(end+1,:) = {'framesPerTime','integer',0}; % load all
defaults(end+1,:) = {'channelOrder',{'alternate','sequential'},'alternate'};
defaults(end+1,:) = {'figHandle','handle',100};
defaults(end+1,:) = {'childHandle','handle',200};
defaults(end+1,:) = {'autoContrastMin','fraction',.3};
defaults(end+1,:) = {'autoContrastMax','fraction',.9999};
defaults(end+1,:) = {'currTime','integer',0};  % current time to show
defaults(end+1,:) = {'showOverlay','boolean',true};
defaults(end+1,:) = {'showXZ','boolean',false};
defaults(end+1,:) = {'xylim','array',[]};
defaults(end+1,:) = {'alignment_file','string','F:\DefaultPars\alignmentData.txt'}; % needed for 2cam files
defaults(end+1,:) = {'chrom_correct_file','string','F:\DefaultPars\tform3D.mat'}; % needed for 2cam files

defaults(end+1,:) = {'fitFile','cell',''}; % cell array of strings
defaults(end+1,:) = {'colormap','colormap',gray(256)};
defaults(end+1,:) = {'colorbar','boolean','true'};
pars = ParseVariableArguments(varargin,defaults,mfilename); 
   
if exist(pars.alignment_file,'file')
    alignC = readtable(pars.alignment_file);
end

% folderCal2 = 'I:\Jude\2023-10-10_400kb_test\';
% daxFile = [folderCal2,'beads_0001_C1.dax'];
% pars.numChannels =2; 

% 
% pars.currTime  = 1;
% pars.framesPerTime = 10;
daxData.infoStruct = [];
daxData.xmlStruct = [];
daxData.offsetTable = [];
im=[];
imMax=[];
imOut=[];
imMaxZ=[];
imOutZ=[];
%% main function
if strcmp(pars.daxType,'Zstack') ||  strcmp(pars.daxType,'2cam')
    im = cell(pars.numChannels,1);
    % load info file and dax file and xml file 
    % info file
    infoFile = ReadInfoFile(daxFile); 
    daxData.infoStruct = infoFile;
    if nargout ==6
        % offset file
        offsetFile = regexprep(daxFile,'.dax','.off');
        offsetFile = regexprep(offsetFile,'_C1',''); % strip camera ID if present (the offset value is for all cameras) 
        daxData.offsetTable = ReadTableFile(offsetFile,'delimiter',' ');
        % read xml file with scope parameters 
        xmlFile = regexprep(daxFile,'.dax','.xml');
        xmlFile = regexprep(xmlFile,'_C1','');
        if exist(xmlFile,'file')
            daxData.xmlStruct = ReadXML(xmlFile);
        end
    end

    maxFrame = infoFile.number_of_frames;
    startFrame = pars.currTime*pars.framesPerTime+1;
    endFrame = min( (pars.currTime+1)*pars.framesPerTime, maxFrame);
    if pars.framesPerTime == 0  % load all frames of the z-stack
        dax = ReadDax(daxFile,'verbose',false);
    elseif pars.framesPerTime > 0  % load N frames of a z-stack from a single time window
        dax = ReadDax(daxFile,'startFrame',startFrame,'endFrame',endFrame,'verbose',false); 
    end
    if  strcmp(pars.channelOrder,'alternate')
        for c=1:pars.numChannels
            im{c} = dax(:,:,c:pars.numChannels:end);
        end
    end
    % remove unselected channels 
    %    (useful to discard fiducial frames if they are unused)
    %    This is slightly slower than specific loading but more flexible.
    if pars.selectChannels ~= 0 
        im = im(pars.selectChannels);
    end

    if strcmp(pars.daxType,'2cam')
        daxFile2 = regexprep(daxFile,'C1','C2');
        im2 = cell(pars.numChannels,1);
        if pars.framesPerTime == 0  % load all frames of the z-stack
            dax = ReadDax(daxFile2,'verbose',false);
        elseif pars.framesPerTime > 0  % load N frames of a z-stack from a single time window
            dax = ReadDax(daxFile2,'startFrame',startFrame,'endFrame',endFrame,'verbose',false); 
        end
        % Apply camera registreation
        dax = ApplyReg(fliplr(dax),alignC,'invert',false);
        if  strcmp(pars.channelOrder,'alternate')
            for c=1:pars.numChannels
                im2{c} = dax(:,:,c:pars.numChannels:end);
            end
        end
        if pars.selectChannels ~= 0 
            im2 = im2(pars.selectChannels);
        end
        % merge channels
        im = cat(1,im,im2);
    end
end

if ~isempty(pars.fitFile)
    nFits = length(pars.fitFile);
    pars.spotXY = cell(nFits,1);
    for c=1:nFits
        data = LoadHD5Fits(pars.fitFile{c},'parameters',pars,'frames',startFrame:endFrame);
        if c==2 && strcmp(pars.daxType,'2cam')
            data = Register2CamFits(data,'alignment_file',pars.alignment_file,'chrom_correct_file',pars.chrom_correct_file);
        end
        xx = cat(1,data.x) + 1;
        yy = cat(1,data.y) + 1;
        hold on; plot(xx,yy,'ro','MarkerSize',16);
        pars.spotXY{c} = [xx,yy];
    end
end

%%
[imMax,imOut] =  ShowXY(im,'parameters',pars);

if pars.showXZ
    [imMaxZ,imOutZ] =  ShowXZ(im,'parameters',pars);
end





