function [spotData,imProps] = LoadDaxFast(daxFileName,spots,varargin)
%% Overview
%  Return a cell array of 4D images, height x width x z-frames x c-colors
%     each image is centered at the x,y position provided by the input
%     array 'spots';
% The aim of this version is to only read a subset of the data
% 
% 
% Note on data organization:
% images were produced by sequentially scanning along the rows, saving one
% line at a time. After all lines in a frame were saved, the laser color
% changed and another frame was collected line by line. After this, the
% z-position was changed and the process repeats.  LoadDaxSpots will unpack
% the data in the same way.
% An xml file accompanies the binary image file and records image metadata,
% including the frame dimensions and number of channels (needed to load the
% data), and key experimental properties such as the physcial pixel size. 
% 
%%  Required inputs: 
% daxFileName - string, provides full file path to a dax file
% spots - N x 2 matrix, provides the x,y positions to crop from the image
%
%%  Outputs:
% spotData - N x 1 cell array, of h x w x z x c 4D images 
% imProps - structure recording image metadata 
% 
%% Optional inputs
% driftCorrectStruct
% 


% TO ADD
% select channels
% remove z-offset
% apply linear drift correction during load by shifting pointers


% daxFileName = 'K:\2021-04-26_GemininIF_HoxA_3Mb\Hyb_001\ConvZscan_01.dax';

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
defaults = cell(0,3);
% crop parameters
defaults(end+1,:) = {'driftCorrectStruct','struct',[]}; % specify x,y,translation
defaults(end+1,:) = {'boxRadius', 'integer',10};
defaults(end+1,:) = {'boxRadiusZ','integer',10};
% load dax parameters
defaults(end+1,:) = {'verbose', 'boolean',true};
defaults(end+1,:) = {'overwrite', 'boolean',false};
defaults(end+1,:) = {'maxProject', 'boolean',false};% not implemented yet in LoadDaxFast   
defaults(end+1,:) = {'channel', 'integer',0};  % not implemented yet in LoadDaxFast  
defaults(end+1,:) = {'skipFirst', 'integer',0};% not implemented yet in LoadDaxFast   
defaults(end+1,:) = {'writeData', 'boolean',true};  %   
pars = ParseVariableArguments(varargin,defaults,mfilename); 
% pars = ParseVariableArguments([],defaults,mfilename); 


%% Parse image metadata
xmlStruct = ReadXML(regexprep(daxFileName,'.dax','.xml'));
[imProps.dataFolder,imProps.daxName] = fileparts(daxFileName);
imProps.nFrames  =       xmlStruct.settings.acquisition.number_frames;
imProps.stageXY  =       str2num(xmlStruct.settings.acquisition.stage_position); %#ok<ST2NM>
imProps.imSize =         [xmlStruct.settings.camera1.x_pixels,xmlStruct.settings.camera1.y_pixels];
try % software_z_scan is a new parameter
    imProps.nChns =          xmlStruct.settings.focuslock.software_z_scan.frames_to_pause;
    imProps.preStartFrames = xmlStruct.settings.focuslock.software_z_scan.deadtime; 
    imProps.zStepSize =      xmlStruct.settings.focuslock.software_z_scan.step_size;  % step size in nm, as recorded by hal
    imProps.zSteps =         2*xmlStruct.settings.focuslock.software_z_scan.range/imProps.zStepSize;  % total scan is +/- scan range  
catch % handle the now py27 version of image xml 
    imProps.nFrames  =       eval(xmlStruct.settings.acquisition.number_frames);
    imProps.nChns =          xmlStruct.settings.focuslock.zscan_frames_to_pause; 
    imProps.preStartFrames = 10; 
    imProps.zStepSize =      xmlStruct.settings.focuslock.zscan_step*1e3;  % step size in nm, as recorded by hal
    imProps.zSteps =         2*xmlStruct.settings.focuslock.zscan_stop*1e3/imProps.zStepSize;  % total scan is +/- scan range  
end
imProps.postImageFrames= imProps.nFrames -imProps.preStartFrames - imProps.zSteps*imProps.nChns;
% auto-detect pixel size
objFields = fields(xmlStruct.settings.mosaic);
objUsed = strcmp(objFields,xmlStruct.settings.mosaic.objective);
magData = strsplit(xmlStruct.settings.mosaic.(objFields{objUsed}),',');
imProps.xy2um = str2double(magData(2));  % xy-pixel size in um
imProps.z2um = imProps.zStepSize/1e3;  % z-pixel size in um 

% auto-detect FOV from the append format on the daxName
nameParts = strsplit(imProps.daxName,'_');
imProps.fov = str2double(nameParts{end});
% auto-detect Hyb number based on the folder name. 
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
shutter = regexprep(shutter,{'shutters_','_series'},{'',''});
imProps.chnNames = strsplit(shutter,'_');
imProps.projection = false;

% a little short hand 
W = imProps.imSize(1);
H = imProps.imSize(2);
Z = imProps.zSteps;

%% Now load the image
N = size(spots,1); 
spotData = cell(N,1);
for s=1:N
    
    % handle drift if requested
    % drift was corrected in a 2 step algorithm.  As the algorithm allows
    % rotation and dilation at each step the translation components xshift
    % and xshift2 are reported separately rather than a single translation,
    % since translation rotation dilation are non-commuting operators. 
    % 
    % here we only correct drift to the nearest pixel.  Sub-pixel drift
    % correction is typically computed using the loaded images and added
    % into the fit. 
    if ~isempty(pars.driftCorrectStruct)
        spots(s,1) = spots(s,1) + pars.driftCorrectStruct.xshift + pars.driftCorrectStruct.xshift2;
        spots(s,2) = spots(s,2) + pars.driftCorrectStruct.yshift + pars.driftCorrectStruct.yshift2;
    end
    % parse out the ROI.  A little shorthand here makes the code below more
    % readable. 
    x1 = spots(s,1)-pars.boxRadius;
    x2 = spots(s,1)+pars.boxRadius;
    y1 = spots(s,2)-pars.boxRadius;
    y2 = spots(s,2)+pars.boxRadius;
    z1 = spots(s,3)-pars.boxRadiusZ;
    z2 = spots(s,3)+pars.boxRadiusZ;
    x1 = max(1,x1);
    x2 = min(W,x2);
    y1 = max(1,y1);
    y2 = min(H,y2);
    z1 = max(1,z1);
    z2 = min(Z,z2);
    
    %%
    binaryFormat = 'l';

    nChns = imProps.nChns;

    fid = fopen(daxFileName);
    if fid < 0
        error(['Invalid file: ' daxFileName]);
    end
    Hout = y2-y1+1;
    Wout = x2-x1+1;
    Zout = z2-z1+1;

    image4D = zeros(Hout,Wout,Zout,nChns,'uint16');
    startFrame=1;
    % read in one ine at a time, skip to the next 
    bits = 16/8; % bits/(bytes per bit) 
    stripSize = Wout;
    skipNframes = 0;
    fseek(fid,0,'bof');  
    % probably would have been faster not to return to the beginning of
    % file each time we switch spot. The loop over spots should really be
    % inside the loop over rows, with the spots pre-sorted by row.
    bytesToSkip = imProps.preStartFrames*H*W ;
    for z=z1:z2
        for c=1:nChns  % should let user select channels
            bytesToSkip = bytesToSkip + y1*W;
            for h=1:Hout % -1
                bytesToSkip = bytesToSkip +x1-1;
                fseek(fid,bytesToSkip*bits,'cof');
                lineIn = fread(fid, stripSize, '*uint16', binaryFormat);
                bytesToSkip = (W-x2); % reset bytesToSkip after readline to be end of line  
                image4D(h,:,z,c) = lineIn;
            end
            % bytesToSkip = bytesToSkip + (H-y2)*W; % skip rest of y in this frame
            bytesToSkip = bytesToSkip + (H-y2-1)*W; % skip rest of y in this frame
            if skipNframes % not currently used. could skip every other frame 
                bytesToSkip = bytesToSkip + skipNframes*H*W;
            end
        end
    end
    fclose(fid);
    spotData{s} = image4D;
end

% 
% figure(1); clf;
% ProjectIm3D(image4D(:,:,:,1));
