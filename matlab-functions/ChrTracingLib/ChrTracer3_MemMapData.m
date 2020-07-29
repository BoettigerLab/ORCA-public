function [fidMapAlign,datMapAlign,fidFlatFrames,datFlatFrames] = ChrTracer3_MemMapData(fileNames,varargin)
% build memory maps for all data
% This is done only for the single field of view, passed in fidNames and
% datNames. It allows for rapid and memory efficient spot fitting of the
% data later on.
% Also loads max-projections of all hybes for fiduical and/or data if
% requested. If these files have already been created (SaveRegData will
% offer to do this), we load them instead for speed. This should also
% accelerate execution.  If they haven't been created, this function will
% create them. 
% 

defaults = cell(0,3);
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'fov','integer',1}; % dummy parameter, processed in ChrTracer3  
defaults(end+1,:) = {'loadFlatData','boolean',true};
defaults(end+1,:) = {'loadFlatFid','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'saveMaxProject','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);



% ------------ Create new memory maps of saved AlignedDax -----------------
tic
numHybes = size(fileNames.datNames,1);
numDataChns = size(fileNames.datNames,2);

fidMapAlign = cell(numHybes,1);
datMapAlign = cell(numHybes,numDataChns);
fidFlatFrames = cell(numHybes,1);
datFlatFrames = cell(numHybes,numDataChns); 

if pars.verbose
    disp('memory mapping data');
end
for h=1:numHybes
    if pars.verbose
        disp(['mapping data from Hyb ',num2str(h),' of ',num2str(numHybes),'...']);
    end
        
    % determine if a max-projection exists and is requested
    [~,maxName]=fileparts(fileNames.fidNames{h});
    maxName = ['max_',maxName,'.inf']; %#ok<AGROW>
    hasMax = exist([pars.saveFolder,filesep,maxName],'file')==2;
    if hasMax && pars.loadFlatFid
        loadFile = false;
        fidFlatFrames{h} = ReadDax([pars.saveFolder,maxName],'verbose',pars.veryverbose);
    elseif ~hasMax && pars.loadFlatFid
        loadFile = true;
    else
        loadFile = false;
    end
    % The main step: Map the data. If requested also load the data now to create a max projection.
    [dax,fidMapAlign{h}] = ReadData(fileNames.fidNames{h},'loadFile',loadFile,'verbose',pars.veryverbose);
    
    % if we created a max projection, let's save a copy for future use (unless requested not to).  
    if loadFile
        fidFlatFrames{h} = max(dax,[],3);
        if pars.saveMaxProject && ~hasMax % write max projection for speed. 
            infoOut = fidMapAlign{h}.infoFile;
            infoOut.number_of_frames= 1;
            infoOut.localName = maxName;
            WriteDAXFiles(fidFlatFrames{h},infoOut,'verbose',pars.veryverbose);
        end
    end
    
    % we also need to do this for the data channel
    for n=1:numDataChns
        [~,maxName]=fileparts(fileNames.datNames{h,n});
        maxName = ['max_',maxName,'.inf']; %#ok<AGROW>
        hasMax = exist([pars.saveFolder,filesep,maxName],'file')==2 ;
        if hasMax && pars.loadFlatData
            loadFile = false;
            datFlatFrames{h,n} = ReadDax([pars.saveFolder,filesep,maxName],'verbose',pars.veryverbose);
        elseif ~hasMax && pars.loadFlatData
            loadFile = true;
        else
            loadFile = false;
        end
        % main step to load data
        [dax,datMapAlign{h,n}] = ReadData(fileNames.datNames{h,n},'loadFile',loadFile,'verbose',pars.veryverbose);
        if loadFile
            datFlatFrames{h,n} = max(dax,[],3);
            if pars.saveMaxProject && ~hasMax % write max projection for speed.
                infoOut = datMapAlign{h,n}.infoFile;
                infoOut.number_of_frames= 1;
                infoOut.localName = maxName;
                WriteDAXFiles(datFlatFrames{h,n},infoOut,'verbose',pars.veryverbose);
            end
        end
    end
end
fidFlatFrames = cat(3,fidFlatFrames{:});
datFlatFrames = cat(3,datFlatFrames{:});

% % data returned by function
% dataOut.fidMapData = fidMapAlign;
% dataOut.datMapData  = datMapAlign;


t = toc;
disp(['data mapped in ',num2str(t/60),' min']);

%%
function [movie,mapData] = ReadData(fileName,varargin)
%  To update:
% I think this was mostly moved into its own function already, if not, it
% probably should live there rather than a subfunction here. It looks
% pretty general. I think it's hacked out of ReadDax's memmap approach to
% reading specified ROI's efficiently, which Alistair developed earlier. 
% 

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'startFrame','integer',1}; % 0 = first frame
defaults(end+1,:) = {'endFrame','integer',0}; % 0 = last frames
defaults(end+1,:) = {'selectFrames','integer',0}; % 0 = all frames
defaults(end+1,:) = {'roi','integer',0}; % 0 = full field
defaults(end+1,:) = {'mapName','string',''};
defaults(end+1,:) = {'loadFile','boolean',true};
pars = ParseVariableArguments(varargin,defaults,'ReadData');
% pars = ParseVariableArguments([],defaults,'ReadData');

% parse info file 
infoFile = ReadInfoFile(fileName, 'verbose', pars.verbose);
framesInDax = infoFile.number_of_frames;
frameDim = [infoFile.frame_dimensions(1)/infoFile.binning(1),...
            infoFile.frame_dimensions(2)/infoFile.binning(2)];
frameDim = uint32(frameDim);
framesInDax = uint32(framesInDax);
if ~isempty( strfind(infoFile.data_type,'little endian') )
    binaryFormat = 'l';
else
    binaryFormat = 'b';
end

if pars.roi == 0
   xi = 1; xe = frameDim(1);
   yi = 1; ye = frameDim(2);
else
    xi = pars.roi(1);
    xe = pars.roi(2);
    yi = pars.roi(3);
    ye = pars.roi(4);
end
if pars.startFrame == 0
    startFrame = 1;
else
    startFrame = pars.startFrame;
end
if pars.endFrame == 0
    endFrame = framesInDax;
else
    endFrame = pars.endFrame;
end
if sum(pars.selectFrames) == 0
    selectFrames = uint32(startFrame):uint32(endFrame);
else
    selectFrames = uint32(pars.selectFrames);
end
framesToLoad = length(selectFrames);


% main function
memoryMap = memmapfile(fileName, ...
            'Format', 'uint16', ...
            'Writable', false, ...
            'Offset', 0); %  (startFrame-1)*frameSize*16/8);  

if pars.loadFile
    [ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);
    inds = sub2indFast([frameDim(2),frameDim(1),framesInDax],ri(:),ci(:),zi(:));
    movie = memoryMap.Data(inds);
    if strcmp(binaryFormat,'b')
        movie = swapbytes(movie);
    end
    xs = xe-xi+uint32(1);
    ys = ye-yi+uint32(1);
    movie = reshape(movie,[ys,xs,framesToLoad]); % figure(7); clf; imagesc(max(movie,[],3)); colormap gray;
else
    movie = [];
end 
mapData.memMap = memoryMap;
mapData.frameDim =frameDim;
mapData.framesInDax = framesInDax;
mapData.selectFrames = selectFrames;
mapData.name = pars.mapName;   
mapData.binaryFormat = binaryFormat;
mapData.infoFile = infoFile;
