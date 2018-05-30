function [movie, infoFile, infoFileRoi,memoryMap,memMapToFrame] = ReadDax(varargin)
%--------------------------------------------------------------------------
% [movie, infoFiles] = ReadDax(fileName, varargin)
% This function loads a STORM movies from the dax file associated with the
% provided .inf file
%--------------------------------------------------------------------------
% Outputs:
% movies/LxMxN array: A 3D array containing the specified movie
% infoFile: infoFile structure for the specified daxfile
% infoFileRoi: modified infoFile corresponding to the daxfile
%--------------------------------------------------------------------------
% Input:
% fileName/string or structure: Either a path to the dax or inf file or an
%   infoFile structure specifying the dax file to load
%
%--------------------------------------------------------------------------
% Variable Inputs:
%
% 'file'/string ([]): A path to the associated .inf file
%
% 'path'/string ([]): Default path to look for .inf files
%
% 'startFrame'/double  (1): first of movie to load.  
%
% 'endFrame'/double ([]): last frame of the movie to load.  If empty will
% be max. 
%
% 'subregion'/double (zeros(4,1)):  [xmin, xmax, ymin, ymax] of the region
% of the dax file to load.  Pixels indexed from upper left, as in images.  
%
% 'infoFile'/info file structure ([]): An info file for
% the files to be loaded.  
%
% 'imageDimensions'/3x1 integer array ([]): The size of the movies to be
% loaded.  
%
% 'verbose'/boolean (true): Display or hide function progress
%
% 'orientation'/string ('normal'): Control the relative rotation of the data
%   structure
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% September 7, 2012
%
% Version 1.1
%-------------------Updates:
% 01/19/13: ANB
% modified to allow arbitrary start and end frame to be specified
% by the user.  Removed 'image_dimension' flag (this was non-functional)
% and removed allFrames (this has become redundant);  
%-----------------------
% 2/14/13: JRM
% Minor fix to dax data type
%-----------------------
% ~12/15/13: ANB
% ReadDax now respects binning options in dax file
% ReadDax also computes how much memory it will take to load the requested
% file and throws a warning if this exceeds a certain max value. Default
% max is 1 Gb.  Warning allows user to continue, reduce frames, or abort.
%-----------------------
% 12/22/13: ANB
% Added 'subregion' feature.  ReadDax now uses memmaps
% -------------------
% ~08/01/15: ANB
% fixed bug: data-type was hard-coded, should use what the info file
% specifies. 
% ---------------------
% 09/20/17: ANB
% converted to modern default parameter parsing
% 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;


%--------------------------------------------------------------------------
% Parse Required Input
%--------------------------------------------------------------------------

if isstruct(varargin{1})
    infoFile = varargin{1};
    fileName = [];
else
    infoFile = [];
    fileName = varargin{1};
end
varin = varargin(2:end);

%--------------------------------------------------------------------------
% Parse Optional Variables
%--------------------------------------------------------------------------


defaults = cell(0,3);
defaults(end+1,:) = {'fileName','string',fileName};
defaults(end+1,:) = {'startFrame','freeType',1};
defaults(end+1,:) = {'endFrame','freeType',[]};
defaults(end+1,:) = {'subregion','array',[0,0,0,0]};
defaults(end+1,:) = {'allFrames','boolean',[]};
defaults(end+1,:) = {'infoFile','struct',infoFile};
defaults(end+1,:) = {'maxMemory','positive',10E9};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'orientation',{'normal','nd2'},'normal'};
defaults(end+1,:) = {'dataPath','string',defaultDataPath};

pars = ParseVariableArguments(varin,defaults,mfilename);
fileName = pars.fileName;
startFrame = pars.startFrame;
endFrame = pars.endFrame;
subregion = pars.subregion;
allFrames = pars.allFrames;
infoFile = pars.infoFile;
maxMemory = pars.maxMemory;
verbose = pars.verbose;
orientation = pars.orientation;
dataPath = pars.dataPath;

%--------------------------------------------------------------------------
% Check parameter consistency
%--------------------------------------------------------------------------
if ~isempty(infoFile) && ~isempty(fileName)
    error('You cannot specify info files and file names');
end

%--------------------------------------------------------------------------
% Load info files if needed
%--------------------------------------------------------------------------
if isempty(infoFile)
    if isempty(fileName)
        infoFile = ReadInfoFile('path', dataPath, 'verbose', verbose);
    else
        infoFile = ReadInfoFile(fileName, 'verbose', verbose);
    end
    
    if isempty(infoFile)
        disp('Canceled');
        movie = {};
        return;
    end 
end

%--------------------------------------------------------------------------
% Load Dax Files
%--------------------------------------------------------------------------

% --------- Determine number of frames to load
framesInDax = infoFile.number_of_frames;

% parse now outdated 'allFrames' for backwards compatability
if ~isempty(allFrames) 
    if allFrames
        endFrame = framesInDax;
    else 
        endFrame = 1;
        startFrame = 1; 
    end
end

% parse startFrame and endFrame
if isempty(endFrame)
    endFrame = framesInDax;
end
if endFrame > framesInDax
    if verbose
        warning('input endFrame greater than total frames in dax_file.  Using all available frames after startFrame');
    end
    endFrame = framesInDax; 
end
framesToLoad = endFrame - startFrame + 1;
frameDim = [infoFile.frame_dimensions(1)/infoFile.binning(1),...
            infoFile.frame_dimensions(2)/infoFile.binning(2)];
frameSize = frameDim(1)*frameDim(2);

% Check memory requirements.  Ask for confirmation if > maxMemory.  
memoryRequired = frameSize*framesToLoad*16/8;
DoThis = 1; 
if memoryRequired > maxMemory
    warning([fileName,'  ', 'is > ',num2str(maxMemory/1E6,3),' Mbs.']);
    DoThis = input(['Requested file requires ',...  
        num2str(memoryRequired/1E6,3),' Mbs. ',...
        'Are you sure you want to load it? ',... 
        '(Filling memory may crash the computer) ',...
        '0 = abort, 1 = continue, n = new end frame  ']);
    if DoThis > 1 
        endFrame = DoThis; % read in the new last frame
        framesToLoad = endFrame - startFrame + 1;
        DoThis = true;
    end
end
  
%--------------------------------------------------------
% Proceed to load specified poriton of daxfile
%--------------------------------------------------------
if DoThis
    fileName = [infoFile.localName(1:(end-4)) '.dax'];
    if verbose
        disp(['Loading ' infoFile.localPath fileName ]);
    end

    if ~isempty( strfind(infoFile.data_type,'little endian') )
        binaryFormat = 'l';
    else
        binaryFormat = 'b';
    end
    
    %----------------------------------------------------------------- 
    % Read all pixels from selected frames
    %-----------------------------------------------------------------
    if sum(subregion)==0  && nargout < 3   
        fid = fopen([infoFile.localPath fileName]);
        if fid < 0
            error(['Invalid file: ' infoFile.localPath fileName]);
        end

        fseek(fid,(frameSize*(startFrame - 1))*16/8,'bof'); % bits/(bytes per bit) 
        dataSize = frameSize*framesToLoad;
        movie = fread(fid, dataSize, '*uint16', binaryFormat);
        fclose(fid);

        try % Catch corrupt files
            if framesToLoad == 1
                movie = reshape(movie, frameDim)';
            else
                switch orientation % Change orientation
                    case 'normal'
                        movie = permute(reshape(movie, [frameDim framesToLoad]), [2 1 3]);
                    case 'nd2'
                        movie = permute(reshape(movie, [fliplr(frameDim) framesToLoad]), [2 1 3]);
                    otherwise
                end
            end
        catch
            disp('Serious error somewhere here...check file for corruption');
            movie = zeros(frameDim);
        end
        infoFileRoi = infoFile; 
    %-----------------------------------------------------------------    
    % Read the Indicated SubRegion using Memory Maps 
    %-----------------------------------------------------------------
    else
        % parse short-hand: xmin = 0 will start at extreme left
        %               ymax = 0 will go from ymin to the bottom
        xi = uint32(subregion(1));
        xe = uint32(subregion(2));
        yi = uint32(subregion(3));
        ye = uint32(subregion(4));
        if xi == 0 
            xi = uint32(1);
        end
        if xe == 0
            xe = uint32(frameDim(1));
        end
        if yi == 0
            yi = uint32(1);
        end
        if ye == 0
            ye = uint32(frameDim(2));
        end
        %------------------ arbitrary region ------------------------
        if verbose
            disp('creating memorymap');
        end
        memoryMap = memmapfile([infoFile.localPath fileName], ...
                'Format', 'uint16', ...
                'Writable', false, ...
                'Offset', (startFrame-1)*frameSize*16/8, ...
                'Repeat', framesToLoad*frameSize);  
  
        [ri,ci,zi] = meshgrid(xi:xe,yi:ye,uint32(1):uint32(framesToLoad));
        inds = sub2indFast([frameDim(1),frameDim(2),framesInDax],...
                        ri(:),ci(:),zi(:));
        movie = memoryMap.Data(inds); 
        movie = swapbytes(movie);
        xs = xe-xi+uint32(1);
        ys = ye-yi+uint32(1);
        movie = reshape(movie,[xs,ys,framesToLoad]);
        memMapToFrame = @(frames) reshape(frames(inds), [xs,ys,framesToLoad]);
        if ~strcmp(orientation,'normal')
         movie = permute(reshape(movie, [xs,ys,framesToLoad]), [2 1 3]);
        end
        infoFileRoi = infoFile; 
        infoFileRoi.hend = xs;
        infoFileRoi.vend = ys;
        infoFileRoi.frame_dimensions = [infoFile.hend,infoFile.vend];
        infoFileRoi.file = [infoFile.localPath,infoFile.localName(1:end-4),'.dax'];
        %--------------------------------------------------
    end
    

    if verbose
        disp(['Loaded ' infoFile.localPath fileName ]);
        disp([num2str(framesToLoad) ' ' num2str(frameDim(1)) ' x ' num2str(frameDim(2)) ...
            ' frames loaded']);
    end
else
    error('User aborted load dax due to memory considerations '); 
end