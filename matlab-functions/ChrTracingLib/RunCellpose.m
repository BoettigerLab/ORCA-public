
function cellMasks = RunCellpose(imFiles,varargin)
%% input
% imFiles can be a cell array of file paths to images, which can be dax
% files, jpg, png or tiff.
% imFiles can also be a cell array of images (2D matrices).
% imFiles can also be the path to single file (as a string) or an
% individual image. 
% 
%% Output
% if the input is a cell array, the output will be a cell array
% if the input is a single image or single filepath, the ouput will be a
% single labeled image.
% 
% This implementation wraps only the 2D function, 3D wrapper will follow
% 
% Overwrite/rerun
%    Be default it saves outputs in the scratch folder, which are rewritten
%    everytime to avoid confustion. 
%    If a different saveFolder is given, it will load that. This option is
%    best when passing a stack of images. 
%    
%% Examples
% imFiles = FindFiles([nucleiFolder,'\',root,'*.dax']);
% RunCellpose(imFiles);
% 
%% Updates
%  11/2/2023 - updated code to read the image size of the dax file from the
%  info file.  When loading pre-saved data, this ensures the output
%  behavior of this function returns a map of the same size as when the
%  function is first called.  Note, it does assume the images are the same
%  size. 
global cellpose_env scratchPath  

defaults = cell(0,3);
defaults(end+1,:) = {'saveFolder','string',''}; % will create a default folder if blank. This is where the cellpose output and processed table will be saved.
defaults(end+1,:) = {'overwrite','boolean',false}; %
defaults(end+1,:) = {'do_3D','boolean',false}; %
defaults(end+1,:) = {'clearScratch','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'useGPU','boolean',true}; 
defaults(end+1,:) = {'veryverbose','boolean',false};
% defaults(end+1,:) = {'hyperstack','array',[]}; 
defaults(end+1,:) = {'imSize','array',[512,512]}; % for cellpose, we downsample images to this size to run faster
defaults(end+1,:) = {'figShowLoadedImages','integer',1};
defaults(end+1,:) = {'diameter','nonnegative',0}; % 0 = autodetect  a fixed diameter is substantially faster, and may improve uniformity. 
defaults(end+1,:) = {'model',{'nuclei','cyto'},'nuclei'};
defaults(end+1,:) = {'deleteWhenDone','boolean',false};
defaults(end+1,:) = {'cellpose_env','string',''};
defaults(end+1,:) = {'device',{'default','none','cuda:1','cuda:2','cuda:3','cpu'},'default'}; % cuda:1
% defaults(end+1,:) = {'env','string','mlab_cellpose'};
% defaults(end+1,:) = {'condaPrompt','string',condaPrompt};
% defaults(end+1,:) = {'activateConda','string',activateConda};
pars = ParseVariableArguments(varargin,defaults,mfilename);


if isempty(cellpose_env) % try to guess most likely place
    cellpose_env = pars.cellpose_env;
    if isempty(cellpose_env)
        warning('global cellpose_env not specified. Please specify how to activate the cellpose environment in your startup');
        warning('something like: global cellpose_env; cellpose_env = "C:\Anaconda3\Scripts\activate.bat && conda activate cellpose";')
        error('cellpose not found');
    end
end


if isempty(scratchPath) && isempty(pars.saveFolder)
    error('must either specify a global scratchPath or explicit "saveFolder" /path/to/saveData');
end
deleteWhenDone = pars.deleteWhenDone;
if isempty(pars.saveFolder)
    saveFolder = [scratchPath,'cellposeScratch\']; 
    deleteWhenDone = true; % scratch is always cleaned up;
else
    saveFolder = [pars.saveFolder,'cellpose\'];
    % note default is not to delete when done if making new folder
end

% allow single files
if ~iscell(imFiles)
    imFiles = {imFiles};
    singleFile = true;
else
    singleFile = false;
end
nFOV = length(imFiles);

% some flags
skipFitting = false;
skipPngs = false;


% create folder if it does not exist. Otherwise clear its contents
% it is essential the folder is empty for this to work.
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
else
    if pars.overwrite
        delete([saveFolder,'*']); 
    end    
    allMasks = FindFiles([saveFolder,'*_cp_masks.tif']);
    allPngs =  FindFiles([saveFolder,'image*.png']);
    if length(allMasks) == length(imFiles)
        skipFitting = true;
    end
    if length(allPngs) == length(imFiles)
        skipPngs = true;
    end
end

% get image size
if strcmp(imFiles{1}(end-2:end),'dax')
    infoData = ReadInfoFile(imFiles{1});
    imHeight = infoData.y_end - infoData.y_start +1;
    imWidth = infoData.x_end - infoData.x_start +1;
else
    [imHeight,imWidth,~] = size(imFiles{1});
end

% save files as 
 % requires a cell array of file paths. If we recieved just a string,
 % convert it to a 1-element array for processing
origSize = repmat([imHeight,imWidth,1],nFOV,1);
if ~skipPngs
    for f=1:nFOV
        if ischar(imFiles{f})
            if strcmp(imFiles{f}(end-2:end),'dax')
                im1 = ReadDax(imFiles{f},'verbose',pars.veryverbose);
            else
                im1 = imread(imFiles{f});
            end
            [~,imName] = fileparts(imFiles{f});
        else
            im1 = imFiles{f};
            imName = ['image_',num2str(f,'%03d')];
        end
        [h,w,z] = size(im1);
        origSize(f,:) = [h,w,z];
        if z>1 && ~pars.do_3D  % can add a 3D option here later, will also need to change save to write a stack. 
             
%             % Logically, I think splitting channels should be handled
%             % upstream of the 
%             if ~isempty(pars.hyperstack)
%                 if pars.hyperstack(3) == 0
%                     pars.hyperstack(3) = z;
%                 end
%                 im1 = im1(:,:,pars.hyperstack(1):pars.hyperstack(2):pars.hyperstack(3));
%             end  
            im1 = max(im1,[],3); 
            z = 1;
        end   
        % convert to 8-bit 512x512
       if z==1 % save 2D data 
           imOut = makeuint(imresize(im1,pars.imSize),8);
           if pars.figShowLoadedImages
                figure(pars.figShowLoadedImages); clf;
                imagesc(imOut); colormap(gray); pause(.01); 
           end
          imwrite(imOut,[saveFolder,imName,'.png']);
       else % save 3D data
            filename =[saveFolder,imName,'.tif'];
            if exist(filename,'file') == 2
                delete(filename); % needed to allow append
            end
            im1 = makeuint(im1,8); % apply to whole 3D image, since this also normalizes contrast
            for k=1:z
                frameOut =  imresize(im1(:,:,k),pars.imSize);
                imwrite(frameOut,filename,'WriteMode','append');
            end
           
       end
    end
end
%%  Run Cellpose 
if ~skipFitting
    if pars.verbose
        disp('running cellpose to ID nuclei, please wait...');
    end
    % this is 2D, since we're running on 2D images based on above 
    %  a fixed diameter is substantially faster, and may improve uniformity. 
    diameter = num2str(pars.diameter);
    if pars.useGPU
        gpu = ' --use_gpu';
    else
        gpu = '';
    end
    if pars.do_3D
       do_3D = ' --do_3D';
    else
        do_3D = '';
    end
    if strcmp(pars.device,'default')
        callPython = ['python -m cellpose --dir ',saveFolder,...
        ' --pretrained_model ',pars.model,' --diameter ',diameter,...
        ' --save_tif --no_npy',gpu,do_3D]; % ,' --resample'
    else
    callPython = ['python -m cellpose --dir ',saveFolder,...
        ' --pretrained_model ',pars.model,' --diameter ',diameter,...
        ' --device ',pars.device,...
        ' --save_tif --no_npy',gpu,do_3D]; % ,' --resample'
    end
    % pars.env = 'mlab_cellpose';
    % setEnv = ['activate ', pars.env,' ']; % '! activate cellpose ';
    % cmdOut = ['! ',pars.condaPrompt, ' ', setEnv,' && ',callPython];
    % cmdOut = ['! ',pars.activateConda,' ',pars.env,' && ',callPython];
    cmdOut = ['! ',cellpose_env,' && ',callPython];
    if pars.verbose
        disp(cmdOut);
    end
    eval(cmdOut); % should move to system Run 
end

%% load results
% ensure that the masks returned are the same size as those submitted
allMasks = FindFiles([saveFolder,'*_cp_masks.tif']);
nFOV = length(allMasks);
cellMasks = cell(nFOV,1);
if ~pars.do_3D
    for f =1:nFOV
         temp = imread(allMasks{f});
         cellMasks{f} = imresize(temp,origSize(f,1:2),'nearest');
    end
else
    for f =1:nFOV
        seg3D = zeros(origSize(f,:));
        for k=1:z
            temp = imread(allMasks{f},k);
            temp = imresize(temp,origSize(f,1:2),'nearest');
            seg3D(:,:,k) = temp;
        end
        cellMasks{f} = seg3D;
    end
end
    
%% clean up any temporary folders 
if deleteWhenDone
    delete([saveFolder,'*']); 
    rmdir(saveFolder); 
end 

if singleFile
    cellMasks = cellMasks{1};
end
