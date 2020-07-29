function [imageTiles,stageXY,pars] = DaxToImageTiles(daxFiles,varargin)
% [imageTiles,stageXY,pars] = DaxToImageTiles(daxFiles,varargin)
% loads daxfiles
% 
% inputs: daxFiles - a cell array of dax files to load
% outputs: imageTiles - a cell array of image tiles 
%
%
%% optional parameters
% defaults(end+1,:) = {'pix_to_mm','positive',6.55};  % 6.55 = scope 1. 6.45 = scope 2
% defaults(end+1,:) = {'trimBorder','nonnegative',0};
% defaults(end+1,:) = {'selectFrame','integer',2};
% defaults(end+1,:) = {'projectAll','boolean',true};
% defaults(end+1,:) = {'transpose','boolean',true}; % see the mosaic parameters in the Hal parameter file used to collect the movie 
% defaults(end+1,:) = {'fliplr','boolean',true};  % see the mosaic parameters in the Hal parameter file used to collect the movie
% defaults(end+1,:) = {'flipud','boolean',false};  % see the mosaic parameters in the Hal parameter file used to collect the movie
% defaults(end+1,:) = {'verbose','boolean',false};
% defaults(end+1,:) = {'saveProjections','boolean',true};  % save the max projected images for quick loading 
% defaults(end+1,:) = {'loadProjections','boolean',true};  % load previousl save the max projected images for quick loading 
% % parameters useful for iterative mosaic overlays
% defaults(end+1,:) = {'offset','float',0}; % this is useful to pass from previous mosaics to have multiple files in same register
% defaults(end+1,:) = {'buffer','nonnegative',1000};

%% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'positionsFile','string',''};  %
defaults(end+1,:) = {'trimBorder','nonnegative',0};
defaults(end+1,:) = {'selectFrame','integer',2}; % show only a single frame (overriden by projectAll)
defaults(end+1,:) = {'projectAll','boolean',true}; % show a projection of x,y data
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'saveProjections','boolean',true};  % save the max projected images for quick loading 
defaults(end+1,:) = {'loadProjections','boolean',true};  % load previously save the max projected images for quick loading 
defaults(end+1,:) = {'loadInfoOnly','boolean',false};  % load only the info files, return only stageXY 
% parameters useful for iterative mosaic overlays
defaults(end+1,:) = {'offset','float',0}; % this is useful to pass from previous mosaics to have multiple files in same register
defaults(end+1,:) = {'buffer','nonnegative',1000};
% scope specific Pars
defaults(end+1,:) = {'scope',{'autoDetect','scope1','scope2','scope3','other'},'autoDetect'};  %
defaults(end+1,:) = {'pix_to_mm','positive',6.55};  % 6.55 = scope 1. 6.45 = scope 2
defaults(end+1,:) = {'transpose','boolean',true}; % see the mosaic parameters in the Hal parameter file used to collect the movie 
defaults(end+1,:) = {'fliplr','boolean',true};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'flipud','boolean',false};  % see the mosaic parameters in the Hal parameter file used to collect the movie
pars = ParseVariableArguments(varargin,defaults,mfilename);

% load data
positionsTable = []; % default to empty
numIms = length(daxFiles);
imageTiles = cell(numIms,1);
stagePos = zeros(numIms,2);
if pars.verbose
    disp('loading data...');
end
for f = 1:numIms
    if pars.projectAll 
        [folder,filename] = fileparts(daxFiles{f});
        maxFile = [folder,filesep,'max_',filename,'.dax'];
        if contains(filename,'Max')
            maxFile = [folder,filesep,filename,'.dax'];
        end
        if exist(maxFile,'file')~=0 && pars.loadProjections  % use previously written
            [dax,infoFile] = ReadDax(maxFile,'verbose',pars.veryverbose);
        elseif pars.loadInfoOnly
            infoFile = ReadInfoFile(maxFile,'verbose',pars.veryverbose);
            dax = [];
        else
            [dax,infoFile] = ReadDax(daxFiles{f},'verbose',pars.veryverbose);
            dax = max(dax,[],3); 
            if pars.saveProjections % save max projections
                infoOut = infoFile;
                infoOut.number_of_frames= 1;
                infoOut.localName = ['max_',infoOut.localName];
                WriteDAXFiles(dax,infoOut,'verbose',pars.verbose);
            end
        end
    else
        [dax,infoFile] = ReadDax(daxFiles{f},'startFrame',pars.selectFrame,'endFrame',pars.selectFrame,'verbose',pars.veryverbose);
    end   
    
    % clean up dax file
    if f==1 % all data will be from same scope
        pars = GetScopeSettings(infoFile,'parameters',pars);
    end
    
    if ~isempty(dax)
        z_i = size(dax,3); % assumes all dax are the same size;
        if pars.transpose
            if z_i ==1
                dax = dax';
            else
                dax = permute(dax,[2,1,3]);
            end
        end
        if pars.fliplr
            dax = fliplr(dax);
        end
        if pars.flipud
            dax = flipud(dax);
        end
        if pars.trimBorder > 0
            dax = dax(pars.trimBorder+1:end-pars.trimBorder,pars.trimBorder+1:end-pars.trimBorder,:);
        end
        imageTiles{f} = dax;
    end
    noStagePos = infoFile.Stage_X == 0 && infoFile.Stage_Y ==0;
    if noStagePos  && isempty(positionsTable) && isempty(pars.positionsFile)
        [positionsFile,positionsFolder] = uigetfile('*position*.txt');
        pars.positionsFile = [positionsFolder,filesep,positionsFile];
        positionsTable = readtable(pars.positionsFile);
    elseif noStagePos  && isempty(positionsTable) && ~isempty(pars.positionsFile)
        positionsTable = readtable(pars.positionsFile);
    end
    if noStagePos
        infoFile.Stage_X = positionsTable{f,1};
        infoFile.Stage_Y = positionsTable{f,2}; 
    end
    stagePos(f,:) = [infoFile.Stage_X,infoFile.Stage_Y]*pars.pix_to_mm;    
end
% [h_i,w_i,~] = size(dax); 
% assumes all dax are the same size;
h_i = infoFile.frame_dimensions(1);
w_i = infoFile.frame_dimensions(2);
    

% use optional offset
if pars.offset == 0
    %  This is helpful to align multiple mosaics to the same stage coordinates
    pars.offset = - min(stagePos) + [h_i/2,w_i/2] + pars.buffer*[1,1];
end
% update stageXY
stageXY = stagePos + pars.offset;