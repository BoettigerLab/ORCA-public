function [fidSpts,datSpts,pars] = ChrTracer3p2_CropSpots(rawDataNames,regData,eTable,spots,varargin)
% 
% -------------------------------------------------------------------------
% Required Inputs
% -------------------------------------------------------------------------
% 'parameters',pars
% A structure containing pars.fiducialFrames, and pars.dataFrames
% Both are 4D matrices, nRows x nCols x nStacks x nHybes 
% 
% -------------------------------------------------------------------------
% Optional inputs
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
% Notes
% -------------------------------------------------------------------------
% called by ChrTracer

% supress some unnecessary warnings. 
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
warning('off','MATLAB:prnRenderer:opengl');

defaults = cell(0,3);
% current spot 
defaults(end+1,:) = {'currentSpot','integer',1};
defaults(end+1,:) = {'boxWidth', 'positive', 16};
% FOV parameters
defaults(end+1,:) = {'fov', 'integer', 0};  % Field of view in experiment,  0 for not known;
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'numParallel', 'integer',1};
pars = ParseVariableArguments(varargin, defaults, mfilename);




%%

% optional parameters
boxWidth = pars.boxWidth; % 16;
s = pars.currentSpot; %  22; % select spot;

numHybes = size(rawDataNames,1);
memMaps = cell(numHybes,1);
for h=1:numHybes
    memMaps{h} = memmapfile(rawDataNames{h}, ...
                'Format', 'uint16', ...
                'Writable', false, ...
                'Offset', 0); %  (startFrame-1)*frameSize*16/8);  
end


% ------ Identify fiducial channel and parameter channel ------
% Determine number of channels based on first hybe 
[isFidChannel,frameChannels,~,~,dataChns] = GetFidChnFromTable(eTable);
nBufFrames = eTable.bufferFrames(1);
totFrames  = eTable.totalFrames(1);
skipFrames = false(1,totFrames);
skipFrames([1:nBufFrames,totFrames-nBufFrames+1:totFrames]) = true;
keepFrames = nBufFrames+1:totFrames-nBufFrames;
numDataChns = length(dataChns);



%% parse info file 
infoFile = ReadInfoFile(rawDataNames{h}, 'verbose', false);
numRows = infoFile.frame_dimensions(1); % check, these dimensions might be swapped! 
numCols = infoFile.frame_dimensions(2);


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


%%


% now we actually load the data;
fidSpts = cell(numHybes,1);
datSpts = cell(numHybes,numDataChns);
fidSptsRaw = fidSpts;
datSptsRaw = datSpts;
for h=1:numHybes   % consider parfor this

    driftX = floor(regData(h).xshift);
    driftY = floor(regData(h).yshift);
    subPixel = struct();
    subPixel.xshift = regData(h).xshift-driftX;
    subPixel.yshift = regData(h).yshift-driftY;
    
    yi = round(max(1,spots(s,2)-boxWidth/2-driftY));
    ye = round(min(spots(s,2)+boxWidth/2-driftY,numRows));
    xi = round(max(1,spots(s,1)-boxWidth/2-driftX));
    xe = round(min(spots(s,1)+boxWidth/2-driftX,numCols));
    
  %   disp([xi,xe,yi,ye]);  % for troubleshooting

    % get fiducial indices
    isFidChannel(skipFrames) = false; % could move out of loop
    selectFrames = uint32(find(isFidChannel)); 
    framesToLoadFid = length(selectFrames);
    [ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);
    fidInds = sub2indFast([frameDim(2),frameDim(1),framesInDax],ri(:),ci(:),zi(:));

    % get data indicies
    datInds = cell(1,numDataChns);
    datFramesToLoad = cell(1,numDataChns);
    for n=1:numDataChns  
        isCurrDat = StringFind(frameChannels,dataChns{n},'boolean',true);
        isCurrDat(skipFrames) = false;
        selectFrames = uint32(find(isCurrDat)); 
        datFramesToLoad{n} = length(selectFrames);
        [ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);
        datInds{n} = sub2indFast([frameDim(2),frameDim(1),framesInDax],ri(:),ci(:),zi(:));
    end


    % these are the same for all hybes;
    xs = xe-xi+uint32(1);
    ys = ye-yi+uint32(1);

    %====== Get Fid Images and apply drift correction ========
    % load the ROI from disk using the memory map
    % 
    % Note for future Dev with rotation:
    % if rotation is allowed, fidInds need to be recomputed based on
    % rotation, so that the correct area is loaded. The resulting image
    % also will need to be rotated.  This is probably best accomplished by
    % using translation rotation matrices to address the rotate-about-point
    % constraint.
    
    % movie = memMaps{h}.Data(fidInds);
    % movie1  = memMaps{h}.Data;
    tic; % this is indexing step is stochastically slow. 
    im =memMaps{h}.Data(fidInds); %  movie1(fidInds); 
    tt = toc;   
    if pars.verbose
        disp(['indexed in ', num2str(tt),' s']);
    end

    
    
    if strcmp(binaryFormat,'b')
        im = swapbytes(im);
    end
    im = reshape(im,[ys,xs,framesToLoadFid]); % figure(7); clf; imagesc(max(movie,[],3)); colormap gray;   
    fidSptsRaw{h} = im;
    movieReg = ApplyReg(im,subPixel); % Apply hyb-specific registration
    fidSpts{h} = movieReg;
    
    %===== Get Data Images and apply drift correction ===
    for n=1:numDataChns
        im = memMaps{h}.Data(datInds{n});
        if strcmp(binaryFormat,'b')
            im = swapbytes(im);
        end
        im = reshape(im,[ys,xs,datFramesToLoad{n}]);
        datSptsRaw{h,n} = im;
        movieReg = ApplyReg(im,subPixel);  % apply hybe specific drift correction
        datSpts{h,n} = movieReg;
    end
end


%%

if pars.showPlots
    datSptsT = datSpts';
    datMat = cat(4,datSptsT{:});
    fidMat = cat(4,fidSpts{:});

    % Projections
    figure(2);  % Fiducial
    stk = cellfun(@(x) max(x,[],3),fidSpts,'UniformOutput',false);
    stk1xy = cat(3,stk{:}); 
    stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSpts,'UniformOutput',false);
    stk1xz = cat(3,stk{:}); 
    nZ = size(stk1xz,1); stk1xz = stk1xz(10:nZ-10,:,:);
    subplot(2,2,1);  Ncolor(2/numHybes*IncreaseContrast(stk1xy)); title('corr. fid x,y aligned');
    subplot(2,2,2); Ncolor(2/numHybes*IncreaseContrast(stk1xz)); title('corr. fid x,z  aligned');

    stk = cellfun(@(x) max(x,[],3),fidSptsRaw,'UniformOutput',false);
    stk1xy = cat(3,stk{:}); 
    stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSptsRaw,'UniformOutput',false);
    stk1xz = cat(3,stk{:}); 
    nZ = size(stk1xz,1); stk1xz = stk1xz(10:nZ-10,:,:);
    subplot(2,2,3);  Ncolor(2/numHybes*IncreaseContrast(stk1xy)); title('corr. fid x,y raw');
    subplot(2,2,4); Ncolor(2/numHybes*IncreaseContrast(stk1xz)); title('corr. fid x,z  raw');

    % Data Projections
    figure(3);  clf;
    stk = cellfun(@(x) max(x,[],3),datSptsT(:),'UniformOutput',false);
    dat1xy = cat(3,stk{:});          
    stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),datSptsT(:),'UniformOutput',false);
    dat1xz = cat(3,stk{:});
    dat1xz(dat1xz==0) = mode(nonzeros(dat1xz(:)));
    subplot(1,2,1); im = Ncolor(dat1xy); imagesc(IncreaseContrast(im)); title('corr. data x,y');
    subplot(1,2,2); im = Ncolor(dat1xz); imagesc(IncreaseContrast(im)); title('corr. data x,z');
    
        
    figure(4); clf;
    PlotProjection4D(fidMat,'fits',[],'projection','xy'); title('fid xy');
    figure(5); clf;
    PlotProjection4D(fidMat,'fits',[],'projection','xz'); title('fid xz');
    figure(6); clf;
    PlotProjection4D(datMat,'fits',[],'projection','xy'); title('data xy');
    figure(7); clf;
    PlotProjection4D(datMat,'fits',[],'projection','xz'); title('data xz');
end
