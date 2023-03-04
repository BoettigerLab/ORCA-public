function [fidSpts,datSpts,pars] = ChrTracer3p3_CropSpots(rawDataNames,regData,eTable,spots,varargin)
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
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'numParallel', 'integer',1};
% advanced / exploratory
defaults(end+1,:) = {'downSampleZ','integer',1}; % 
defaults(end+1,:) = {'checkOffset','boolean',false}; % 
defaults(end+1,:) = {'troubleshoot','boolean',false}; % 
pars = ParseVariableArguments(varargin, defaults, mfilename);




%%

% optional parameters
boxWidth = pars.boxWidth; % 16;
s = pars.currentSpot; %  22; % select spot;

numHybes = size(rawDataNames,1);

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
% here we try to acclerate things by reading only select bits rather than whole 3D dax file into RAM
% in practice it doesn't look like this accelerates things.
h=1;
infoFile = ReadInfoFile(rawDataNames{h}, 'verbose', pars.veryverbose);
% numRows = infoFile.frame_dimensions(1); % check, these dimensions might be swapped! 
% numCols = infoFile.frame_dimensions(2);


framesInDax = infoFile.number_of_frames;
frameDim = [infoFile.frame_dimensions(1)/infoFile.binning(1),...
            infoFile.frame_dimensions(2)/infoFile.binning(2)];
frameDim = uint32(frameDim);

numRows = frameDim(2); % swapped, 12/2/20 
numCols = frameDim(1);
framesInDax = uint32(framesInDax);
if ~isempty( strfind(infoFile.data_type,'little endian') )
    binaryFormat = 'l';
else
    binaryFormat = 'b';
end


%%


% =============  check zoffset
if pars.checkOffset    % check and correct off set to hybe 1
    zshift = ComputeHybeZoffset(rawDataNames,'parameters',pars,'troubleshoot',false);
end
%=====================


% now we actually load the data;
fidSpts = cell(numHybes,1);
datSpts = cell(numHybes,numDataChns);
fidSptsRaw = fidSpts;
datSptsRaw = datSpts;
for h=1:numHybes   % consider parfor this

    % we shift the center by the integer pixel amount to decide what region
    % of the dax to load. we shift that resulting image by the remaining
    % subpixel amount for the final image;
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
    [ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);  % kept [R,C] = mesh(X,Y), though FitPsf3D is [Y,X]=meshgrid(ys,xs)
    fidInds = sub2indFast([frameDim(1),frameDim(2),framesInDax],ri(:),ci(:),zi(:)); % also [H,W]. swapped dim1 and 2, despite annotation in sub2indFast
    
    % get data indicies
    datInds = cell(1,numDataChns);
    datFramesToLoad = cell(1,numDataChns);
    for n=1:numDataChns  
        isCurrDat = StringFind(frameChannels,dataChns{n},'boolean',true);
        isCurrDat(skipFrames) = false;
        selectFrames = uint32(find(isCurrDat)); 
        datFramesToLoad{n} = length(selectFrames);
        [ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);  % meshgrid should be [Y,X] = meshgrid(ys,xs)
        datInds{n} = sub2indFast([frameDim(1),frameDim(2),framesInDax],ri(:),ci(:),zi(:));
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
    % dax = ReadDax(rawDataNames{h},'verbose',false);    
    fid = fopen(rawDataNames{h});
    dataSize = frameDim(1)*frameDim(2)*framesInDax;
    dax = fread(fid, dataSize, '*uint16', binaryFormat);
    fclose(fid);
    tt = toc;
    tic
    im = dax(fidInds); 
    % im =memMaps{h}.Data(fidInds); %  movie1(fidInds); 
    ttt = toc;   
    if pars.veryverbose
        disp(['loaded in ', num2str(tt),' s']);
        disp(['indexed in ', num2str(ttt),' s']);
    end
    
    
    if strcmp(binaryFormat,'b')
        im = swapbytes(im);
    end
    im = reshape(im,[ys,xs,framesToLoadFid]); % figure(7); clf; imagesc(max(im,[],3)); colormap gray; 
    if pars.downSampleZ ~= 1
       im = im(:,:,1:pars.downSampleZ:end);
    end
    fidSptsRaw{h} = im;
    movieReg = ApplyReg(im,subPixel); % Apply hyb-specific registration
    fidSpts{h} = movieReg;
    
    %===== Get Data Images and apply drift correction ===
    for n=1:numDataChns
        % im = memMaps{h}.Data(datInds{n});
        im = dax(datInds{n});
        if strcmp(binaryFormat,'b')
            im = swapbytes(im);
        end
        im = reshape(im,[ys,xs,datFramesToLoad{n}]);
        if pars.downSampleZ ~= 1
           im = im(:,:,1:pars.downSampleZ:end);
        end
        datSptsRaw{h,n} = im;
        movieReg = ApplyReg(im,subPixel);  % apply hybe specific drift correction
        datSpts{h,n} = movieReg;
    end
end

%%
% fidSpts, datSpts
% fidSptsO = fidSpts; datSptsO = datSpts
% fidSpts = fidSptsO; datSpts = datSptsO;
%--------------
if pars.checkOffset
   for h=1:numHybes
       imF = fidSpts{h};
       bkdFid = quantile(imF(:),.001);
       [h_i,w_i,z_i] = size(imF);
       if zshift(h,1) > 0 % shift hybe down towards data
           blnkF = bkdFid*ones(h_i,w_i,zshift(h,1),class(bkdFid));
           z_e = z_i-zshift(h,1); % +1?
           newImF = cat(3,blnkF,imF(:,:,1:z_e));
           newImD = cell(1,numDataChns);
           for n=1:numDataChns
                imD = datSpts{h,n};
                bkdDat = quantile(imD(:),.001);
                blnkD = bkdDat*ones(h_i,w_i,zshift(h,1),class(bkdFid));
                newImD{1,n} = cat(3,blnkD,imD(:,:,1:z_e));
           end
       else  % shift Hybe up
           blnkF = bkdFid*ones(h_i,w_i,-zshift(h,1),class(bkdFid)); 
           z_s = 1-zshift(h,1); % zshift is neg
           newImF = cat(3,imF(:,:,z_s:end),blnkF);
           newImD = cell(1,numDataChns);
           for n=1:numDataChns
                imD = datSpts{h,n};
                bkdDat = quantile(imD(:),.001);
                blnkD =  bkdDat*ones(h_i,w_i,-zshift(h,1),class(bkdFid));
                newImD{1,n} = cat(3,imD(:,:,z_s:end),blnkD);
           end 
       end
       fidSpts(h) = {newImF};
       datSpts(h,:) = newImD;
   end
end
%---------



%%

if pars.showPlots
    try
    datSptsT = datSpts';  % This and next line, sort in alternating order.
    datMat = cat(4,datSptsT{:}); % sort in alternating order.
    fidMat = cat(4,fidSpts{:});

    % Projections
    figure(2);  % Fiducial
    stk = cellfun(@(x) max(x,[],3),fidSpts,'UniformOutput',false);
    stk1xy = cat(3,stk{:}); 
    stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSpts,'UniformOutput',false);
    stk1xz = cat(3,stk{:}); 
    nZ = size(stk1xz,1); stk1xz = stk1xz(10:nZ-10,:,:);
    subplot(2,2,1);  Ncolor(2/numHybes*IncreaseContrast(stk1xy)); title('corr. fid x,y global sub-pixel aligned');
    subplot(2,2,2); Ncolor(2/numHybes*IncreaseContrast(stk1xz)); title('corr. fid x,z  global sub-pixel aligned');

    stk = cellfun(@(x) max(x,[],3),fidSptsRaw,'UniformOutput',false);
    stk1xy = cat(3,stk{:}); 
    stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),fidSptsRaw,'UniformOutput',false);
    stk1xz = cat(3,stk{:}); 
    nZ = size(stk1xz,1); stk1xz = stk1xz(10:nZ-10,:,:);
    subplot(2,2,3);  Ncolor(2/numHybes*IncreaseContrast(stk1xy)); title('corr. fid x,y global pixel aligned');
    subplot(2,2,4); Ncolor(2/numHybes*IncreaseContrast(stk1xz)); title('corr. fid x,z  global pixel aligned');

    % Data Projections
    figure(3);  clf;
    stk = cellfun(@(x) max(x,[],3),datSptsT(:),'UniformOutput',false);
    dat1xy = cat(3,stk{:});          
    stk = cellfun(@(x) max(permute(x,[3,2,1]),[],3),datSptsT(:),'UniformOutput',false);
    dat1xz = cat(3,stk{:});
    dat1xz(dat1xz==0) = mode(nonzeros(dat1xz(:)));
    nHybs = size(dat1xy,3);
    subplot(1,2,1); im = Ncolor(2/nHybs*dat1xy); % someday we will write saturation controls into Ncolor 
    imagesc(IncreaseContrast(im)); title('corr. data x,y');
    subplot(1,2,2); im = Ncolor(2/nHybs*dat1xz); 
    imagesc(IncreaseContrast(im)); title('corr. data x,z');
    
          
% get tile labels
    if ~isempty(pars.eTable)
         [tileLabels_fid,tileLabels_dat] =TileLabelsFromEtable(pars.eTable);
    else
        tileLabels_dat = cellstr(num2str( (1:numHybes*nChns)')) ;
        tileLabels_fid = cellstr(num2str( (1:numHybes)')) ;
    end

    % 4D tiles
    figure(44); clf;
    PlotProjection4D(fidMat,'fits',[],'projection','xy','tileLabels',tileLabels_fid); title('fid xy');
    figure(5); clf;
    PlotProjection4D(fidMat,'fits',[],'projection','xz','tileLabels',tileLabels_fid); title('fid xz');
    figure(6); clf;
    PlotProjection4D(datMat,'fits',[],'projection','xy','tileLabels',tileLabels_dat); title('data xy');
    figure(7); clf;
    PlotProjection4D(datMat,'fits',[],'projection','xz','tileLabels',tileLabels_dat); title('data xz');
    
    % max dataspots

    figure(10); clf;
    for d=1:numDataChns
        maxDat = cellfun(@(x) quantile(x(:),.999),datSpts(:,d));
        minDat = cellfun(@(x) quantile(x(:),.02),datSpts(:,d));
        subplot(numDataChns,1,d);
        bar(maxDat); hold on; bar(minDat);
        title(['channel ',num2str(d)]);
        legend('data spot','data background'); 
    end
    catch er
        warning(er.getReport);
        cprintf([1 .5 0],'error trying to plot figures');
        disp('place debug here');
    end
end
