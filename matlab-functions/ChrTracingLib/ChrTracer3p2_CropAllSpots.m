function [fidiSpots,dataSpots,pars] = ChrTracer3p2_CropAllSpots(rawDataNames,regData,eTable,spots,varargin)
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
defaults(end+1,:) = {'numParallel', 'integer',12};
pars = ParseVariableArguments(varargin, defaults, mfilename);




%%
tic;

% optional parameters
boxWidth = pars.boxWidth; % 16;

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
isFidChannel(skipFrames) = false; % could move out of loop
numDataChns = length(dataChns);



%% parse info file 
infoFile = ReadInfoFile(rawDataNames{1}, 'verbose', false);
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
numSpots = size(spots,1);
fidiSpots = cell(numSpots,1);
dataSpots = cell(numSpots,1);

if pars.numParallel > 1   
    p = gcp('nocreate');
    if pars.numParallel > 1 && isempty(p)
        disp('creating parallel pool ');
        tic
        parpool(pars.numParallel);
        p = gcp('nocreate');
        t = toc;
        disp([num2str(p.NumWorkers),' core pool created in ',num2str(t,3),' s']);     
    end

    %% main parallel processing

    parfor s=1:numSpots

        % now we actually load the data;
        fidSpts = cell(numHybes,1);
        datSpts = cell(numHybes,numDataChns);
        fidSptsRaw = fidSpts;
        datSptsRaw = datSpts;
        for h=1:numHybes

            driftX = floor(regData(h).xshift); %#ok<PFBNS>  (% this is a tiny structure)
            driftY = floor(regData(h).yshift);
            subPixel = struct();
            subPixel.xshift = regData(h).xshift-driftX;
            subPixel.yshift = regData(h).yshift-driftY;

            yi = round(max(1,spots(s,2)-boxWidth/2-driftY));
            ye = round(min(spots(s,2)+boxWidth/2-driftY,numRows));
            xi = round(max(1,spots(s,1)-boxWidth/2-driftX)); %#ok<PFBNS>
            xe = round(min(spots(s,1)+boxWidth/2-driftX,numCols));

            % compute pads
            % pad image so all crops are the same size
            needPad = false;
            xshift = 0;
            yshift = 0;
             %positive adds to the front, negative adds to the back
            if yi == 1
                yshift = 1-(spots(s,2)-boxWidth/2-driftY);
                needPad = true;
            end
            if ye == numRows
                yshift = -(spots(s,2)+boxWidth/2-driftY-numRows);
                needPad = true;
            end
            if xi == 1
                xshift = 1-(spots(s,1)-boxWidth/2-driftX);
                needPad = true;
            end
            if xe == numCols
                xshift = -(spots(s,1)+boxWidth/2-driftX-numCols);
                needPad = true;
            end
            

            
            % disp([xi,xe,yi,ye]);  % for troubleshooting

            % get fiducial indices
            
            selectFrames = uint32(find(isFidChannel)); 
            framesToLoadFid = length(selectFrames);
            [ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);
            fidInds = sub2indFast([frameDim(2),frameDim(1),framesInDax],ri(:),ci(:),zi(:)); %#ok<PFBNS>

            % get data indicies
            datInds = cell(1,numDataChns);
            datFramesToLoad = cell(1,numDataChns);
            for n=1:numDataChns  
                isCurrDat = StringFind(frameChannels,dataChns{n},'boolean',true); %#ok<PFBNS>
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
            movie = memMaps{h}.Data(fidInds); %#ok<PFBNS> % Necessary overhead
            if strcmp(binaryFormat,'b')
                movie = swapbytes(movie);
            end
            movie = reshape(movie,[ys,xs,framesToLoadFid]); % figure(7); clf; imagesc(max(movie,[],3)); colormap gray;   
            fidSptsRaw{h} = movie;
            movieReg = ApplyReg(movie,subPixel); % Apply hyb-specific registration
            % make sure all cropped movies are same dimension 
            if needPad
                movieReg = TranslateImage(movieReg,xshift,yshift,'trim',false) ;  % positive pads left/pre, negative pads right/post 
            end
            fidSpts{h} = movieReg;

            %===== Get Data Images and apply drift correction ===
            for n=1:numDataChns
                movie = memMaps{h}.Data(datInds{n});
                if strcmp(binaryFormat,'b')
                    movie = swapbytes(movie);
                end
                movie = reshape(movie,[ys,xs,datFramesToLoad{n}]);
                datSptsRaw{h,n} = movie;
                movieReg = ApplyReg(movie,subPixel);  % apply hybe specific drift correction
                % make sure all cropped movies are same dimension 
                if needPad
                    movieReg = TranslateImage(movieReg,xshift,yshift,'trim',false) ; % positive pads left/pre, negative pads right/post 
                end
                datSpts{h,n} = movieReg;
            end
        end
        
        fidiSpots{s} = fidSpts;
        dataSpots{s} = datSpts;
    end
end

t=toc;
disp(['Finished Crop All in ',num2str(t/60),' min.']);
