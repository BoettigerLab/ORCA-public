function output = ProcessTracesGUI(gui_step,pars,data,varargin)
% 
%% see ProcessTracesApp for a mini-Step GUI that wraps this. 
%
%
% 
%% Output
% pntArray Matrix nTraces x nFrames x 6 data-values:
%    data-types = x_pix, y_pix, z_nm, t_frame, brightness, background 
%    Outputs two point arrays     
% 
%% Updates
%    rewriting DaoFitZ to improve z-fits. This code was developed on the
%    100 nm step data and needs some parameter optimization and testing
%    still for the 1 um step data.
%        We should probably modify it to allow a pre-computed Z-PSF to be
%        given, rather than learning the PSF from the data in each spot.
% 
% 
%    LinkLiveStruct - updated based on z-projection linkage to change the
%    fraction of frames used (was hard-coded to use only the first half). 
%
%
% 

% Defaults
npp = 108; % nm per pixel


%% Main function

%% stp 1, load data

if gui_step == 1


    if isempty(pars)
        disp('leave saveFolder blank to use dataFolder');
        disp('leave daxName1 blank to use open a file browser');

        defaults = cell(0,3);
        % data loading
        defaults(end+1,:) = {'saveFolder','string',''}; 
        defaults(end+1,:) = {'daxFile1','string',''}; 
        defaults(end+1,:) = {'alignment_file','string',''}; % 
        defaults(end+1,:) = {'chrom_correct_file','string',''}; % 
        defaults(end+1,:) = {'bin_tag','string','_2d_iters'}; % 
        defaults(end+1,:) = {'framesToLoad','integer',0}; % 0 for load all frames
        defaults(end+1,:) = {'zDepth','integer',0}; 
        defaults(end+1,:) = {'verbose','boolean',true}; 
        pars = ParseVariableArguments(varargin,defaults,mfilename);
        pars.directions = 'Load data';
        output = pars;
    else
        disp('running...')
        % putting the uigetfiles in the Execute is better
        if isempty(pars.daxFile1)
            [daxName,dataFolder] = uigetfile('*.dax','select a dax movie to load');
            daxFile1 = [dataFolder,daxName];
            
        else
            daxFile1 = pars.daxFile1;
            [dataFolder,~] = fileparts(daxFile1);
        end
        % === Get source folder and calibraion files
        if isempty(pars.alignment_file)
            folderCal2 = 'M:\2023-08-26_130kb_rad21_dTag\';  % update me
            pars.alignment_file = [folderCal2,'alignmentData.txt'];
        else
            if pars.verbose
                disp(['using camera alignment file: ',pars.alignment_file]);
            end
        end
        if isempty(pars.chrom_correct_file)
            folderCal2 = 'M:\2023-08-26_130kb_rad21_dTag\';
            pars.chrom_correct_file = [folderCal2,'tform3D.mat'];
        else
            if pars.verbose
                disp(['using camera alignment file: ',pars.chrom_correct_file]);
            end
        end
        alignT = readtable(pars.alignment_file);
        alignS = table2struct(alignT);
        
        % error checking -- make sure files exist
        if ~exist(pars.alignment_file,'file')
            error(['unable to find ',pars.alignment_file, ', check filepath']);
        end
        if ~exist(pars.chrom_correct_file,'file')
            error(['unable to find ',pars.chrom_correct_file, ', check filepath']);
        end

        % create a save folder if none was passed
        if isempty(pars.saveFolder)
            pars.saveFolder = [dataFolder,'Analysis\'];
            if ~exist(pars.saveFolder,'dir')
                mkdir(pars.saveFolder)
            end
        end

        if pars.zDepth == 0
            offFile =  regexprep(daxFile1,'_C1.dax','.off');
            offTable = ReadTableFile(offFile,'delimiter',' ');
            stageZ = offTable.offset;
            % auto determine cycle
                currZ = 0; oldZ = -inf; 
                [~,z0] =min(stageZ(1:100));
                z = z0;
                while currZ > oldZ
                    oldZ = stageZ(z);
                    z=z+1;
                    currZ = stageZ(z);
                end
                pars.zDepth = z-z0;
           disp(['system determines the stack height = ',num2str(pars.zDepth)]);
           figure(10); clf; 
           subplot(2,1,1); plot(stageZ(1:100),'.-');
            title(['stackHeight = ',num2str(pars.zDepth)'])
            subplot(2,1,2); plot(stageZ,'.-');
        end
        
        % daxName1 = [dataFolder,filesep,daxName,'.dax'];
        daxFile2 = regexprep(daxFile1,'C1','C2'); %  [dataFolder,'36mW_0001_C2.dax'];
        binName1 = regexprep(daxFile1,'.dax',[pars.bin_tag,'.hdf5']); % [dataFolder,'36mW_0001_C1_2d_iters.hdf5'];  % 2d no iters
        binName2 = regexprep(daxFile2,'.dax',[pars.bin_tag,'.hdf5']);
        


        % === load the images 
        f = 1; 
        [im1f,info1] = ReadDax(daxFile1,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);  % note the +3 since its a 4 step series 
        [im2f,info2] = ReadDax(daxFile2,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);
        im1 = max(im1f,[],3);
        im2 = max(im2f,[],3);
        im2 = ApplyReg(fliplr(im2),alignS);

        overlayImage = cat(3,im1,im2); % higher contrast 
        overlayImage = IncreaseContrast(overlayImage,'high',.9999,'low',.1);
        figure(1); clf; Ncolor(overlayImage);  axis image;
        [h,w] = size(im1);
        

        % === Load the fit data
        if pars.framesToLoad <= 0 
            nFrames = info1.number_of_frames;  %
        else
            nFrames = pars.framesToLoad;
        end
        % nFrames = 200; % shorten for testing
        if pars.verbose
            disp('loading channel 1 fits:')
        end
        fits1 = LoadHD5Fits(binName1,'nFrames',nFrames,'verbose',pars.verbose);    
        if pars.verbose
            disp('loading channel 2 fits:')
        end
        fits2 = LoadHD5Fits(binName2,'nFrames',nFrames,'verbose',pars.verbose);
        
        % === Apply camera alignment and chromatic correction
        fits2c = Register2CamFits(fits2,'alignment_file',pars.alignment_file,'chrom_correct_file',pars.chrom_correct_file);

        % Load an overlay movie 
        T = 10; % number of frames
        selFrames = floor(linspace(1,nFrames,T));
        [h,w,~]=size(overlayImage);
        overlayMovie = zeros(h,w,2,T,'uint16');
        overlayMovie(:,:,:,1) = overlayImage;
        if pars.verbose
            disp('loading movie frames...')
        end
        for i=1:T-1
            if pars.verbose
                disp([num2str(100*(1+i)/T),'% complete']);
            end
            f=selFrames(1+i);
            im1 = ReadDax(daxFile1,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);  % note the +3 since its a 4 step series 
            im2 = ReadDax(daxFile2,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);
            im1 = max(im1,[],3);
            im2 = max(im2,[],3);
            im2 = ApplyReg(fliplr(im2),alignS);
            overlayMovie(:,:,1,1+i) = IncreaseContrast(im1,'high',.99999,'low',.1);
            overlayMovie(:,:,2,1+i) = IncreaseContrast(im2,'high',.99999,'low',.1);
        end
        
        for f=1:T  % play movie backwards (to end on brightest frame); 
            figure(1); clf;
             Ncolor(overlayMovie(:,:,:,T-f+1));  axis image;
             pause(.01);
        end

        output.daxFile1 = daxFile1;
        output.fits1 = fits1;
        output.fits2 = fits2;
        output.fits2c = fits2c;
        output.alignS = alignS;
        output.overlayImage = overlayImage;
        output.overlayMovie = overlayMovie;
        output.saveFolder = pars.saveFolder;
        output.dataFolder = dataFolder; 
        disp('step completed');
    end
end

%% ID pairs
% (may want to split up Match / link / Match
if gui_step == 2
    % not the step to be eliminating doublets?
    %   current linking scheme is exclusive - any spot can only belong to one trace.  
    %   maybe we should just consider doublets as separate and possibly colliding spots and we allow individual spots to be counted in multiple traces? 

    if isempty(pars)
        defaults = cell(0,3);
        % Trace identification 
        defaults(end+1,:) = {'seedThresh','float',0.996}; % autoselect threshold for grouping localizations in z. 
        defaults(end+1,:) = {'seedBinResolution','positive',6}; % autoselect threshold downsampling for grouping localizations in z. 
        defaults(end+1,:) = {'maxSep','positive',25}; % max separation for spots to be considered a pair (in pixels)
        defaults(end+1,:) = {'removeDots','positive',[]}; % max separation for spots to be considered a pair (in pixels)
        pars = ParseVariableArguments(varargin,defaults,mfilename);
        pars.directions = 'Find spot pairs';
        output = pars;
    else

        % process inputs
        fits1 = data.fits1;
        fits2c = data.fits2c;
        overlayImage = data.overlayImage;

        % === auto-select centers for dancing spots 
        xy1 = TraceCenters({fits1.x},{fits1.y},...
            'autoSelectThreshold',pars.seedThresh,...
            'binResolution',pars.seedBinResolution);  % ,'showPlots',true
        xy2 = TraceCenters({fits2c.x},{fits2c.y},...
            'autoSelectThreshold', pars.seedThresh,...
            'binResolution',pars.seedBinResolution);
        

        
        % === pair selected spots
        % keep only dancing spots that have partners in the other chromatic channel
        [matched,cost] = MatchPoints(xy1,xy2,'maxDist',pars.maxSep);
        m =matched(cost<pars.maxSep,:);      
        nMatched = size(m,1);
        xy1p = xy1(m(:,1),1:2);
        xy2p = xy2(m(:,2),1:2); % indexed in order of starting point
            
        % ---- remove pairs if requested (based on previous run of step)
        if ~isempty(pars.removeDots)
            xy1p(pars.removeDots,:) = [];
            xy2p(pars.removeDots,:) = [];
            nMatched = size(xy1p,1);
        end

        % ---- view
        figure(1); clf; Ncolor(overlayImage);  axis image;
        hold on; plot(xy1p(:,1),xy1p(:,2),'y+');
        hold on; plot(xy2p(:,1),xy2p(:,2),'gs');
        sNum = cellstr(num2str((1:nMatched)'));
        text(xy1p(:,1),xy1p(:,2),sNum,'color','r'); hold on;
        % --------

        % view all points in an image array with slider 
        nP = size(xy1p,1);
        xf = xy1p(:,1);
        yf = xy1p(:,2);
        r= 40;
        [h,w,~,T] = size(data.overlayMovie); % (this one is big, we keep it in 'data' to avoid duplication)   
        imTile1 = zeros(r,r,nP,T,'uint16');
        imTile2 = zeros(r,r,nP,T,'uint16');
        for p=1:nP
            x1 = max([1,xf(p)-r]);
            x2 = min([w,xf(p)+r]);
            y1 = max([1,yf(p)-r]);
            y2 = min([h,yf(p)+r]);
            xL = x2-x1+1;
            yL = y2-y1+1;
            if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                disp(([y1,y2,x1,x2]))
                disp(([h,w]));
                disp('debug here')
            end
            im1 = squeeze(data.overlayMovie(y1:y2,x1:x2,1,:));
            im2 = squeeze(data.overlayMovie(y1:y2,x1:x2,2,:));
            im1 = IncreaseContrast(im1,'high',.99999,'low',.1);
            im2 = IncreaseContrast(im2,'high',.99999,'low',.1);
            imTile1(1:yL,1:xL,p,:) = im1;
            imTile2(1:yL,1:xL,p,:) = im2;
        end
        % the display is the slow
        t=1;
        [~,labelOffsets,imT1] = TileImage(imTile1(:,:,:,t),'colormap',gray(256),'numRows',10);
        [~,~,imT2] = TileImage(imTile2(:,:,:,t),'colormap',gray(256),'numRows',10);            
        imO = Ncolor(cat(3,imT1,imT2));
        nums = cellstr(num2str( (1:nP)' ));

        [h,w,c] = size(imO);
        movieIn = zeros(h,w,c,T,'uint16');
        movieIn(:,:,:,1) = imO;
        for t = 2:T
            [~,~,imT1] = TileImage(imTile1(:,:,:,t),'colormap',gray(256),'numRows',10);
            [~,~,imT2] = TileImage(imTile2(:,:,:,t),'colormap',gray(256),'numRows',10);            
            movieIn(:,:,:,t) = Ncolor(cat(3,imT1,imT2));
        end
        options.labels.x = labelOffsets(:,1);
        options.labels.y = labelOffsets(:,2)+6;
        options.labels.text = nums;
        options.labels.color = 'w';
        [fig,ax] = MovieSlider(movieIn,'firstFrame',1,'options',options);

        fig.Position(1:2) = [0,0];

    % === outputs
    output.fig = fig;
    output.movieIn = movieIn;
    output.movieOptions = options;
    output.seedPointsClr1 = xy1p;
    output.seedPointsClr2 = xy2p;

    end
end

%% Link points
if gui_step == 3
    if isempty(pars)
        defaults = cell(0,3);
        % Linking
        defaults(end+1,:) = {'maxStep','positive',6}; % max step since last observed frame to consider linking (in pixels)
        defaults(end+1,:) = {'maxDistToSeed','positive',16}; % max distance from seed point to link in a trace (in pixels)
         
        % Merging traces
        defaults(end+1,:) = {'maxTraceSep','positive',30};  %       maxTraceSep = 30;     
        defaults(end+1,:) = {'maxOverlapFrac','nonnegative',.01}; %  maxOverlapFrac  = 0.01;
        defaults(end+1,:) = {'maxStepVarIncrease','positive',2};  % fold change increase in step-size variation after trace merge maxStepVarIncrease = 2; 
        defaults(end+1,:) = {'maxSep','positive',25};  % max separation for relinking pairs
        defaults(end+1,:) = {'minPointsPerTrace','positive',100};  %  minPointsPerTrace = 100;
        % smoothing
        defaults(end+1,:) = {'smooth','boolean',true}; % use smoothing?
        defaults(end+1,:) = {'maxGap','integer',10}; %  max gap to attempt linear fill on (in frames)
        defaults(end+1,:) = {'smoothWindow','integer',4}; %  moving average smoothing window size (in frames) (should be 1 larger than the number of z-frames, otherwise we have nothing to average between)
        pars = ParseVariableArguments(varargin,defaults,mfilename);
        pars.directions = 'Link spots into traces';
        output = pars;
    else
         % process inputs
        fits1 = data.fits1;
        fits2c = data.fits2c;
        overlayImage = data.overlayImage;
        xy1p = data.seedPointsClr1;
        xy2p = data.seedPointsClr2;
        nPts = size(xy1p,1);
        sNum = cellstr(num2str( (1:nPts)' ));

        % 
  

        %===== main function
        nFrames = length(fits1);
        % === Link spots into traces 
        pntArray1 = LinkLiveStruct(fits1(1:nFrames),'seedPoints',xy1p,...
            'maxStep',pars.maxStep,'maxDistToSeed',pars.maxDistToSeed,'uniqueMatch',true); 
        pntArray2 = LinkLiveStruct(fits2c(1:nFrames),'seedPoints',xy2p,...
            'maxStep',pars.maxStep,'maxDistToSeed',pars.maxDistToSeed,'uniqueMatch',true); 
         
        if ~exist(data.saveFolder,'dir')
            mkdir(data.saveFolder)
        end

        %------ Show Results
        f4 = figure(4); clf;
        f4.Name = 'All Traces';
        nFrames = size(pntArray1,2);
        frameT = (1:nFrames);
        frameT1 = repmat(frameT,nPts,1);
        Ncolor(overlayImage);  axis image; hold on;
        %--- convert fits to 
        fit1s = fits1;  
        fit2s = fits2c; 
        for s=1:length(fit1s)
            fit1s(s).frame = repmat(fit1s(s).frame,length(fit1s(s).x),1);
            fit2s(s).frame = repmat(fit2s(s).frame,length(fit2s(s).x),1);
        end
        scatter(cat(1,fit1s.x),cat(1,fit1s.y),[],cat(1,fit1s.frame),'.','SizeData',1); hold on;
        scatter(cat(1,fit2s.x),cat(1,fit2s.y),[],cat(1,fit2s.frame),'+','SizeData',1);
        xx = pntArray1(:,:,1); xx=xx(:);
        yy = pntArray1(:,:,2); yy=yy(:);
        tt = frameT1(:);
        hold on; scatter(xx,yy,[],tt,'o','SizeData',6); colormap('jet');
        xx = pntArray2(:,:,1); xx=xx(:);
        yy = pntArray2(:,:,2); yy=yy(:);
        hold on; scatter(xx,yy,[],tt,'s','SizeData',6);
        colorbar;
        text(xy1p(:,1),xy1p(:,2),sNum,'color','y','FontSize',15); hold on;
        text(xy2p(:,1),xy2p(:,2),sNum,'color','w','FontSize',14); hold on;
        f4.Position(1:2) = [0,0];
        %------

        % --- Merge traces
        % should look at pntArrays for proximity and temporal complimentarity vs.
        % overlap. 
        %   we do this with pdist to allow for arbitary goup size (unlike knn)
        %   and to avoid potential binning splits of adjacent groups. 
        dMap = squareform(pdist(xy1p));
        dMap = triu(dMap);
        dMap(dMap==0) = nan; % ignore self
        dMap(dMap>pars.maxTraceSep) = nan;
        % figure(7); clf; imagesc(dMap); colorbar;  clim([0,10]);
        dMap(dMap<=pars.maxTraceSep) = 1;
        dMap = ~isnan(dMap);
        % figure(7); clf; imagesc(dMap); colorbar; % clim([0,1]);
        hasData = find(sum(dMap,2)>0);
        nGrps = length(hasData);
        grpList = cell(nGrps,1);
        for n=1:nGrps
            currSpot = hasData(n);
            sameGroup = find(dMap(hasData(n),:));
            grpList{n} = [currSpot,sameGroup];
            % remove from map to avoid duplicates
            dMap(sameGroup,:) = false;
            dMap(:,sameGroup) = false;
        end
        grpSize = cellfun(@length,grpList);
        grpList(grpSize<=1) = [];
        


        % Merge traces
        toRemove = false(nPts,2);
        pntArray = {pntArray1,pntArray2};
        pntArrayA = {pntArray1,pntArray2};
        for c=1:2 % loop over colors
            figure(7+c); clf;
            % figure(8); clf;
            nD = length(grpList);
            for d=1:nD
                figure(7+c); subplot(nD,1,d);
                for g=1:length(grpList{d})
                    plot( squeeze(pntArray{c}(grpList{d}(g),:,1)),'.-' ); hold on; 
                    ylabel(grpList{d}(:))
                end
                legend();
                % look for overlap
                traceOlap = pntArray{c}(grpList{d},:,1) > 0;
                % figure(8); subplot(nD,1,d); imagesc(traceOlap);
                % --- merge traces 
                % criteria: 
                %      to merge, the enteries should not overlap more than 1%
                %      both traces should have more than 1% of total
                %      Also, post merge, the magnitude of the step-size should not
                %      increase substantially
                obsPerTrace = sum(traceOlap,2);
                totPoints = sum(obsPerTrace);
                uniqueFrames = sum(traceOlap,1)==1;
                totUnique = sum(uniqueFrames);
                if totUnique/totPoints >(1-pars.maxOverlapFrac) &&  min(obsPerTrace/totPoints) > pars.maxOverlapFrac
                    for g=2:length(grpList{d})
                        mergeArray = cat(4,pntArrayA{c}(grpList{d}(1),uniqueFrames,:),  pntArrayA{c}(grpList{d}(g),uniqueFrames,:));
                        newEntry = nansum(mergeArray,4);  
                        % look at step size variation to decide if we keep this merge.    
                        stp1 = nanmean(abs(diff(pntArrayA{c}(grpList{d}(1),uniqueFrames,1))));
                        stpG = nanmean(abs(diff(pntArrayA{c}(grpList{d}(g),uniqueFrames,1))));
                        stpNew = nanmean(abs(diff(newEntry(1,:,1))));  % stepsize mean of merged distribution 
                        stpOrig = nanmean([stp1,stpG]);    % stepsize mean of original distributions 
                        if stpNew/stpOrig < pars.maxStepVarIncrease
                            pntArrayA{c}(grpList{d}(1),uniqueFrames,:) = newEntry;
                            toRemove(grpList{d}(g),c) = true;
                            figure(7+c); subplot(nD,1,d); 
                            plot( squeeze(pntArrayA{c}(grpList{d}(1),:,1))+2,'.-' ); hold on; 
                            title(['spot channel ',num2str(c)]);
                        end
                    end
                end
            end
        end

 

        % [find(toRemove(:,1));0;find(toRemove(:,2))]
        
        % also remove traces with too few points:
        tooFew = false(nPts,2);
        for c=1:2
            tooFew(:,c) = sum(~isnan(pntArrayA{c}(:,:,1)),2)<pars.minPointsPerTrace; 
        end
        toRemove = toRemove | tooFew;
        
        % actually do the removal
        for c=1:2
            pntArrayA{c}(toRemove(:,c),:,:) = [];
        end
        % now we need to re-pair points
        %    some of the traces may have only merged in one channel
        %    some of the traces may have only been too sparse in one channel

% pntArrayB is re-paired version of pntArrayA
        % === Pair / match traces 
        %  previous step matched based on all spots. This step re-matches the data
        %  using only the linked spots. It is also necessary because the Linking
        %  alogrithm is not guarenteed to preserve order [though if it doesn't drop
        %  unlinked starting points, I think I could force it to prserve order].
        xy1 = squeeze(nanmean(pntArrayA{1}(:,:,1:2),2)); % average over frames
        xy2 = squeeze(nanmean(pntArrayA{2}(:,:,1:2),2)); % average over frames
        [matched,cost] = MatchPoints(xy1,xy2,'maxDist',pars.maxSep);
        m = matched(cost<pars.maxSep,:);      
        nMatched = size(m,1);
        pntArrayB{1} = pntArrayA{1}(m(:,1),:,:);
        pntArrayB{2} = pntArrayA{2}(m(:,2),:,:);
        xy1b = xy1(m(:,1),:);
        xy2b = xy2(m(:,2),:);  
        nPts = size(xy2b,1);
        sNum = cellstr(num2str( (1:nPts)' ));


        % 
        % figure(3); clf; Ncolor(overlayImage);  axis image;
        % hold on; text(xy2b(:,1),xy2b(:,2),sNum,'color','w');
        % hold on; plot(pntArrayB{1}(:,:,1),pntArrayB{1}(:,:,2),'y+');
        % hold on; plot(pntArrayB{2}(:,:,1),pntArrayB{2}(:,:,2),'gs');

        % === noise filter
        % Reduce time resolution to improve inference.  
        % maxGap = 10;
        % smoothWindow = 4; % 3 would integrates up and down 1 point, so the time resolution is effectively cut in 1/2. 
                         % 4 matches the z-scan duration of this dataset 
        pntArrayC = pntArrayB;     
        if pars.smooth
            for p=1:size(pntArrayC{c},1) 
                for d =1:2
                    for c=1:2
                        x1 = squeeze(pntArrayB{c}(p,:,d));
                        xx1 = fillmissing(x1,'linear','maxGap',pars.maxGap);
                        xs1 = smooth(x1,pars.smoothWindow,'moving');
                        xs1(isnan(xx1)) = nan;
                        pntArrayC{c}(p,:,d) = xs1;
                    end
                end
            end
        end
        % (maybe we should do this filtering all all pntArrays)? 
        % === re-order spots by data depth
        isGood = ~isnan(pntArrayC{1}(:,:,1) - pntArrayC{2}(:,:,1));
        [nGood,idx] = sort(sum(isGood,2),'descend');
        pntArrayC{1} = pntArrayC{1}(idx,:,:);
        pntArrayC{2} = pntArrayC{2}(idx,:,:);
        % remove empty traces
        pntArrayC{1}(nGood<pars.minPointsPerTrace,:,:) = [];
        pntArrayC{2}(nGood<pars.minPointsPerTrace,:,:) = [];    
        % apply to other data as well
        xy1c = xy1b(idx,:);
        xy2c = xy2b(idx,:);
        xy1c(nGood<pars.minPointsPerTrace,:) = [];
        xy2c(nGood<pars.minPointsPerTrace,:) = [];
        pntArrayB{1} = pntArrayB{1}(idx,:,:);
        pntArrayB{2} = pntArrayB{2}(idx,:,:);
        pntArrayB{1}(nGood<pars.minPointsPerTrace,:,:) = [];
        pntArrayB{2}(nGood<pars.minPointsPerTrace,:,:) = [];       
        
        figure(2); clf;
        ax1 = subplot(2,2,1); imagesc(pntArrayB{1}(:,:,1)); title('chn 1 detected')
        ax2 = subplot(2,2,2); imagesc(pntArrayB{2}(:,:,1)); title('chn 2 detected')
        ax3 = subplot(2,2,3); imagesc(pntArrayC{1}(:,:,1)); title('chn 1 inferred')
        ax4 = subplot(2,2,4); imagesc(pntArrayC{2}(:,:,1)); title('chn 2 inferred')
        linkaxes([ax1,ax2,ax3,ax4],'xy')

        [nCells,nFrames,~]=size(pntArrayB{1});
        ts = 1:nFrames;
        isChn1 =  ~isnan(pntArrayB{1}(:,ts,1));
        isChn2 = ~isnan(pntArrayB{2}(:,ts,1));
        cnts =[sum(isChn1,2),sum(isChn2,2),sum(isChn1 & isChn2,2)];
        figure(1); clf;
        subplot(2,1,1); bar(cnts); legend('chn1','chn2','both')
        
        isChn1 =  ~isnan(pntArrayC{1}(:,ts,1));
        isChn2 = ~isnan(pntArrayC{2}(:,ts,1));
        cnts =[sum(isChn1,2),sum(isChn2,2),sum(isChn1 & isChn2,2)];
        subplot(2,1,2);  bar(cnts); legend('chn1-inf','chn2-inf','both-inf')


        % % ====== Linked and matched
        % f4 = figure(4); clf;
        % Ncolor(overlayImage); axis image; hold on;
        % plot(cat(1,fits1.x),cat(1,fits1.y),'y+','MarkerSize',3); hold on;
        % plot(cat(1,fits2c.x),cat(1,fits2c.y),'gs','MarkerSize',3);
        % hold on; plot(pntArrayC{1}(:,:,1)',pntArrayC{1}(:,:,2)','m.-');
        % hold on; plot(pntArrayC{2}(:,:,1)',pntArrayC{2}(:,:,2)','c.-');
        % text(xy1b(:,1),xy1b(:,2),sNum,'color','r'); hold on;
        % text(xy2b(:,1),xy2b(:,2),sNum,'color','c'); hold on;
        % savefig(f4,[data.saveFolder,'FOV_Link_Overview.fig'],'compact');
        % xlim([1340,1480]); ylim([260,390]); % zoom in for fun
        % xlim([0,2304]); ylim([0,2304]);
        


        %------ Show Results
        f5 = figure(5); clf;
        f5.Name = 'Merged traces';
        nFrames = size(pntArrayC{1},2);
        frameT = (1:nFrames);
        
        Ncolor(overlayImage);  axis image; hold on;
        %--- convert fits to 
        fit1s = fits1;  
        fit2s = fits2c; 
        for s=1:length(fit1s)
            fit1s(s).frame = repmat(fit1s(s).frame,length(fit1s(s).x),1);
            fit2s(s).frame = repmat(fit2s(s).frame,length(fit2s(s).x),1);
        end
        scatter(cat(1,fit1s.x),cat(1,fit1s.y),[],cat(1,fit1s.frame),'.','SizeData',1); hold on;
        scatter(cat(1,fit2s.x),cat(1,fit2s.y),[],cat(1,fit2s.frame),'+','SizeData',1);
        xx = pntArrayC{1}(:,:,1); xx=xx(:);
        yy = pntArrayC{1}(:,:,2); yy=yy(:);
        nPts = size(xy1c,1);
        frameT1 = repmat(frameT,nPts,1);
        tt = frameT1(:);
        hold on; scatter(xx,yy,[],tt,'o','SizeData',6); colormap('jet');
        xx = pntArrayC{2}(:,:,1); xx=xx(:);
        yy = pntArrayC{2}(:,:,2); yy=yy(:);
        hold on; scatter(xx,yy,[],tt,'s','SizeData',6);
        colorbar;
        sNum = cellstr(num2str( (1:nPts)' ));
        text(xy1c(:,1),xy1c(:,2),sNum,'color','y','FontSize',15); hold on;
        text(xy2c(:,1),xy2c(:,2),sNum,'color','w','FontSize',14); hold on;
        f4.Position(1:2) = [0,0];
        %------


        output.FinalCenterXY = {xy1c,xy2c};
        output.pntArrayB = pntArrayB; % raw traces 
        output.pntArrayC = pntArrayC; % smoothed traces
        
    end
end

%% Fit-Z
if gui_step == 4
    if isempty(pars)
        defaults = cell(0,3);
        % z-fitting
        defaults(end+1,:) = {'nmPerStep','positive',1000};  % width
        defaults(end+1,:) = {'stackHeight','integer',0}; % autocompute
        pars = ParseVariableArguments(varargin,defaults,mfilename);
        pars.directions = 'Fit-Z';
        output = pars;
    else
% ==== z calc
        % maybe it would be better to fit 3D voxel, rather than rely on the
        % fitting the spots independently 
        offFile =  regexprep(data.daxFile1,'_C1.dax','.off');
        offTable = ReadTableFile(offFile,'delimiter',' ');
        pntArrayB = data.pntArrayB;
        pntArrayC = data.pntArrayC;
        [zMat_1b,stackHeight] = FitZTrace(pntArrayB{1},offTable,'parameters',pars); % 
        pntArrayB{1}(:,:,3) = zMat_1b;
        pntArrayC{1}(:,:,3) = zMat_1b;
        zMat_2b = FitZTrace(pntArrayB{2},offTable,'parameters',pars);
        pntArrayB{2}(:,:,3) = zMat_2b;
        pntArrayC{2}(:,:,3) = zMat_2b;
        % save outputs
        output.pntArrayB = pntArrayB; 
        output.pntArrayC = pntArrayC; 
        output.stackHeight = stackHeight;
        disp('step complete')
        % channel2 is mostly not detected in 2 sequential z-steps at this
        % threshold
    end

end

%% Save data
% (should split off save table from save movies
% ==== save table
if gui_step ==5 
     if isempty(pars)
        defaults = cell(0,3);
        % saveData 
        defaults(end+1,:) = {'saveFolder','string',''}; %  save folder
        defaults(end+1,:) = {'saveTables','boolean',true}; %  
        defaults(end+1,:) = {'verbose','boolean',true}; %  
        defaults(end+1,:) = {'npp','positive',108}; %  nm per pixel xy
        pars = ParseVariableArguments(varargin,defaults,mfilename);
        pars.directions = 'Save spot tables';
        output = pars;
    else
        disp(['writing data tables']);
        if isempty(pars.saveFolder)
            pars.saveFolder = [data.saveFolder];
        end
        [dataFolder,daxName] = fileparts(data.daxFile1);
        daxName = regexprep(daxName,'_C1.dax','');
        
        if pars.saveTables
            pntArrayB = data.pntArrayB;
            [nCells,nFrames,nMeasurements] = size(pntArrayB{1});            
            for c=1:nCells
                x1_nm =  pntArrayB{1}(c,:,1)'*npp;
                y1_nm =  pntArrayB{1}(c,:,2)'*npp;
                z1_nm =  pntArrayB{1}(c,:,3)'; % already in nm
                h1 =  pntArrayB{1}(c,:,5)';
                bkd1 = pntArrayB{1}(c,:,6)';
                t1 = (1:length(x1_nm))';
                x2_nm =  pntArrayB{2}(c,:,1)'*npp;
                y2_nm =  pntArrayB{2}(c,:,2)'*npp;
                z2_nm =  pntArrayB{2}(c,:,3)'; % already in nm
                h2 =  pntArrayB{2}(c,:,5)';
                bkd2 = pntArrayB{2}(c,:,6)';
                t2 = (1:length(x2_nm))';
                distXY_nm = sqrt(sum(  (pntArrayB{1}(c,:,1:2) - pntArrayB{2}(c,:,1:2)).^2, 3))'.*npp;
                ctable = table(x1_nm,y1_nm,z1_nm,h1,bkd1,t1,x2_nm,y2_nm,z2_nm,h2,bkd2,t2,distXY_nm);
                tableName = [pars.saveFolder,daxName,'_traj_',num2str(c,'%03d'),'.txt'];
                writetable(ctable,tableName);
                if pars.verbose
                    disp(['wrote ',tableName]);
                end
            end
            % filtered
            pntArrayB = data.pntArrayC;
            for c=1:nCells
                x1_nm =  pntArrayB{1}(c,:,1)'*npp;
                y1_nm =  pntArrayB{1}(c,:,2)'*npp;
                z1_nm =  pntArrayB{1}(c,:,3)'; % already in nm
                h1 =  pntArrayB{1}(c,:,5)';
                bkd1 = pntArrayB{1}(c,:,6)';
                t1 = (1:length(x1_nm))';
                x2_nm =  pntArrayB{2}(c,:,1)'*npp;
                y2_nm =  pntArrayB{2}(c,:,2)'*npp;
                z2_nm =  pntArrayB{2}(c,:,3)'; % already in nm
                h2 =  pntArrayB{2}(c,:,5)';
                bkd2 = pntArrayB{2}(c,:,6)';
                t2 = (1:length(x2_nm))';
                distXY_nm = sqrt(sum(  (pntArrayB{1}(c,:,1:2) - pntArrayB{2}(c,:,1:2)).^2, 3))'.*npp;
                ctable = table(x1_nm,y1_nm,z1_nm,h1,bkd1,t1,x2_nm,y2_nm,z2_nm,h2,bkd2,t2,distXY_nm);
                tableName = [pars.saveFolder,daxName,'_Filtered_Traj_',num2str(c,'%03d'),'.txt'];
                writetable(ctable,tableName);
                if pars.verbose
                    disp(['wrote ',tableName]);
                end
            end
        end
        output.saveRoot = [pars.saveFolder,daxName];
        output.pntArrayA = pntArrayB; % overwriteA; 
     end
end

if gui_step == 6
    % save movies
    % this version uses a moving box 
    if isempty(pars)
        defaults = cell(0,3);
        % saveData 
        defaults(end+1,:) = {'saveFolder','string',''}; %  save folder
        defaults(end+1,:) = {'verbose','boolean',true}; %  save folder
        defaults(end+1,:) = {'saveMovies','boolean',true}; %  save folder
        defaults(end+1,:) = {'cropRadius','integer',15}; % Radius of image around located spot to crop for movies.
        pars = ParseVariableArguments(varargin,defaults,mfilename);
        pars.directions = 'Save Movies';
        output = pars;
    else
        %% assemble movies
        % load the images
        
        if pars.saveMovies
            disp(['loading raw movie data...']);
            daxName1 = data.daxFile1;
            daxName2 = regexprep(daxName1,'C1.dax','C2.dax');
            alignS = data.alignS;
            pntArrayA = data.pntArrayA;
            [nCells,nFrames,~] = size(pntArrayA{1});
            xf = nan(nCells,1);
            yf = nan(nCells,1);
            r = pars.cropRadius;
            pars.stackHeight = data.stackHeight;
            nT = ceil(nFrames/pars.stackHeight);
            im5D = zeros(2*r+1,2*r+1,pars.stackHeight,2,nT,nCells,'uint16');
            t=0;
            figure(1); clf; 
            k = 1; % counter of each frame in a batch
            b = 1; % counter for the batch number
            batchSize = 100;
           for f = 1:1:nFrames
                    % read movie in chunks.
                    %   when running in parpool this should reduce the number of 
                    %   simultaneous read/write commands trying to pull from the 
                    %   same disk.  
                    if k==1
                        ff = batchSize*(b-1)+k; 
                        fe = min(ff+batchSize-1,nFrames);
                        % disp(ff)
                        im1_dax = ReadDax(daxName1,'startFrame',ff,'endFrame',fe,'verbose',false); 
                        im2_dax = ReadDax(daxName2,'startFrame',ff,'endFrame',fe,'verbose',false);
                        if pars.verbose
                            disp([num2str(100*ff/nFrames),'% movie loaded'])
                        end
                    else
                        if k==batchSize
                            k=0; % reset
                            b=b+1;
                        end
                    end
                     % disp([f,k,b,ff])
                    im1 = im1_dax(:,:,k+1);
                    im2 = im2_dax(:,:,k+1);
        %             im1 = ReadDax(daxName1,'startFrame',f,'endFrame',f,'verbose',false); 
        %             im2 = ReadDax(daxName2,'startFrame',f,'endFrame',f,'verbose',false);
                    im2 = ApplyReg(fliplr(im2),alignS);
                    im3 = cat(3,im1,im2);
                    [h,w] = size(im1);
                    % im3 = IncreaseContrast(im3,'high',.99995,'low',.001); % constant contrast
                    
                    % figure(1); clf; Ncolor(im3);
                o = rem(f,pars.stackHeight); 
                if o==0
                    o=pars.stackHeight; % keeping indexes straight
                elseif o==1
                    t=t+1; % starting new series
                end
                for s=1:nCells % s =4
                    try
                   % update centering
                        x = round(pntArrayA{1}(s,f,1)); % in pixels :)
                        y = round(pntArrayA{1}(s,f,2));
                        if ~isnan(x)
                            xf(s) = x;
                            yf(s) = y;
                        end
                    % plot
                    if ~isnan(xf(s))
                        x1 = max([1,xf(s)-r]);
                        x2 = min([w,xf(s)+r]);
                        y1 = max([1,yf(s)-r]);
                        y2 = min([h,yf(s)+r]);
                        imS = im3(y1:y2,x1:x2,:);
                        xL = x2-x1+1;
                        yL = y2-y1+1;
                        if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                            disp(([y1,y2,x1,x2]))
                            disp(([h,w]))
                            disp('debug here')
                        end
                        im5D(1:yL,1:xL,o,:,t,s) = imS; 
                    end
                    catch er
                        warning(['error processing spot ',num2str(s), ' frame ',num2str(f)])
                        disp(([y1,y2,x1,x2]))
                        warning(er.message);
                        warning(er.getReport);
                        disp('debug here')
                    end
                end
                % set(gcf,'color','k');
                % pause(.01);
                k=k+1;
            end
            
        
            %% save files
            [dataFolder,daxName] = fileparts(data.daxFile1);
            daxName = regexprep(daxName,'_C1.dax','');
            saveMovieFolder = [data.saveFolder,daxName,'_movies\'];
            disp(['saving cropped movies as i5d...'])
            if ~exist(saveMovieFolder,'dir')
                mkdir(saveMovieFolder);
            end          
            for s=1:nCells
                im5Dout = im5D(:,:,:,:,:,s);
                % [dH,dW,dZ,dC,dT] = size(im5Dout);
                 WriteImage5D(im5Dout,[saveMovieFolder,'spot_',num2str(s,'%03d'),'.i5d'],'lociX',xf(s),'lociY',yf(s));
            end    
        end
        output.saveMoiveFolder = saveMovieFolder;
    end

end
%% example on plotting movie

%  %% Load and View Saved Files
%  s=4;
% [im5,im5_info] = ReadImage5D([saveFolder,'spot_',num2str(s,'%03d'),'.i5d']);
% for t=1:dT  % t=1
%     figure(1); clf; 
%     for z=1:pars.zFrames
%         subplot(1,pars.zFrames,z);
%         im = squeeze(im5(:,:,z,:,t));
%         Ncolor(im); title(t); 
%         hold on;
%         xf(s) = im5_info.lociX;
%         yf(s) = im5_info.lociY;
%         % some troubleshooting plots
%       plot([fits1(f).x]-xf(s)+r+2,[fits1(f).y]-yf(s)+r+2,'ys'); hold on; % all spots 1
%     plot([fits2c(f).x]-xf(s)+r+2,[fits2c(f).y]-yf(s)+r+2,'bs'); hold on; % all spots 2  (to verify linking)
%     plot(pntArray1a(:,f,1)-xf(s)+r+2,pntArray1a(:,f,2)-yf(s)+r+2,'yo'); hold on; % linked spots from all arrays (for cross talk checking)
%     % the key plots 
%     plot(pntArray1a(s,f,1)-xf(s)+r+2,pntArray1a(s,f,2)-yf(s)+r+2,'w+'); hold on; % linked spots 1
%     plot(pntArray2a(s,f,1)-xf(s)+r+2,pntArray2a(s,f,2)-yf(s)+r+2,'b+'); hold on; % linked spots 2
%         title(['f',num2str(f),'  s',num2str(s),'  rb',num2str(pntArray1a(s,f,5),4),'  gb',num2str(pntArray2a(s,f,5),4)],'color','w')
%     end
%     pause(.01); 
% end


