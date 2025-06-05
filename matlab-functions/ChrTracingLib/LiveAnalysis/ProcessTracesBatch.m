% function output = ProcessTracesBatch(varargin)
% 
% v7 is a 'quiet' batch version of v6
% 
%% see ProcessTracesApp 
clc; close all; clear all; startup; % clean start
NAS02_Vol4 = 'L:\';
% dataFolder = [NAS02_Vol4,'Jude\2023-11-28_16kb\'];  bin_tag = '_pars_240207';  % this one ran completely (did we overwrite)?  
% dataFolder = [NAS02_Vol4,'Jude\2023-11-26_16kb\'];  bin_tag = '_pars_240205';%  ran completely for all  
%  dataFolder = [NAS02_Vol4,'Jude\2023-11-22_76kb\'];   bin_tag = '_pars_240131'; 
%  dataFolder = [NAS02_Vol4,'Jude\2023-11-24_76kb\'];  bin_tag = '_pars_240129';   % 
% dataFolder = '[NAS02_Vol4,'Jude\2023-12-01_50kb\]';  bin_tag = '_pars_240210';  
% dataFolder = '[NAS02_Vol4,'Jude\2023-11-30_50kb\]';  bin_tag = ''; % not analyzed 


NAS02_Vol2 = 'J:\';
% dataFolder = [NAS02_Vol2,'Jude\2023-12-04_130kb\'];  bin_tag = '_pars_240205'; % complete  
% dataFolder = [NAS02_Vol2,'Jude\2023-12-06_130kb\'];  bin_tag = '_pars_240208'; % complete    
% dataFolder = [NAS02_Vol2,'Jude\2023-12-18_255kb\'];  bin_tag = '_pars_240211'; % complete    
% dataFolder = [NAS02_Vol2,'Jude\2023-12-20_255kb\']; bin_tag = '_pars_240213'; % complete   
% dataFolder = [NAS02_Vol2,'Jude\2023-12-26_400kb\'];  bin_tag = '_pars_240209'; % complete    
dataFolder = [NAS02_Vol2,'Jude\2023-12-28_400kb\'];  bin_tag = '_pars_240209'; % running    

saveOverlayTraces = false;
showPlots = false;

t_tot = tic;
sampleDax = FindFiles([dataFolder,'sample*C1.dax']);
movie_index = 1;

for movie_index = 1:length(sampleDax)
    saveFolder = SetFigureSavePath([dataFolder,'Analysis_v2\','sample',num2str(movie_index,'%02d')],'makeDir',true);
    
   trajFiles = FindFiles([saveFolder,'*traj*.txt']);
   if length(trajFiles) > 10
       continue
   end

   disp(['analyzing ',sampleDax{movie_index}])

%%
varIn = {'daxFile1',sampleDax{movie_index},...   % this changes for every FOV  
        'bin_tag',bin_tag,...
        'saveFolder', saveFolder,...
        'alignment_file', [dataFolder,'alignmentData.txt'],...
        'chrom_correct_file',[dataFolder,'tform3D_chromatic.mat']};

%  UPDATE US! 

for morePars = 1  % shorthand, allows this to be collapsed for easier reading 

    defaults = cell(0,3);
    % extra figures to monitor
    defaults(end+1,:) = {'Fig_zLink','integer',0}; %#ok<*SAGROW> % 0 for off
    defaults(end+1,:) = {'Fig_zFit','integer',0}; % 0 for off
    defaults(end+1,:) = {'fig_mergeT','integer',0}; % 0 for off
    defaults(end+1,:) = {'saveMovies','boolean',false}; %  save folder
    % figures for xy-link
    defaults(end+1,:) = {'movAveFig','integer',0}; % per trace
    defaults(end+1,:) = {'spotMapFig','integer',0}; % per trace
    defaults(end+1,:) = {'traceFig','integer',0}; %  per trace

    % ID pairs / seed points 
    defaults(end+1,:) = {'maxSep','nonnegative',20}; % fraction of frames to use 
    defaults(end+1,:) = {'seedBinResolution','integer',7}; % 
    defaults(end+1,:) = {'minSpotsPerBin','integer',50}; % 
    defaults(end+1,:) = {'minSpotsPerTrace','integer',800}; % 
    defaults(end+1,:) = {'imSize','array',[]}; %    
    defaults(end+1,:) = {'removeDots','positive',[]}; % manually ID some spots to remove by index

    % xy-Linking
defaults(end+1,:) = {'parallel','boolean',true}; % 
defaults(end+1,:) = {'movAveSteps','positive',10}; % Number of windows in which to split the trace in when computing the moving average. Note, multiplied by z-depth    
defaults(end+1,:) = {'movAveStepMaxFold','positive',3}; % max fold change in step size relative to local average, used by Remove Jumps
defaults(end+1,:) = {'moveAveMaxPixStep','positive',7}; % max absolute step size (in pixels) for moving avearage trace, used by Remove Jumps
defaults(end+1,:) = {'movAveLocalStep','positive',6}; % window size to determine local jumps  (both rounds use this filter)
% parameters for coarse linking of whole field
defaults(end+1,:) = {'movAveMaxGap','positive',4}; % Downsampling scale to determine spot clusters (spots more than this number of pixels apart will be considered as separate clusters)    
defaults(end+1,:) = {'minObsPerAve','positive',10}; % minimum number of observations in the moving average window to be counted as a valid obs in the coarse downsampling
defaults(end+1,:) = {'maxDistToSeed','positive',35}; % maximum distance to seed point to be counted in the coarse linking
% parameters for refined measurement of moving average
defaults(end+1,:) = {'windowR','positive',35}; % 'radius' of window around the original seed point in which valid localizations may occur
defaults(end+1,:) = {'maxDistFromRoughAve','positive',5}; % Distance in pixels from the rough moving average trace in constructing the refined trace (no binning this time, just a restricted xy window)
defaults(end+1,:) = {'minObsPerAve2','positive',6}; % minimum number of observations in the moving average window to be counted as a valid obs
defaults(end+1,:) = {'coarseInterpWindow','positive',20}; % window size of points used for interpolation filter (allows recovery in fine trace of time blocks not detected in coarse trace)
% parameters for trace assembly
defaults(end+1,:) = {'maxDistFromMovingAve','positive',4}; % Max distance an acceptable spot can be from its local-time-moving average centroid reported as fold-change relative to moving-average stepsize
defaults(end+1,:) = {'maxStep','positive',1.5}; % distance in pixels, the max distance from the last localization to be added to the trace (permisive)   
defaults(end+1,:) = {'maxStepFrame1','positive',30}; % distance in pixels, the max distance from the starting point to start counting (starting point is first non-nan in the time average)


   % merge sister traces xyz/C 
   pars.figSisLabel = 0; % plot demultiplexed sisters  5

       % Merging traces over time
    defaults(end+1,:) = {'maxTraceSep','positive',30};  %       maxTraceSep = 30;     
    defaults(end+1,:) = {'maxOverlapFrac','nonnegative',.01}; %  maxOverlapFrac  = 0.01;
    defaults(end+1,:) = {'maxStepVarIncrease','positive',2};  % fold change increase in step-size variation after trace merge maxStepVarIncrease = 2; 
    defaults(end+1,:) = {'minPointsPerTrace','positive',100};  %  minPointsPerTrace = 100;
    
   
    % save data-tables
    defaults(end+1,:) = {'saveTables','boolean',true}; %  
    defaults(end+1,:) = {'npp','positive',108}; %  nm per pixel xy
    
    % save Movies
    defaults(end+1,:) = {'cropRadius','integer',15}; % Radius of image around located spot to crop for movies.

    % filepath
    defaults(end+1,:) = {'saveFolder','string',''}; 
    defaults(end+1,:) = {'daxFile1','string',''}; 
    defaults(end+1,:) = {'alignment_file','string',''}; % 
    defaults(end+1,:) = {'chrom_correct_file','string',''}; % 
    defaults(end+1,:) = {'bin_tag','string','_2d_iters'}; % 
    defaults(end+1,:) = {'framesToLoad','integer',0}; % 0 for load all frames
    defaults(end+1,:) = {'zDepth','integer',0};  % 0 to autodetect   5
    defaults(end+1,:) = {'loadOverlayMovie','boolean',false};
    defaults(end+1,:) = {'verbose','boolean',true};
    defaults(end+1,:) = {'veryverbose','boolean',false};
    
    
   
    
    pars = ParseVariableArguments(varIn,defaults,mfilename);
    pars.saveOverlayTraces = saveOverlayTraces;

end

%% Main function

%% stp 1, load data
for stp1 = 1   % shorthand, allows this to be collapsed for easier reading
        disp('running...')
        % putting the uigetfiles in the Execute is better
        if isempty(pars.daxFile1)
            [daxName,dataFolder] = uigetfile('*.dax','select a dax movie to load');
            daxFile1 = [dataFolder,daxName];
            
        else
            daxFile1 = pars.daxFile1;
            [dataFolder,~] = fileparts(daxFile1);
        end
        if ~strcmp(dataFolder(end),filesep)
            dataFolder = [dataFolder,filesep]; %#ok<AGROW>
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

        % Auto-compute the z-scan cycle from the off-set file
       
        offFile =  regexprep(daxFile1,'_C1.dax','.off');
        offTable = ReadTableFile(offFile,'delimiter',' ');
        stageZ = offTable.offset;
            % auto determine cycle
         if pars.zDepth == 0
                currZ = 0; oldZ = -inf; 
                [~,z0] =min(stageZ(1:100));
                z = z0;
                nRounds = 100;
                zDep = nan(nRounds,1);
                for n=1:nRounds
                    while currZ > oldZ
                        oldZ = stageZ(z);
                        z=z+1;
                        currZ = stageZ(z);
                    end
                    zDep(n)= z-z0;
                    z0 = z; currZ = 0; oldZ = -inf;  % reset
                end
                pars.zDepth = median(zDep);
           disp(['system determines the stack height = ',num2str(pars.zDepth)]);
           % figure(10); clf; 
           % subplot(2,1,1); plot(stageZ(1:100),'.-');
           %  title(['stackHeight = ',num2str(pars.zDepth)'])
           %  subplot(2,1,2); plot(stageZ,'.-');
        end


        daxFile2 = regexprep(daxFile1,'C1','C2'); %  [dataFolder,'36mW_0001_C2.dax'];
        binName1 = regexprep(daxFile1,'.dax',[pars.bin_tag,'.hdf5']); % [dataFolder,'36mW_0001_C1_2d_iters.hdf5'];  % 2d no iters
        binName2 = regexprep(daxFile2,'.dax',[pars.bin_tag,'.hdf5']);
       

        % === load the images 
        % note that instead of using the first image, we use 
        f = floor(height(offTable)/2)+1; % 1
        [im1f,info1] = ReadDax(daxFile1,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);  % note the +3 since its a 4 step series 
        [im2f,info2] = ReadDax(daxFile2,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);
        im1 = max(im1f,[],3);
        im2 = max(im2f,[],3);
        im2 = ApplyReg(fliplr(im2),alignS);

        overlayImage0 = cat(3,im1,im2); % higher contrast 
        overlayImage0 = IncreaseContrast(overlayImage0,'high',.9999,'low',.1);
        figure(1); clf; Ncolor(overlayImage0);  axis image;

        % === Load the fit data
        if pars.framesToLoad <= 0 
            nFrames = info1.number_of_frames;  %
        else
            nFrames = pars.framesToLoad;
        end
        pars.nFrames = nFrames; % save for letter

        % nFrames = 200; % shorten for testing
        if pars.verbose
            disp('loading channel 1 fits:')
        end    
        fits1 = LoadDaoFits(binName1,'verbose',pars.verbose)

        if pars.verbose
            disp('loading channel 2 fits:')
        end

        fits2 = LoadDaoFits(binName2,'verbose',pars.verbose) 


        % === Apply camera alignment and chromatic correction
        fits2c = Register2CamFitTable(fits2,'alignment_file',pars.alignment_file,'chrom_correct_file',pars.chrom_correct_file);

        % ==== Load movie for overlays
        %   This loads T frames which the next step will crop to and
        %   assemble as a movie with a slider bar.  This can be a useful
        %   way to browse the data, but is unrelated to the analysis. 
        for p=1:double(pars.loadOverlayMovie) %  (made an -if statement into a loop as matlab script editor has a neat function to hide loops, but not if statements).
            % Load an overlay movie 
            T = 50; % number of frames
            selFrames = floor(linspace(1,nFrames-pars.zDepth,T));
            [h,w,~]=size(overlayImage0);
            overlayMovie = zeros(h,w,2,T,'uint16');
            if pars.verbose
                disp('loading movie frames...')
            end
            for i=1:T
                if pars.verbose
                    disp([num2str(100*(i)/T),'% complete']);
                end
                f=selFrames(i);
                im1 = ReadDax(daxFile1,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);  % note the +3 since its a 4 step series 
                im2 = ReadDax(daxFile2,'startFrame',f,'endFrame',f+pars.zDepth,'verbose',false);
                im1 = max(im1,[],3);
                im2 = max(im2,[],3);
                im2 = ApplyReg(fliplr(im2),alignS);
                overlayMovie(:,:,1,i) = IncreaseContrast(im1,'high',.99999,'low',.1);
                overlayMovie(:,:,2,i) = IncreaseContrast(im2,'high',.99999,'low',.1);
            end
            
            if showPlots
                for f=1:T  % play movie backwards (to end on brightest frame); 
                    figure(1); clf;
                    Ncolor(overlayMovie(:,:,:,T-f+1));  axis image;
                    pause(.01);
                end
            end
        end
end 

%%
% fig1= MovieSlider(overlayMovie);

%% step 2: ID pairs
% 
% ID pairs
for stp2 = 1
        [h_im,w_im,~]=size(overlayImage0);
        % === auto-select centers for dancing spots 
        [xy1,ptMap1] = TraceCenters({fits1.x},{fits1.y},...
            'minSpotsPerBin',pars.minSpotsPerBin,...
            'minSpotsPerTrace',pars.minSpotsPerTrace,...
            'binResolution',pars.seedBinResolution,...
             'imSize',[h_im,w_im],...
            'showPlots',0);  % ,'showPlots',true
        [xy2,ptMap2] = TraceCenters({fits2c.x},{fits2c.y},...
            'minSpotsPerBin',pars.minSpotsPerBin,...
            'minSpotsPerTrace',pars.minSpotsPerTrace,...
            'binResolution',pars.seedBinResolution,...
            'imSize',[h_im,w_im],...
            'showPlots',0);

        if showPlots
            figure(4); clf;
            Ncolor(IncreaseContrast(cat(3,ptMap1,ptMap2),'high',.999));
            hold on; plot(xy1(:,1)/pars.seedBinResolution,xy1(:,2)/pars.seedBinResolution,'ro','MarkerSize',20)
            hold on; plot(xy2(:,1)/pars.seedBinResolution,xy2(:,2)/pars.seedBinResolution,'co','MarkerSize',20)
        end

         
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
            if exist('overlayMovie','var')
                tLapseMov = max(overlayMovie,[],4);
            else
                tLapseMov = overlayImage0;
            end
            f1 = figure(1); clf; Ncolor(tLapseMov);  axis image;
            hold on; plot(xy1(:,1),xy1(:,2),'ro','MarkerSize',20);
            hold on; plot(xy2(:,1),xy2(:,2),'bs','MarkerSize',20);
    
            hold on; plot(xy1p(:,1),xy1p(:,2),'yo','MarkerSize',20);
            hold on; plot(xy2p(:,1),xy2p(:,2),'gs','MarkerSize',20);
    
            
            sNum = cellstr(num2str((1:nMatched)'));
            text(xy1p(:,1)+25,xy1p(:,2),sNum,'color','w'); hold on;
            f1.Position = [0,0,1200,1200]; pause(.1);
    
            SetFigureSavePath(pars.saveFolder,'makeDir',true); 
            [~,movieName] = fileparts(sampleDax{movie_index});
            SaveFigure(f1,'name',['initPairs_',movieName],'formats',{'png'},'overwrite',true);
        
        % --------

end


%% Step 3 - LinkXY

for stp3 = 1    
    % the same analysis is applied to both color channels
    % for transparency and troubleshooting, I've written this out
    t_link = tic;
    SetFigureSavePath(pars.saveFolder);
      [spotTrace1,spotTrace2,spotsPerStep] = LinkXYdoublets(fits1,fits2c,xy1p,xy2p,pars.zDepth,...
       'nFrames',nFrames,'parameters',pars);
      % 'selSpots',0,'tLapseMov',tLapseMov, ,'pause',false,'saveFigs',true ...
      %  'parallel',false, 'movAveFig',0,'traceFig',0,'parallel',true
    t_link = toc(t_link);
    disp(['xy link complete in ',num2str(t_link),' seconds'])
end
[nSpots,tObs,dDims] = size(spotTrace1);
zDepth = pars.zDepth;


%% Step 4: Assign sisters in XYZ and match sisters across color channels
 spotTrace1_out = spotTrace1;
 spotTrace2_out = spotTrace2;
  % view the traces
   dims = {'x','y'};
   sisRatios = zeros(nSpots,2);
   for s=1:nSpots  % s=3 % 88 s=91  s=31   s=4  s=33  s=48
    % detect sisters
        sisRatio_c1 = DetectSisters(squeeze(spotTrace1(s,:,:,:)),'zDepth',pars.zDepth);
        sisRatio_c2 = DetectSisters(squeeze(spotTrace2(s,:,:,:)),'zDepth',pars.zDepth);
        sisRatios(s,:) = [sisRatio_c1,sisRatio_c2];

        if showPlots
            figure(1); clf; % s=101 % - clear doublet
            for d=1:2
                subplot(2,1,d);
                plot(squeeze(spotTrace1(s,:,d,1)),'.'); hold on;
                plot(squeeze(spotTrace1(s,:,d,2)),'.')
                plot(squeeze(spotTrace2(s,:,d,1)),'.'); hold on;
                plot(squeeze(spotTrace2(s,:,d,2)),'.'); ylabel(dims{d})
                legend('1a','1b','2a','2b'); title(['doublet trace spot ',num2str(s), ' sis rs=',num2str(sisRatio_c1,3),' ',num2str(sisRatio_c2,3)]);
                pause(.01);
            end
            if sisRatio_c1 > 0.2 || sisRatio_c2 > 0.2
                pause(1);
                [spotTrace1_out(s,:,:,:),spotTrace2_out(s,:,:,:)] = AssignSistersXYZC(...
                    squeeze(spotTrace1(s,:,:,:)),squeeze(spotTrace2(s,:,:,:)),'traceNum',s,...
                    'figSisLabel',pars.figSisLabel) ;
            end
        end
   end


%% Step 5: Fit Z
% this is a bit slow 
spotStacks = {spotTrace1_out,spotTrace2_out};
for stp5 = 1
    t_zfit = tic;
    pntArray = cell(1,2);
    for c=1:2
        pntArray_cs = cell(nSpots,1); % a cell array makes this easier / more mem-efficient to parfor loop  
        spotStacks_s = cell(nSpots,1);
        for s=1:nSpots % prep the parfor loop
            spotStacks_s{s} =  squeeze(spotStacks{c}(s,:,:,1)); % only take sister 1
        end
        parfor s=1:nSpots %  parfor this  % s=33
            % sort array into z-stack
            currStack = nan(nFrames/zDepth,8,zDepth);  %  a 3D matrix, t_Obs x d_Dims x z_Depth. 
            for z =1:zDepth
                currStack(:,:,z) =spotStacks_s{s}(z:zDepth:end,:) ; 
            end
            % fitStack = currStack;
            fitStack = FitZ(currStack,'spotNum',s,'showFig',0); % pars.Fig_zFit
            % sort z-stack back into an array
            pntArray_cs{s} = nan(nFrames,8);
            for z=1:zDepth
                pntArray_cs{s}(z:zDepth:end,:) = fitStack(:,:,z);
            end
        end
        % reorganize into array
        pntArray{c} = permute(cat(3,pntArray_cs{:}),[3,1,2]);
    end
    t_zfit = toc(t_zfit);
    disp(['z link complete in ',num2str(t_zfit),' seconds']);
end


%% Step 6: ID traces to merge in time
%  pntArray are the original traces
%  pntArrayA contains the the merged traces
for stp6 = 1
    % figure(3); clf; Ncolor(tLapseMov);  axis image;
    % sNum = cellstr(num2str( (1:size(xy1p,1))' ));
    % text(xy1p(:,1),xy1p(:,2)+5,sNum,'color','y','FontSize',12); hold on;
  

    % --- Merge traces
    % should look at pntArrays for proximity and temporal complimentarity vs.
    % overlap. 
    %   we do this with pdist to allow for arbitary goup size (unlike knn)
    %   and to avoid potential binning splits of adjacent groups. 
    % find seed points that are nearby
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
    toRemove = false(nSpots,2);
    pntArrayA = pntArray;
    for c=1:2 % loop over colors
        if pars.fig_mergeT
            figure(7+c); clf;
        end
        % figure(8); clf;
        nD = length(grpList);
        for d=1:nD  % d=19
            if pars.fig_mergeT
                figure(7+c); clf; subplot(2,1,1); % subplot(nD,1,d);
                for g=1:length(grpList{d})
                    plot( squeeze(pntArray{c}(grpList{d}(g),:,1)),'.-' ); hold on; 
                    ylabel(grpList{d}(:))
                end
                legend();
            end
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
                        if pars.fig_mergeT
                            figure(7+c);  subplot(2,1,2); 
                            plot( squeeze(pntArrayA{c}(grpList{d}(1),:,1)),'.-' ); hold on;   % +2
                            title(['spot channel ',num2str(c)]); pause(1);
                        end
                    end
                end
            end
        end
    end
    disp(['number of merged traces found for chn1 & chn2 = ',num2str(sum(toRemove))])
            
    
    %--- show merged traces
    f2= figure(2); clf; Ncolor(tLapseMov);  axis image;
    sNum = cellstr(num2str( (1:size(xy1p,1))' ));
    keep = ~toRemove(:,1);
    text(xy1p(keep,1),xy1p(keep,2),sNum(keep),'color','w','FontSize',15); hold on;
    text(xy1p(~keep,1),xy1p(~keep,2),sNum(~keep),'color','r','FontSize',15); hold on;
    SaveFigure(f2,'name',['removedPairs_',movieName],'formats',{'png'},'overwrite',true);
     
% for trouble-shooting
   % add overlay of colored spots; 
    xx = pntArray{1}(:,:,1); xx=xx(:);
    yy = pntArray{1}(:,:,2); yy=yy(:);
    tt = (1:length(xx))';
    scatter(xx,yy,[],tt,'o','SizeData',6); colormap('jet'); colorbar;

    xx = pntArray{2}(:,:,1); xx=xx(:);
    yy = pntArray{2}(:,:,2); yy=yy(:);
    scatter(xx+15,yy+15,[],tt,'s','SizeData',6);
    


    % also remove traces with too few points:
    tooFew = false(nSpots,2);
    for c=1:2
        tooFew(:,c) = sum(~isnan(pntArrayA{c}(:,:,1)),2)<pars.minPointsPerTrace; 
    end
    disp(['removing from chn1 & chn2 ' , num2str(sum(tooFew)), ' traces with fewer than ',num2str(pars.minPointsPerTrace),' observations']);
    toRemove = toRemove | tooFew;
    
    % % remove pairs        
    % f2= figure(2); clf; Ncolor(tLapseMov);  axis image;
    % sNum = cellstr(num2str( (1:size(xy1p,1))' ));
    % keep = ~toRemove(:,1);
    % text(xy1p(keep,1),xy1p(keep,2),sNum(keep),'color','c','FontSize',15); hold on;
    % text(xy1p(~keep,1),xy1p(~keep,2),sNum(~keep),'color','r','FontSize',15); hold on;
    % SaveFigure(f2,'name',['remvoedPairs_',movieName],'formats',{'png'});

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
end

%% Step 7: fill missing (optional?)
for stp7 = 1
    pntArrayC = pntArrayB;     
    for p=1:size(pntArrayC{c},1) 
        for d =1:2
            for c=1:2
                x1 = squeeze(pntArrayB{c}(p,:,d));
                xx1 = fillmissing(x1,'linear','maxGap',1*pars.zDepth);
                pntArrayC{c}(p,:,d) =xx1;
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
end       

%% Step 8: Save data tables
% ==== save table
for stp8 = 1

  %==== run step
  npp= pars.npp;
        disp(['writing data tables']);
        [~,daxName] = fileparts(daxFile1);
        daxName = regexprep(daxName,'_C1.dax','');
        
        if pars.saveTables
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
                if pars.veryverbose
                    disp(['wrote ',tableName]);
                end
            end
            % filtered

            for c=1:nCells
                x1_nm =  pntArrayC{1}(c,:,1)'*npp;
                y1_nm =  pntArrayC{1}(c,:,2)'*npp;
                z1_nm =  pntArrayC{1}(c,:,3)'; % already in nm
                h1 =  pntArrayC{1}(c,:,5)';
                bkd1 = pntArrayC{1}(c,:,6)';
                t1 = (1:length(x1_nm))';
                x2_nm =  pntArrayC{2}(c,:,1)'*npp;
                y2_nm =  pntArrayC{2}(c,:,2)'*npp;
                z2_nm =  pntArrayC{2}(c,:,3)'; % already in nm
                h2 =  pntArrayC{2}(c,:,5)';
                bkd2 = pntArrayC{2}(c,:,6)';
                t2 = (1:length(x2_nm))';
                distXY_nm = sqrt(sum(  (pntArrayC{1}(c,:,1:2) - pntArrayC{2}(c,:,1:2)).^2, 3))'.*npp;
                ctable = table(x1_nm,y1_nm,z1_nm,h1,bkd1,t1,x2_nm,y2_nm,z2_nm,h2,bkd2,t2,distXY_nm);
                tableName = [pars.saveFolder,daxName,'_Filtered_Traj_',num2str(c,'%03d'),'.txt'];
                writetable(ctable,tableName);
                if pars.veryverbose
                    disp(['wrote ',tableName]);
                end
            end
        end

end
 


%%  Step9: save some image stats
if saveOverlayTraces
    for stp9 =  1
    
        % ==== create & save a tiled 2-color image array of the spot-data 
        %   using the the middle frame 
        T = size(overlayMovie,4);
        tLapseMov = overlayMovie(:,:,:,floor(T/2)); % middle frame
        for p=0:double(pars.loadOverlayMovie) %   (made an -if statement into a loop as matlab script editor has a neat function to hide loops, but not if statements).  
            nP = size(xy1c,1);
            xf = round(xy1c(:,1));
            yf = round(xy1c(:,2));
            r= pars.spotDisplayRadius;
            [h,w,~] = size(tLapseMov); % (this one is big, we keep it in 'data' to avoid duplication)   
            imTile1 = zeros(r,r,nP,'uint16');
            imTile2 = zeros(r,r,nP,'uint16');
            for p=1:nP % loop over points
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
                im1 = squeeze(tLapseMov(y1:y2,x1:x2,1));
                im2 = squeeze(tLapseMov(y1:y2,x1:x2,2));
                im1 = IncreaseContrast(im1,'high',.99999,'low',.1);
                im2 = IncreaseContrast(im2,'high',.99999,'low',.1);
                imTile1(1:yL,1:xL,p) = im1;
                imTile2(1:yL,1:xL,p) = im2;
            end
            % the display is the slow
            t=1;
            [~,labelOffsets,imT1] = TileImage(imTile1(:,:,:),'colormap',gray(256),'numRows',10);
            [~,~,imT2] = TileImage(imTile2(:,:,:),'colormap',gray(256),'numRows',10);            
            imO = Ncolor(cat(3,imT1,imT2));
            f1 = figure(1); clf; imagesc(imO);
            nums = cellstr(num2str( (1:nP)' ));
            text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
        end
        f1.Position = [0 0 1200 1200];  pause(.1); % keep a constant size for export;
        SaveFigure(f1,'name',['singlePointRC_',movieName],'formats',{'png'},'overwrite',true);
    
    
        % === create a tiled 2-color superposition of all time as an array
        tLapseMov = max(overlayMovie,[],4);
        for p=0:double(pars.loadOverlayMovie) %   (made an -if statement into a loop as matlab script editor has a neat function to hide loops, but not if statements).  
            nP = size(xy1c,1);
            xf = round(xy1c(:,1));
            yf = round(xy1c(:,2));
            r= pars.spotDisplayRadius;
            [h,w,~] = size(tLapseMov); % (this one is big, we keep it in 'data' to avoid duplication)   
            imTile1 = zeros(r,r,nP,'uint16');
            imTile2 = zeros(r,r,nP,'uint16');
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
                im1 = squeeze(tLapseMov(y1:y2,x1:x2,1));
                im2 = squeeze(tLapseMov(y1:y2,x1:x2,2));
                im1 = IncreaseContrast(im1,'high',.99999,'low',.1);
                im2 = IncreaseContrast(im2,'high',.99999,'low',.1);
                imTile1(1:yL,1:xL,p) = im1;
                imTile2(1:yL,1:xL,p) = im2;
            end
            % the display is the slow
            t=1;
            [~,labelOffsets,imT1] = TileImage(imTile1(:,:,:),'colormap',gray(256),'numRows',10);
            [~,~,imT2] = TileImage(imTile2(:,:,:),'colormap',gray(256),'numRows',10);            
            imO = Ncolor(cat(3,imT1,imT2));
            f1 = figure(1); clf; imagesc(imO);
            nums = cellstr(num2str( (1:nP)' ));
            text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
        end
        f1.Position = [0 0 1200 1200];  pause(.1); % keep a constant size for export;
        SaveFigure(f1,'name',['timeLapseRC_',movieName],'formats',{'png'},'overwrite',true);
        
        % === create a tiled array of the 2D trajectories of each point 
        %   in 2-colors as a simple plot
        [nPts,tObs,dDat] = size(pntArrayB{1});
        gap = 50;
        f2 = figure(2); clf;
        for n=1:nPts
            [rowNum,colNum] = ind2sub([10,ceil(nPts/10)],n);
            x1 = pntArrayB{1}(n,:,1) - nanmean(pntArrayB{1}(n,:,1) ) + labelOffsets(n,1); %#ok<*NANMEAN>
            y1 = pntArrayB{1}(n,:,2) - nanmean(pntArrayB{1}(n,:,2) ) + labelOffsets(n,2);% rowNum*gap;
            x2 = pntArrayB{2}(n,:,1) - nanmean(pntArrayB{2}(n,:,1) ) + labelOffsets(n,1);% + colNum*gap;
            y2 = pntArrayB{2}(n,:,2) - nanmean(pntArrayB{2}(n,:,2) ) + labelOffsets(n,2);% + rowNum*gap;
            plot(x1(~isnan(x1)),y1(~isnan(x1)),'.-','color',[0 .5 0],'linewidth',.1); hold on;
            plot(x2(~isnan(x2)),y2(~isnan(x2)),'.-','color',[.5 0 .5],'linewidth',.1); hold on;
            set(gca,'YDir','reverse')
            % pause(.3);
        end
        f2.Position = [0 0 1200 1400];  pause(.1);  % keep a constant size for export;
        SaveFigure(f2,'name',['traceDataMG_',movieName],'formats',{'png'},'overwrite',true);
    
    
    % === Create 2-pannel tiled array of color-as-time of all the spots
        for p=0:double(pars.loadOverlayMovie)
            nP = size(xy1c,1);
            xf = round(xy1c(:,1));
            yf = round(xy1c(:,2));
            r= pars.spotDisplayRadius;
            [h,w,~] = size(tLapseMov); % (this one is big, we keep it in 'data' to avoid duplication)   
            imTile1 = cell(nP,1); % zeros(r,r,nP,T,'uint16');
            imTile2 = cell(nP,2); % zeros(r,r,nP,T,'uint16');
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
                im3D = zeros(2*r+1,2*r+1,T,'uint16');
                im1 = squeeze(overlayMovie(y1:y2,x1:x2,1,:));
                [hi,wi,~] = size(im1);
                im3D(1:hi,1:wi,:) = IncreaseContrast(im1,'high',.9999,'low',.1);
                imTile1{p} = Ncolor(im3D,'colorMax',true,'colormap','cool');
        
                im3D = zeros(2*r+1,2*r+1,T,'uint16');
                im2 = squeeze(overlayMovie(y1:y2,x1:x2,2,:));
                [hi,wi,~] = size(im2);
                im3D(1:hi,1:wi,:) = IncreaseContrast(im2,'high',.9999,'low',.1);
                imTile2{p} = Ncolor(im3D,'colorMax',true,'colormap','cool');
        
            end
            % figure(3); clf; TileImageStack(cat(4,imTile1{:}));
            imStk1 = cat(4,imTile1{:}); % x y RGB cell
            imStk2 = cat(4,imTile2{:}); % x y RGB cell
        
            imRGB1 = cell(3,1); 
            imRGB2 = cell(3,1); 
            for c=1:3
            [~,labelOffsets,imRGB1{c}] = TileImage(squeeze(imStk1(:,:,c,:)),'colormap',gray(256),'numRows',10);
            [~,labelOffsets,imRGB2{c}] = TileImage(squeeze(imStk2(:,:,c,:)),'colormap',gray(256),'numRows',10);
            end
            rgb1 = cat(3,imRGB1{:});
            rgb2 = cat(3,imRGB2{:});
            nums = cellstr(num2str( (1:nP)' ));
        end
        f3 = figure(3); clf; 
        f3.Position = [0 0 3000 1000];  pause(.1);
        subplot(1,2,1); imagesc(rgb1); text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
        subplot(1,2,2); imagesc(rgb2); text(labelOffsets(:,1),labelOffsets(:,2)+5,nums,'color','w');
        set(gca,'color','k');  pause(.1);
        SaveFigure(f3,'name',['timeLapseColor_',movieName],'formats',{'png'},'overwrite',true);
    
    
    % === For all points, create a time-overlay trace
    if pars.saveOverlayTraces
        for stp9b=1
            % this is a bit slow
            [nPts,tObs,dDat] = size(pntArrayB{1});
            frameT = (1:tObs);
            f3= figure(3); clf;
             r = 20;
            xf = round(xy1c(:,1));
            yf = round(xy1c(:,2));
            for n=1:nPts   
                    x1 = max([1,xf(n)-r]);
                    x2 = min([w,xf(n)+r]);
                    y1 = max([1,yf(n)-r]);
                    y2 = min([h,yf(n)+r]);
                    xL = x2-x1+1;
                    yL = y2-y1+1;
                    if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                        disp(([y1,y2,x1,x2]))
                        disp(([h,w]));
                        disp('debug here')
                    end
                    im3D = zeros(2*r+1,2*r+1,T,'uint16');
                    im1 = squeeze(overlayMovie(y1:y2,x1:x2,1,:));
                    [hi,wi,~] = size(im1);
                    im3D(1:hi,1:wi,:) = IncreaseContrast(im1,'high',.9999,'low',.1);
                    rgb1 = Ncolor(im3D,'colorMax',true,'colormap','cool');
            
                    im3D = zeros(2*r+1,2*r+1,T,'uint16');
                    im2 = squeeze(overlayMovie(y1:y2,x1:x2,2,:));
                    [hi,wi,~] = size(im2);
                    im3D(1:hi,1:wi,:) = IncreaseContrast(im2,'high',.9999,'low',.1);
                    rgb2 = Ncolor(im3D,'colorMax',true,'colormap','cool');
                x1 = pntArrayB{1}(n,:,1) - nanmean(pntArrayB{1}(n,:,1) ) + r; % labelOffsets(n,1); %#ok<*NANMEAN>
                y1 = pntArrayB{1}(n,:,2) - nanmean(pntArrayB{1}(n,:,2) ) + r; %  labelOffsets(n,2);% rowNum*gap;
                x2 = pntArrayB{2}(n,:,1) - nanmean(pntArrayB{2}(n,:,1) ) + r; % labelOffsets(n,1);% + colNum*gap;
                y2 = pntArrayB{2}(n,:,2) - nanmean(pntArrayB{2}(n,:,2) ) + r; %labelOffsets(n,2);% + rowNum*gap;
                f3 = figure(3); clf;
                subplot(1,2,1); imagesc(rgb1); hold on;
                plot(x1(~isnan(x1)),y1(~isnan(x1)),'-','color',[0 .5 0],'linewidth',.1); hold on;
                plot(x2(~isnan(x2)),y2(~isnan(x2)),'-','color',[.5 0 .5],'linewidth',.1); hold on;
                scatter(x1(:),y1(:),[],frameT(:),'o','SizeData',6); colormap('jet'); hold on; set(gca,'YDir','reverse')
                subplot(1,2,2); imagesc(rgb2); hold on;
                plot(x1(~isnan(x1)),y1(~isnan(x1)),'-','color',[0 .5 0],'linewidth',.1); hold on;
                plot(x2(~isnan(x2)),y2(~isnan(x2)),'-','color',[.5 0 .5],'linewidth',.1); hold on;
                scatter(x2(:),y2(:),[],frameT(:),'o','SizeData',6); colormap('jet');  hold on; set(gca,'YDir','reverse')
               f3.Position = [0 0 2000 800];   pause(.1); % keep a constant size for export;
               SaveFigure(f3,'name',['traceDataColorCode_',movieName,'_spot',num2str(n)],'formats',{'png'},'overwrite',true,'verbose',false);
               pause(.1);
            end
        end
    
    end
    end
end
%% Step 10, save mini-movies (optional)
for stp10 =1
        if pars.saveMovies
            disp(['loading raw movie data...']);
            t_sav = tic;
            daxName1 = daxFile1;
            daxName2 = regexprep(daxName1,'C1.dax','C2.dax');
            [nCells,nFrames,~] = size(pntArrayB{1});
            xf = nan(nCells,1);
            yf = nan(nCells,1);
            r = pars.cropRadius;
            nT = ceil(nFrames/pars.zDepth);
            im5D = zeros(2*r+1,2*r+1,pars.zDepth,2,nT,nCells,'uint16');
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
                o = rem(f,pars.zDepth); 
                if o==0
                    o=pars.zDepth; % keeping indexes straight
                elseif o==1
                    t=t+1; % starting new series
                end
                for s=1:nCells % s =4
                    try
                   % update centering
                        x = round(pntArrayB{1}(s,f,1)); % in pixels :)
                        y = round(pntArrayB{1}(s,f,2));
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
            [~,daxName] = fileparts(daxFile1);
            daxName = regexprep(daxName,'_C1.dax','');
            saveMovieFolder = [pars.saveFolder,daxName,'_movies\'];
            disp(['saving cropped movies as i5d...']);
            if ~exist(saveMovieFolder,'dir')
                mkdir(saveMovieFolder);
            end          
            for s=1:nCells
                im5Dout = im5D(:,:,:,:,:,s);
                % [dH,dW,dZ,dC,dT] = size(im5Dout);
                 WriteImage5D(im5Dout,[saveMovieFolder,'spot_',num2str(s,'%03d'),'.i5d'],'lociX',xf(s),'lociY',yf(s));
            end    
            t_sav = toc(t_sav)
            disp(['spent ',num2str(t_sav/60,3),' minutes saving movies']);
        end

end


end
t_tot = toc(t_tot);
disp(['completed in ',num2str(t_tot/60,3),' minutes']);



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