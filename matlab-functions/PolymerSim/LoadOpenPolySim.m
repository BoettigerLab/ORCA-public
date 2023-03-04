function [polyReps,mapReps,lefData,loopData,parsTable] = LoadOpenPolySim(folder,varargin)
%% Outputs:
% polyReps - cell array of replicates, each nMonomers x 3-dims x timesteps
% mapReps - cell array of replicates, each nMonomers x nMonomers x timesteps
% lefData - cell array of replicates x timesteps, each nLEFs x 3-dims x
%            2-heads.  Tracks actual 3D positions of each head through time.  
% loopData - cell array replicates x timesteps x LEFs
% parsTable   - table of simulation properties
% 
%% optional Inputs
%
% 


defaults = cell(0,3); 
% load properties
defaults(end+1,:) = {'polyStep','integer',5}; % downsample the polymer (in units of monomers)
defaults(end+1,:) = {'mapStep','integer',5}; % downsample the distance map (in units of polymer steps, applied after polySteps)
defaults(end+1,:) = {'timeStep','integer',1}; 
defaults(end+1,:) = {'tScale','integer',10}; 
% less used properties
defaults(end+1,:) = {'lenPoly','integer',0}; % should load from file
defaults(end+1,:) = {'rptPoly','integer',1}; % should load from file
defaults(end+1,:) = {'maxBlocks','integer',inf};
defaults(end+1,:) = {'showInfo','boolean',false}; % 

% just for plotting
defaults(end+1,:) = {'showPlots','boolean',false}; % 
defaults(end+1,:) = {'threshold','integer',[]}; % monomer radii for contact map (empty for auto compute from polyStep mapStep)
defaults(end+1,:) = {'distMapFig','integer',11}; % Figure to show distance map in
defaults(end+1,:) = {'contMapFig','integer',12}; % Figure to show contact map in
defaults(end+1,:) = {'TADs','freeType',[]}; % downsample the polymer (in units of monomers)
defaults(end+1,:) = {'trajFig','integer',13}; % figure to show LEF trajectory image
defaults(end+1,:) = {'avePosFig','integer',14}; % figure to show average LEF position (ChIP track)
pars = ParseVariableArguments(varargin,defaults,mfilename);

% try to load parameter tables
if pars.lenPoly==0
    try
        parsTable = readtable([folder,'simPars.txt']);
        pars.TADs = str2num(parsTable{1,5}{1}); %#ok<ST2NM>
        pars.rptPoly = parsTable.replicates(1);
        pars.lenPoly = parsTable.monomers(1);
    catch er
        warning('did not detect a parameter table, using defaults');
    end
else
    parsTable = [];
end

% find data blocks and determine block size. 
blockFiles = FindFiles([folder,'*blocks_*.h5']);
if isempty(blockFiles)
    error(['no data blocks found in ',folder]);
end
[~,blockName] = fileparts(blockFiles{1});
blockVal = strsplit(blockName,'_');
blockSize = -(eval(blockVal{2})-1); % turn 0-9 into 10

% short-hand parameters for easy of use
polyStep = pars.polyStep;
mapStep  = pars.mapStep;
rptPoly  = pars.rptPoly; % 10;
lenPoly = pars.lenPoly;
nB = min([pars.maxBlocks,length(blockFiles)]); % number of blocks detected.  for speed lets not use all of them, allow truncation to  maxBlocks
if length(blockFiles) > pars.maxBlocks % show warning if truncating
    disp(['detected ', num2str(length(blockFiles)) ' maxBlocks ',num2str(pars.maxBlocks)]);
end

%% Populate the data array 
timeStep = min(pars.timeStep,blockSize); % can't be larger than block size
totTime = nB*blockSize/timeStep;
lenTotPoly = lenPoly*pars.rptPoly; % 40e3; % total size of polymer in simulation
lenDispPoly = lenPoly/polyStep;
subset = 1:polyStep:lenTotPoly;
stps = ceil(length(subset)/rptPoly);
polyReps = cell(rptPoly,1);
for p=1:rptPoly
    polyReps{p} = zeros(stps,3,totTime);
end
t=0;
tt=0;
for b=1:nB
    blockStart = num2str(  blockSize*(0+b-1) );
    blockEnd   = num2str(  blockSize*(0+b-1) + (blockSize-1));
    file = [folder,'blocks_',blockStart,'-',blockEnd,'.h5'];
    for e=1:ceil(blockSize/timeStep)
        t=t+timeStep;
        tt=tt+1;
        data = h5read(file,['/',num2str(t-1),'/pos']); % there should be an easy way to load subsets of the data with h5read, just need to look into it
        polyData = data(:,subset)';
        for p=1:rptPoly
            m1 = round((p-1)*lenDispPoly+1);  % 
            m2 = round(p*lenDispPoly);         
            temp = double(polyData(m1:m2,:));
            polyReps{p}(:,:,tt) = temp;
        end
    end
end

%% sort polymers into maps
if nargout > 1
    mapSteps = 1:mapStep:stps;
    nBins = length(mapSteps);
    mapReps = cell(rptPoly,1);
    for p=1:rptPoly
        mapReps{p} = zeros(nBins,nBins,totTime);
        for t=1:totTime
           mapReps{p}(:,:,t) =  squareform(pdist(polyReps{p}(mapSteps,:,t)));
        end
    end
    mapsOut = cat(3,mapReps{:});
end


%% LEF position analysis
% no need to load LEF positions if not required
if nargout > 2

    try
     % requires the actual polymer length
    parsTable = readtable([folder,'simPars.txt']);
    pars.TADs = str2num(parsTable{1,5}{1}); %#ok<ST2NM>
    pars.rptPoly = parsTable.replicates(1);
    lenPoly = parsTable.monomers(1);
    lenPolyLoaded = size(polyReps{1},1);
    catch
    end

    % --- LEF positions
    lefFile =[folder,'LEFPos','.h5'];
    info = h5info(lefFile);
    if pars.showInfo
        disp('Attribute Names:');
        disp({info.Attributes.Name});
        disp('Attribute Values:');
        disp({info.Attributes.Value});
        disp(['Dataset name: ',info.Datasets.Name]);
    end
    lefPos = h5read(lefFile,'/positions');  % 2 x N-extruders x T-timesteps
    numExtruders = info.Attributes(2).Value; % total LEFs across all replicates
    %%
    tStepsSimulated = size(lefPos,3);
    tScale = pars.tScale; % saves per batch in python
    simTimeStep = timeStep*tScale;
    %  tStepsOut = round(tStepsSimulated/simTimeStep);
    
    simTimeStep = tStepsSimulated./size(polyReps{1},3)  % (trying to upgrade to not guessing tscale 


    extruders = 1:numExtruders; % which of the extruders to plot

    % 
    lefPos = double(lefPos);
    ts = simTimeStep:simTimeStep:tStepsSimulated;
    tStepsOut = length(ts);
    lefData = cell(rptPoly,tStepsOut);
    loopData= cell(rptPoly,tStepsOut);
    tt=0;
    try
    for t = ts % t=ts(1);
        tt=tt+1;
        nn = 0;
        for n=1:numExtruders
           rL = floor((lefPos(1,n,t))/lenPoly)+1;% which replicate
           mL = rem(lefPos(1,n,t),lenPoly); % which monomer
           mL = floor(mL/polyStep)+1;
           
           rR = floor((lefPos(2,n,t))/lenPoly)+1;% which replicate
           mR = rem(lefPos(2,n,t),lenPoly); % which monomer
           mR = floor(mR/polyStep)+1;
           if mL > lenPolyLoaded || mR > lenPolyLoaded
               continue;
           end
           if rR ~= rL
              if rR > rL && mR < mL && mL~=1
                  rR = rR-1;
                  mR = lenPoly/polyStep;
               else
                  disp([rL,rR,mL,mR]);
                  disp(lefPos(:,n,t));
              end
           end
           % disp([rR,rL,mL,mR])
               curL = polyReps{rL}(mL,:,tt); % xyz of left extruder
               curR = polyReps{rR}(mR,:,tt); % xyz of right extruder
               newData = cat(3,curL,curR);
               lefData{rR,tt} = cat(1,lefData{rR,tt},newData); % the given replicate has a random number of extruders, much less than total extruder 
               nn = size(loopData{rR,tt},1)+1;
               loopData{rR,tt}{nn,1} = mL:mR ; % polyReps{rL}(mL:mR,:,tt);
        end
    end
    catch er
        warning(er.getReport);
        warning('error parsing LEF data, returning raw positions instead');
        lefData = lefPos; 
    end
    
end

if pars.showPlots
    % images of extrusion
    imLeft = false(tStepsOut,lenTotPoly);
    imRight = false(tStepsOut,lenTotPoly);
    for t=1:tStepsOut
        leftSteps = squeeze(round(lefPos(1,extruders,t)/pars.polyStep));
        leftSteps(leftSteps==0) = lenTotPoly; % just for plotting
        imLeft(t,leftSteps) = true; 
        rightSteps = squeeze(round(lefPos(2,extruders,t)/pars.polyStep));
        rightSteps(rightSteps==0) = lenTotPoly; 
        imRight(t,rightSteps) = true; 
    end
    imTraj = imLeft + 2*imRight; % colors L and R extruders differently.

    % sort repeats into separate blocks
    imLef = cell(rptPoly,1);
    for k=1:rptPoly
        m1 = (k-1)*lenPoly/pars.polyStep+1;  % 
        m2 = k*lenPoly/pars.polyStep;
        imLef{k} = imTraj(:,m1:m2);
    end
end

%% Display results
if pars.showPlots && nargout > 1
    pars.ctcf = ceil(pars.TADs/mapStep/polyStep);
    if pars.distMapFig == 0
        pars.f1 = figure(); 
    else
        pars.f1 = figure(pars.distMapFig);
    end
    clf; 
    imagesc( nanmedian( mapsOut,3)); colorbar;
    hold on; 
    plot(pars.ctcf,pars.ctcf,'ko','MarkerSize',10);
    if isempty(pars.threshold)
        pars.threshold = ceil(mapStep*polyStep*2/3);
    end
    if pars.contMapFig == 0
        pars.f2 = figure(); 
    else
        pars.f2 = figure(pars.contMapFig);
    end
    clf; 
    imagesc(log2(ContactFrac( mapsOut,'threshold',pars.threshold))); 
    colorbar;
    hold on; 
    plot(pars.ctcf,pars.ctcf,'ko','MarkerSize',5);
end
