function [polyReps,mapReps,pars] = LoadPolySimDynamics(folder,varargin)
% Outputs:
% polyReps - cell array of replicates, each nMonomers x 3 x timesteps
% mapReps - cell array of replicates, each nMonomers x nMonomers x timesteps

defaults = cell(0,3); 
defaults(end+1,:) = {'polyStep','integer',5}; % downsample the polymer (in units of monomers)
defaults(end+1,:) = {'mapStep','integer',5}; % downsample the distance map (in units of polymer steps, applied after polySteps)
defaults(end+1,:) = {'timeStep','integer',1}; 
defaults(end+1,:) = {'TADs','freeType',[]}; % downsample the polymer (in units of monomers)
defaults(end+1,:) = {'blockSize','integer',50}; % should determine from file name
defaults(end+1,:) = {'threshold','integer',[]}; % monomer radii for contact map (empty for auto compute from polyStep mapStep)
defaults(end+1,:) = {'lenPoly','integer',0}; % should load from file
defaults(end+1,:) = {'rptPoly','integer',1}; % should load from file
defaults(end+1,:) = {'showPlots','boolean',false}; % number of examples per block to load
defaults(end+1,:) = {'distMapFig','integer',0}; % Figure to show distance map in
defaults(end+1,:) = {'contMapFig','integer',0}; % Figure to show distance map in
defaults(end+1,:) = {'maxBlocks','integer',inf};
defaults(end+1,:) = {'readPars','boolean',false};
defaults(end+1,:) = {'offset','integer',0}; % offset polymer index by this number of monomers (avoids getting only integer ctcf sites)
pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.lenPoly == 0
    simPars = [folder,'simPars.txt'];
    if exist(simPars,'file')
        parsTable = readtable(simPars);
        pars.lenPoly = parsTable.monomers;
        pars.rptPoly = parsTable.replicates;
        pars.TADs = parsTable.ctcfSites;
    else
        pars.lenPoly = 4e3; % original default size
    end
end

blockFiles = FindFiles([folder,'*blocks_*.h5']);
% determine block size from name
if isempty(blockFiles)
    error(['no data blocks found in ',folder]);
end
[~,blockName] = fileparts(blockFiles{1});
blockVal = strsplit(blockName,'_');
blockSize = -(eval(blockVal{2})-1); % turn 0-9 into 10

% blockFiles = cellstr(ls([folder,'*blocks_*.h5']));

polyStep = pars.polyStep;
mapStep = pars.mapStep;

% 0..9900
nB = min([pars.maxBlocks,length(blockFiles)]); % 99; % number of blocks.  for speed lets not use all of them;
if length(blockFiles) > pars.maxBlocks
    disp(['detected ', num2str(length(blockFiles)) ' maxBlocks ',num2str(pars.maxBlocks)]);
end
% exPerBlock = pars.exPerBlock; % 2;
% blockSize =pars.blockSize; % 100;
rptPoly = pars.rptPoly; % 10;
timeStep = min(pars.timeStep,blockSize); % can't be larger than block size
totTime = ceil(nB*blockSize/timeStep);
pars.N = pars.lenPoly*pars.rptPoly; % 40e3; % total size of polymer in simulation
lenPolyOrig = pars.lenPoly;
lenPolyDs = lenPolyOrig/polyStep;
subset = (1+pars.offset):polyStep:pars.N;
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
            m1 = round((p-1)*lenPolyDs+1);  % 
            m2 = round(p*lenPolyDs);         
            temp = double(polyData(m1:m2,:));
            polyReps{p}(:,:,tt) = temp;
        end
    end
end


if nargout > 1
    mapSteps = 1:mapStep:stps;
    nBins = length(mapSteps);
    mapReps = cell(rptPoly,1);
    for p=1:rptPoly
        mapReps{p} = nan(nBins,nBins,totTime);
        for t=1:totTime
           mapReps{p}(:,:,t) =  squareform(pdist(polyReps{p}(mapSteps,:,t)));
        end
    end
    mapsOut = cat(3,mapReps{:});
end

if pars.showPlots && nargout > 1
    pars.ctcf = ceil(pars.TADs/mapStep/polyStep);
    if pars.distMapFig == 0
        pars.f1 = figure(); 
    else
        pars.f1 = figure(pars.distMapFig);
    end
    clf; imagesc( nanmedian( mapsOut,3)); colorbar;
    hold on; plot(pars.ctcf,pars.ctcf,'ko','MarkerSize',10);
    if isempty(pars.threshold)
        pars.threshold = ceil(mapStep*polyStep*2/3);
    end
    if pars.contMapFig == 0
        pars.f2 = figure(); 
    else
        pars.f2 = figure(pars.contMapFig);
    end
    clf; imagesc( log2(ContactFrac( mapsOut,'threshold',pars.threshold))); colorbar;
    hold on; plot(pars.ctcf,pars.ctcf,'ko','MarkerSize',5);
end