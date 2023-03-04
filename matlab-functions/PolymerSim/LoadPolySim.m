function [polysOut,mapsOut,pars] = LoadPolySim(folder,varargin)

% folder ='T:\2020-06-21_PolychromSims\Simulation1\trajectory\'
% folder ='T:\2020-06-21_PolychromSims\Simulation2\trajectory\' % stiff=0, full LE
% notes on simulation speed - ~1 min per block (20 min per block)
% 24 hours for 999 blocks  


defaults = cell(0,3); 
defaults(end+1,:) = {'polyStep','integer',5}; % downsample the polymer (in units of monomers)
defaults(end+1,:) = {'mapStep','integer',5}; % downsample the distance map (in units of polymer steps, applied after polySteps)
defaults(end+1,:) = {'TADs','freeType',[]}; %  (in units of monomers)
defaults(end+1,:) = {'blockSize','integer',50}; % new default in simulation is block size 50
defaults(end+1,:) = {'exPerBlock','integer',0}; % by default let's load all the blocks there are. We'll make this default to blocksize unless otherwise specified
defaults(end+1,:) = {'threshold','integer',[]};
defaults(end+1,:) = {'lenPoly','integer',1e3}; % should load from file
defaults(end+1,:) = {'rptPoly','integer',10}; % should load from file
defaults(end+1,:) = {'showPlots','boolean',true}; % number of examples per block to load
defaults(end+1,:) = {'distMapFig','integer',0}; % Figure to show distance map in
defaults(end+1,:) = {'contMapFig','integer',0}; % Figure to show distance map in
defaults(end+1,:) = {'maxBlocks','integer',inf};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.exPerBlock == 0
    pars.exPerBlock = pars.blockSize;
end

blockFiles = cellstr(ls([folder,'*blocks_*.h5']));

polyStep = pars.polyStep;
mapStep = pars.mapStep;

% 0..9900
nB = min([pars.maxBlocks,length(blockFiles)]); % 99; % number of blocks.  for speed lets not use all of them;
if length(blockFiles) > pars.maxBlocks
    disp(['detected ', num2str(length(blockFiles)) ' maxBlocks ',num2str(pars.maxBlocks)]);
end
exPerBlock = pars.exPerBlock; % 2;
blockSize =pars.blockSize; % 100;
rptPoly = pars.rptPoly; % 10;
totCells = nB*exPerBlock*rptPoly;
pars.N = pars.lenPoly*pars.rptPoly; % 40e3; % total size of polymer in simulation
lenPolyOrig = pars.lenPoly;
lenPolyDs = lenPolyOrig/polyStep;
subset = 1:polyStep:pars.N;
stps = ceil(length(subset)/rptPoly);
polysOut = zeros(stps,3,totCells);
c=0;
for b=1:nB
    blockStart = num2str(  blockSize*(0+b-1) );
    blockEnd   = num2str(  blockSize*(0+b-1) + (blockSize-1));
    file = [folder,'blocks_',blockStart,'-',blockEnd,'.h5'];
    for e=1:exPerBlock
        ex = num2str(  blockSize*(0+b-1) + (e-1)*(blockSize/exPerBlock) );
        data = h5read(file,['/',ex,'/pos']); % there should be an easy way to load subsets of the data with h5read, just need to look into it
        polyData = data(:,subset)';
        for p=1:rptPoly
            m1 = round((p-1)*lenPolyDs+1);  % 
            m2 = round(p*lenPolyDs);
            c=c+1;           
            temp = double(polyData(m1:m2,:));
            polysOut(:,:,c) = temp;
        end
    end
end

if nargout > 1 | pars.showPlots
    mapSteps = 1:mapStep:stps;
    nBins = length(mapSteps);
    mapsOut = zeros(nBins,nBins,totCells);  % don't make maps at full res!
    for c=1:totCells
        mapsOut(:,:,c) = squareform(pdist(polysOut(mapSteps,:,c)));
    end
end

if pars.showPlots
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