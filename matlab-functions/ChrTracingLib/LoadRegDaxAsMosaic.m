function [imageTiles,stageXY,pars] = LoadRegDaxAsMosaic(maxNameInHybFov,varargin)
% Inputs:
%   nHybs x nFOV x nDataChns cell array of max dax file names.  
%       - Each hyb may contain multiple data channels. The first channel is 
%         assumed to be fiduical data.  There may be multiple data channels
%         for each fiducial/hyb.  Currently all data channels will be
%         loaded.
% FIXME:  maxNames already specify whether they are 'dat' or 'fid' there is
% no need for assumptions.
% 
%
%
%
% Outputs:
%   nFOV x nChns cell array of image tiles
%   nFOV x 2 matrix of XY stage positions
%      - These outputs are ready to pass to MosaicViewerGUI
%      - Hint: if you loop over channels, you may pass these to
%      MosaicViewerRender as well to get mosaics without the GUI. 
% Optional Inputs
% 'backgroundCorrect' 'flatten
%    - currently does lower decile of each pixel -- this is an attempt to
%    find the pixels from the image-data that don't overlap cells/legit
%    signal but represent the actual image bacground.
%    - this background signal is straight *subtracted* from the data
%    images. Division would be better -- if illumination is 2x in the
%    middle vs. the edge, the difference in signal should be 2x the bright
%    stuff, not just 2x the noise. But division has problems with dividing
%    low numbers by low-numbers. Maybe we need a model S = N_a + N_m*D
% 

defaults = cell(0,3);
% default pars (for parsing file names)
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'transpose','boolean',true};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'fovs','freeType',inf};
defaults(end+1,:) = {'hybs','freeType',inf};
defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'data'};
% ---- Scope specific Pars
defaults(end+1,:) = {'scope',{'autoDetect','scope1','scope2','scope3','other'},'autoDetect'};  %
defaults(end+1,:) = {'transpose','boolean',true}; % see the mosaic parameters in the Hal parameter file used to collect the movie 
defaults(end+1,:) = {'fliplr','boolean',true};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'flipud','boolean',false};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'pix_to_mm','positive',6.55};  % 6.55 = scope 1. 6.45 = scope 2
%---- DaxToImageTiles Pars
defaults(end+1,:) = {'saveProjections','boolean',true};  % save the max projected images for quick loading 
defaults(end+1,:) = {'loadProjections','boolean',true};  % load previously save the max projected images for quick loading 
defaults(end+1,:) = {'trimBorder','nonnegative',0};
defaults(end+1,:) = {'selectFrame','integer',2}; % show only a single frame (overriden by projectAll)
defaults(end+1,:) = {'projectAll','boolean',true}; % show a projection of x,y data
%---flatten Background Pars
defaults(end+1,:) = {'smoothBox','positive',1};
defaults(end+1,:) = {'bkd','freeType',[]}; % map of illumination background
defaults(end+1,:) = {'camBkd','freeType','auto'}; 
defaults(end+1,:) = {'backgroundCorrect',{'file','median','medianEdge','none','removeData'},'median'};
pars  = ParseVariableArguments(varargin,defaults,'MosaicBuilder ValidateMosaic');


% try to guess saveFolder if none is passed
if isempty(pars.saveFolder)
   saveFolder = [fileparts(maxNameInHybFov{1}),filesep,'..',filesep,'Analysis',filesep];
else
   saveFolder = pars.saveFolder; 
end

[numHybs,numFovs,numDatas] = size(maxNameInHybFov);
% handle 'all hybs' request with inf
if isinf(pars.hybs)
    hybs=1:numHybs;
else
    hybs = pars.hybs;
end
% handle 'all fovs' request with inf
if isinf(pars.fovs)
    fovs = 1:numFovs;
else
    fovs = pars.fovs;
end
numHybs = min(numHybs,length(hybs));
numFovs = min(numFovs,length(fovs));
% load requested data chns
% Note: imTiles is nFOV x nChns for MosaicViewer GUI

% if LoadDaxFromEtable was asked for data only, it
if strcmp(pars.dataType,'data')
    chns = 1:numDatas; % 2:numDatas;
    % numDatas = numDatas - 1;
    imTileNames = maxNameInHybFov(hybs,fovs,chns);
    imTileNames = reshape(permute(imTileNames,[2,3,1]),numFovs,numHybs*numDatas);
    chns = 1:numHybs*(numDatas); % in the new imTiles strcture, these are the hybs
elseif strcmp(pars.dataType,'fiducial')
    chns = 1;
    numDatas = 1;
    imTileNames = maxNameInHybFov(hybs,fovs,1)';
elseif strcmp(pars.dataType,'all')
    fidTiles = maxNameInHybFov(1,fovs,1)';
    chns = 2:numDatas;
    numDatas = numDatas - 1;
    imTileNames = maxNameInHybFov(hybs,fovs,chns);
    imTileNames = reshape(permute(imTileNames,[2,3,1]),numFovs,numHybs*numDatas);
    imTileNames = cat(2,imTileNames,fidTiles);
    chns = 1:numHybs*(numDatas)+1; % in the new imTiles strcture, these are the hybs
end
regDatas = cell(numFovs,1);
imageTiles = imTileNames; 
% load per-hyb registration data for all FOVs 

for f=fovs
    try
        regTable = readtable([saveFolder,'fov',num2str(f,'%03d'),'_regData.csv']);
        regDatas{f} = table2struct(regTable);
    catch 
        if numHybs > 1
           cprintf([1 .5 0],'Warning: Multiple hybs were passed, but program did not find "fovNNN_regData.csv" to correct Hyb-to-hyb drift'); 
           cprintf([1 .5 0],'Data returned is not drift corrected.');
        end
    end
end
try
for h=chns
    disp(['loading tiles from chn ',num2str(h)]);
    if h==chns(1) || isempty(regDatas{1}) || h>length(regDatas{1}) % Load the first hyb
        [chnTiles,stageXY,pars] = DaxToImageTiles(imTileNames(:,h),'parameters',pars);
        % this will also load the scope specific parameters 
    else % all other hybs must also be drift corrected
        chnTiles = DaxToImageTiles(imTileNames(:,h),'parameters',pars);
        for f=fovs
            r = ceil((h-1)/numDatas); % 
            chnTiles{f} =  ApplyReg(chnTiles{f},regDatas{f}(r)); % register tile
            % just testing commenting this out to see if it creates the horrible border glare 
        end
    end
    % Background correction
    [chnTiles,bkd] = FlattenBackground(chnTiles,'parameters',pars);   
    % figure(10); clf; imagesc(bkd); colorbar; % for troubleshooting
    % save data
    imageTiles(:,h) = chnTiles;
end
catch er
    disp(er.getReport);
    disp('error here');
end
if pars.verbose
    disp('loading mosaic complete')
end