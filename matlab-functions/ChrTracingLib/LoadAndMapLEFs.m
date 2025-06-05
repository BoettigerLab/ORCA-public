function [extMapLR,lefPos] = LoadAndMapLEFs(lefFolder,varargin)
% Load LEFs
% create LEF maps
% 
% Outputs
% extMapLR T x N x reps x 2  map of LEF positions
%  with the "stackReps" flag "true", extMapLR is T*reps x N x 2
% 
% Example applications: 
% Show an image of the LEF maps from replicate r
%   extMapLR = LoadAndMapLEFs(lefFolder)
%   figure(1); clf; imagesc(extMapLR(:,:,r,1)+ 2*(extMapLR(:,:,r,2));
% 
% Make a "ChIPSeq" like plot
%     extMapLR = LoadAndMapLEFs(lefFolder,'stackReps',true)
%     figure(2); clf; plot(sum(sum(extMapLR,3),1));
% 
% Fraction of time at TAD boundaries
%     extMapLR = LoadAndMapLEFs(lefFolder,'stackReps',true)
%     extMap = sum(extMap,3); 
%     ctcf = extMap(:,tads+1);
%     other = extMap;
%     other(:,tads+1) = 0;
%     sum(ctcf(:))./sum(extMap(:))
%     sum(other(:))./sum(extMap(:))

defaults = cell(0,3); 
defaults(end+1,:) = {'parsFile','string',''}; % downsample the polymer step
defaults(end+1,:) = {'lenPoly','integer',[]}; % should load from file
defaults(end+1,:) = {'rptPoly','integer',1}; % should load from file
defaults(end+1,:) = {'stackReps','boolean',false}; % should load from file
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);

% lefFolder = 'T:\2023-05-05_Pitx1Sims\Model_v4_mid\'; % input 

if strcmp(lefFolder(end),filesep)
    lefFile = [lefFolder,'LEFPos.h5'];
else 
    lefFile = lefFolder;
end

if isempty(pars.parsFile)
    pars.parsFile = FindFiles([lefFolder,'simPars.txt']);
end
if ~isempty(pars.parsFile)
    parsTable = readtable(pars.parsFile{1});
    pars.rptPoly = parsTable.replicates;
    pars.TADs = str2num(parsTable.ctcfSites{1});
    pars.lenPoly = parsTable.monomers;
end

info = h5info(lefFile);
lefPos = h5read(lefFile,'/positions'); % 2 x N-extruders x T-timesteps
[~,nExtruders,tSteps] = size(lefPos);

if ~isempty(pars.lenPoly)
nMono = pars.lenPoly;
else
    nMono = max(lefPos(:));
end
nReps = pars.rptPoly;
extMapLR = zeros(tSteps,nMono,nReps,2,'logical');
for t=1:tSteps
    mono = rem(lefPos(1,:,t),nMono);
    rep = floor(lefPos(1,:,t)/nMono);
    for b=1:max(rep)
        isB = rep==b;
        extMapLR(t,mono(isB)+1,b,1) = true;
    end
    mono = rem(lefPos(2,:,t),nMono);
    rep = floor(lefPos(2,:,t)/nMono);
    for b=1:max(rep)
        isB = rep==b;
        extMapLR(t,mono(isB)+1,b,2) = true;
    end
end
if pars.stackReps
    extMapLeftT = reshape(permute(extMapLR(:,:,:,1),[1,3,2]),tSteps*nReps,nMono);
    extMapRightT = reshape(permute(extMapLR(:,:,:,2),[1,3,2]),tSteps*nReps,nMono);
    extMapLR = cat(3,extMapLeftT,extMapRightT);
end

% tads = [200,400,600,800]
% figure(1); clf; imagesc(extMapLeftT+ 2*(extMapRightT));
% figure(2); clf; plot(sum(extMapLeftT));
% % tads = pars.TADs;
% ctcf = extMapLeftT(:,tads+1);
% other = extMapLeftT;
% other(:,tads+1) = 0;
% sum(ctcf(:))./sum(extMapLeftT(:))


% numExtruders = info.Attributes(2).Value; % total LEFs across all replicates
% 
% 
% a = zeros(10)
% mono = rem(15:22,10)
% rep = floor( (15:22)/10)
% for b=1:max(rep)
%         isB = rep==b;
%         a(mono(isB)+1,b) = true;
% end
% a



