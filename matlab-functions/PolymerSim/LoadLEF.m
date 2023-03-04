function [imTraj2,pars,lefPos] = LoadLEF(lefFile,varargin)
% load lef file
defaults = cell(0,3); 
defaults(end+1,:) = {'parsFile','string',''}; % downsample the polymer step
defaults(end+1,:) = {'polyStep','integer',5}; % downsample the polymer step
defaults(end+1,:) = {'timeStep','integer',50}; % downsample the time step
defaults(end+1,:) = {'TADs','freeType',[300,800,1500,2300,2900,3400]}; % downsample the polymer (in units of monomers)
defaults(end+1,:) = {'lenPoly','integer',4e3}; % should load from file
defaults(end+1,:) = {'rptPoly','integer',10}; % should load from file
defaults(end+1,:) = {'showInfo','boolean',true}; % number of examples per block to load
defaults(end+1,:) = {'trajFig','integer',2}; % number of examples per block to load
defaults(end+1,:) = {'avePosFig','integer',3}; % number of examples per block to load
pars = ParseVariableArguments(varargin,defaults,mfilename);

% pars.parsFile='Z:\Alistair\2022-02-04_PolychromWholeNuc\NT_rad21_v1\rep028\simPars.txt'     
    
%% LEF position analysis
if ~isempty(pars.parsFile)
    parsTable = readtable(pars.parsFile);
    pars.rptPoly = parsTable.replicates;
    pars.TADs = str2num(parsTable.ctcfSites{1});
    pars.lenPoly = parsTable.monomers;
end

% --- LEF positions
% lefFile =[folder,'LEFPositions','.h5']
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
lenPoly = pars.lenPoly;
rptPoly = pars.rptPoly;
TADs = pars.TADs/pars.polyStep;
nMonomers = pars.lenPoly*pars.rptPoly;
tSteps = round(size(lefPos,3)/pars.timeStep);
nMon = round(nMonomers/pars.polyStep);
imLeft = false(tSteps,nMon);
imRight = false(tSteps,nMon);
extruders = 1:numExtruders; % which of the extruders to plot
for t=1:tSteps
    leftSteps = squeeze(round(lefPos(1,extruders,t)/pars.polyStep));
    leftSteps(leftSteps==0) = nMon; % just for plotting
    imLeft(t,leftSteps) = true; 
    rightSteps = squeeze(round(lefPos(2,extruders,t)/pars.polyStep));
    rightSteps(rightSteps==0) = nMon; 
    imRight(t,rightSteps) = true; 
end
imTraj = imLeft + 2*imRight; % colors L and R extruders differently.

% % for troubleshooting
% figure(1); clf; imagesc(imTraj); hold on;
% tadClr = GetColorMap('hsvCut',length(TADs));
% for k=1:rptPoly
%     x = (k-1)*lenPoly/pars.polyStep+repmat(TADs',1,2); 
%     y = repmat([0,tSteps],length(TADs),1);
%     plot(x',y','r--');
%     % color code TADs
%     for c=1:length(TADs)
%         x = (k-1)*lenPoly/pars.polyStep+repmat(TADs(c)',1,2); 
%         y = [0,tSteps];
%         plot(x',y','.-','color',tadClr(c,:));
%     end
% end
% figure(2); clf; subplot(1,2,1); plot(sum(imTraj,1));
% subplot(1,2,2); plot(sum(imTraj,2));


% fold matrix back on itself based on rpt Poly
imTraj2 = cell(rptPoly,1);
for k=1:rptPoly
    m1 = (k-1)*lenPoly/pars.polyStep+1;  % 
    m2 = k*lenPoly/pars.polyStep;
    imTraj2{k} = imTraj(:,m1:m2);
end
imTraj2 = cat(1,imTraj2{:});

if pars.trajFig ~= 0
    uni = ones(length(TADs),1);
    t = [0*uni,tSteps*rptPoly*uni,nan*uni]';t= t(:);
    xt = [TADs; TADs; TADs]; xt=xt(:);
    
%     figure(pars.trajFig); clf; 
%     imagesc(imTraj2); colorbar;
%     hold on; plot(xt,t,'r--');
%     ylim([0,max(t)]);
    % line plot version
    [y,x] = find(imTraj2==1);
    figure(3); clf; plot(x,y,'r.','MarkerSize',1);
    [y,x] = find(imTraj2==2);
    figure(3); hold on; plot(x,y,'b.','MarkerSize',1);
    hold on; plot(xt,t,'k.-');
    
    
end

if pars.trajFig ~= 0
    figure(pars.trajFig); clf;
    imagesc(imTraj2);
end

if pars.avePosFig ~= 0
    figure(pars.avePosFig); clf;
    plot(sum( imTraj2>0) ,'k');
end