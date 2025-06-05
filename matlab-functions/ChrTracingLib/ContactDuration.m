function [contactAve,contactSEM,contactAll] = ContactDuration(cellTraj,varargin)


defaults = cell(0,3); % 
% shared parameters for moving average
defaults(end+1,:) = {'theta','positive',50}; % separation distance considered contact. also used for window size at d_start
defaults(end+1,:) = {'minConsecutive','positive',1}; %
defaults(end+1,:) = {'skipLast', 'boolean', true};  % if last column records cell-cycle state, we don't want to confuse it
pars = ParseVariableArguments(varargin,defaults,mfilename);


% cellTraj - cell array, nS data types (typically 2, +/- cohesin) x nG
% genomic gaps.

[nS,nG] = size(cellTraj);
[~,T] = size(cellTraj{1,1});
if pars.skipLast
    T=T-1;
end
contactDuration = cell(nG,nS);
contactAll = cell(nG,nS);
contactAve = nan(nG,nS);
contactStd = nan(nG,nS);
contactSEM  = nan(2,nG,nS);
nObs = nan(nG,nS);
for d=1:nS  % loop over data types
    for e=1:nG  % loop over experiments (distance intervals)
        disM = cellTraj{d,e}(:,1:T);
        nC = size(disM,1);  % (always 10 here)
        contactDuration{e,d} = cell(nC,1);
        for c=1:nC % loop over traces  / cells; 
            inContact = disM(c,:)<pars.theta;
            noData = isnan(disM(c,:));
            noContact = disM(c,:)>pars.theta;
            inContact(~noData);
            inCont = bwareaopen(inContact+noData >0, pars.minConsecutive); % at least minConsecutive frames of contact
            noCont = bwareaopen(noContact+noData >0, pars.minConsecutive); % at least  minConsecutive frames of no contact
            inCont(noCont) = 0;
            regs = regionprops(inCont,inContact,'PixelIdxList','Area','PixelValues');
            regs(cellfun(@sum,{regs.PixelValues}) < pars.minConsecutive/2) = []; % toss blocks that don't have at least half data points  
            contactIntervals = {regs.PixelIdxList}; % the timepoints of all the contact intervals
            contactDuration{e,d}{c} = ([regs.Area]);
        end
        contactAll{e,d} = cat(2,contactDuration{e,d}{:});
        contactAve(e,d) = nanmean(cat(2,contactDuration{e,d}{:}));
        contactStd(e,d) = nanstd(cat(2,contactDuration{e,d}{:})); %#ok<*NANSTD>
        nObs(e,d) =  sum(~isnan((cat(2,contactDuration{e,d}{:}))));
    end
     contactSEM(:,:,d) =  [(contactAve(:,d)-contactStd(:,d)./sqrt(nObs(:,d)))'; (contactAve(:,d)+contactStd(:,d)./sqrt(nObs(:,d)))'];  
    
    %       figure(5); 
    % semilogx(g_dists,contactAve(:,d),'.','color',cmapRB((cc-1)*2+d,:));   
    % hold on; plot([g_dists; g_dists],er,'-','color',cmapRB((cc-1)*2+ d,:));
end