function [passageTime_median,passageTime_CI,passageTimesMerged,passageFreq,passageTimesPerTrace] =  CompFirstPassage(traj_list,varargin)
% Compute First Passage Time from median distance to contact
%  
%% Inputs
% traj_list - cell array of n-experiments, each containing a traj_matrix
%      traj_matrix = C x T array of distances between two points in C-cells
%      observed at T timepoints.  Pad with NaNs for missing data.
% 
%% Optional Inputs
% d_contact - threshold for 'proximity'/contact. default = 0
% d_start - multiplier of the median distance, used for starting separation
%  
% 
%% Outputs
% passageTime_median,
% passageTime_CI - confidence interval (95%) 
% passageTimes - cell array, full distribution of passage times
% 
%% Notes
%
% larger values of d_contact slow down computation considerably. 


defaults = cell(0,3); % 
% shared parameters for moving average
defaults(end+1,:) = {'d_contact','positive',0}; % separation distance considered contact. also used for window size at d_start
defaults(end+1,:) = {'d_start','positive',1}; % multiplier of the median distance. used for start separation;  
defaults(end+1,:) = {'iters', 'positive', 1000};
defaults(end+1,:) = {'cI', 'positive', 100*(1 - (.05)^(1/2))};
defaults(end+1,:) = {'maxSamples', 'positive', 100};
defaults(end+1,:) = {'method', {'median','mean'},'median'};
defaults(end+1,:) = {'verbose', 'boolean', false};
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);

if pars.d_contact == 0
    allTraj = cat(1,traj_list{:});
    pars.d_contact = quantile(allTraj(:),.01);
end
d_int = pars.d_contact;  % nm  (window size)
d_low = pars.d_contact;  % nm



T = size(traj_list{1},2); % total time points 
nExp = length(traj_list);
passageTime_median = nan(nExp,1);
passageTime_CI = nan(nExp,2);
pt_gm = nan(nExp,1);            % mean  (not currently returned)
pt_gse = nan(nExp,1);           % standard error (not currently returned)
passageTimesPerTrace = cell(nExp,1);
passageTimesMerged = cell(nExp,1);
passageFreq = nan(nExp,1);
% figure(6); clf;
for m=1:nExp
    traj = traj_list{m};
    d_high = pars.d_start*nanmedian(traj(:)); % 400; % nm
    nT = size(traj,1);
    passageTimesPerTrace{m} = cell(nT,1);
    timeTotal = T*ones(nT,1);
    for t=1:nT  % t=2 (loop over traces)
        isA = traj(t,:) > d_high & traj(t,:) < d_high + d_int;
        isB = traj(t,:) > 0 & traj(t,:) < d_low;
        % figure(3); clf; imagesc([isA; isB]);
        as = find(isA);
        bs = find(isB);
        tt = 0;  % counter over time
        k=0;
        gaps = zeros(length(as),1);
        while tt<T && ~isempty(as)  && ~isempty(bs) && k<pars.maxSamples
            try
            tt = as(1); %
            bs(bs<tt) = [];
            if ~isempty(bs)
                bt = bs(1);
                at = max(as(as<=bs(1)));
                k=k+1;
                gaps(k) = bt-at;
                as(as<bt) = [];
            end
            catch er
                warning(er.getReport);
                disp('debug here');
            end
        end
        if k==pars.maxSamples 
            timeTotal(t) = tt;
            if pars.verbose
            disp(['m=',num2str(m),' ',num2str(pars.maxSamples),' samples reached. stopping'])
            end
        end
        
        gaps(gaps==0) = [];
        if tt == T || isempty(gaps)
            % gaps = T;  % 
            if pars.verbose
                disp('no passage events')
            end
        end
        if ~isempty(gaps)
            passageTimesPerTrace{m}{t} = gaps;
        end
    end


    try
        pt = cat(1,passageTimesPerTrace{m}{:});
        passageTime_median(m) = nanmedian(pt);
        if strcmp(pars.method,'median')
        [passageTime_median(m), passageTime_CI(m,:)] = MedianWithCI(pt,'cI',pars.cI,'iters',pars.iters);
        elseif strcmp(pars.method,'mean')
            passageTime_median(m) = nanmean(pt);
            passageTime_CI(m,1) = nanmean(pt) - nanstd(pt)./sum(~isnan(pt(:))); %#ok<*NANSTD>
            passageTime_CI(m,2) = nanmean(pt) + nanstd(pt)./sum(~isnan(pt(:)));
        end
        pt_gm(m) = nanmean(pt);
        pt_gse(m) = nanstd(pt)/sqrt(-1+sum(~isnan(pt)));
        passageTimesMerged{m} = pt;
        passageFreq(m) = length(pt)/(sum(timeTotal));
        % figure(6); subplot(1,7,m); hist(pt,0:30:dMax); xlim([0,dMax]);
        % ylabel(nanmedian(pt));
        % title(dis_g(m))
        % pause(.1);
    catch
    end
end
