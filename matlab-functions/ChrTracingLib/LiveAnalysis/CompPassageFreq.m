function [passageFreq,passageFreqCI,passageFreqRep] =  CompPassageFreq(traj_list,varargin)
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
defaults(end+1,:) = {'iters', 'positive', 10};
defaults(end+1,:) = {'cI', 'positive', 100*(1 - (.05)^(1/2))};
defaults(end+1,:) = {'maxSamples', 'positive', 100};
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'skipLast', 'boolean', true}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);

if pars.d_contact == 0
    allTraj = cat(1,traj_list{:});
    pars.d_contact = quantile(allTraj(:),.01);
end
d_int = pars.d_contact;  % nm  (window size)
d_low = pars.d_contact;  % nm

nBoots = pars.iters; % could be parameter

T = size(traj_list{1},2); % total time points 
if pars.skipLast
    T=T-1;
end
nExp = length(traj_list);
passageFreqRep = nan(nExp,nBoots);
% figure(6); clf;
for b=1:nBoots 
    for m=1:nExp
        traj = traj_list{m}(:,1:T);
        d_high = pars.d_start*nanmedian(traj(:)); % 400; % nm
        nT = size(traj,1);
       
        resample = randi(nT,nT,1);
        traj = traj(resample,:); 
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
            passageFreqRep(m,b) = length(pt)/(sum(timeTotal));
            % figure(6); subplot(1,7,m); hist(pt,0:30:dMax); xlim([0,dMax]);
            % ylabel(nanmedian(pt));
            % title(dis_g(m))
            % pause(.1);
        catch
        end
    end
end
passageFreq = nanmedian(passageFreqRep,2);
passageFreqCI(:,1) = quantile(passageFreqRep,1-pars.cI/100,2);
passageFreqCI(:,2) = quantile(passageFreqRep,pars.cI/100,2);

