function [passageTime_median,passageTime_CI,passageTimes,passageTimes_censored] =  CompContactWaitingKM(traj_list,varargin)
% Compute First Passage Time from median distance to contact
%  Modified to report censored trajectories 
%      - front-censored: reached contact, but didn't observe the start at
%      average separation. 
%      - end-censored: left-average separation and got closer, but didn't
%      return to average or to proximity before the end of the window. 
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
% Updated to account for censoring 
% larger values of d_contact slow down computation considerably. 


defaults = cell(0,3); % 
% shared parameters for moving average
defaults(end+1,:) = {'d_contact','positive',0}; % separation distance considered contact. also used for window size at d_start
defaults(end+1,:) = {'d_start','positive',1}; % multiplier of the median distance. used for start separation;  
defaults(end+1,:) = {'iters', 'positive', 20};
defaults(end+1,:) = {'cI', 'positive', 100*(1 - (.05)^(1/2))};
defaults(end+1,:) = {'boot_CI','boolean',true};
defaults(end+1,:) = {'maxSamples', 'positive', 1e4};
defaults(end+1,:) = {'method', {'median','mean'},'median'};
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'skipLast', 'boolean', true};  % if last column records cell-cycle state, we don't want to confuse it
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);

if pars.d_contact == 0
    allTraj = cat(1,traj_list{:});
    pars.d_contact = quantile(allTraj(:),.01);
end
d_int = pars.d_contact;  % nm  (window size)
d_low = pars.d_contact;  % nm


try

T = size(traj_list{1},2); % total time points 
if pars.skipLast
    T=T-1;
end
nExp = length(traj_list);
passageTime_median = nan(nExp,3); % 3 different medians:  1 = orig, full-passage only; 2= raw median (both censored and uncensored); 3= KM corrected median
passageTime_CI = nan(nExp,2,3); % 3 different medians:  1 = orig, full-passage only; 2= raw median (both censored and uncensored); 3= KM corrected median
passageTimes = cell(nExp,5);
passageTimes_censored = cell(nExp,1);
passageTimesPerTrace = cell(nExp,1);
% figure(6); clf;
for m=1:nExp
    traj = traj_list{m}(:,1:T);
    d_high = pars.d_start*nanmedian(traj(:)); % 400; % nm
    nT = size(traj,1);
    passageTimesPerTrace{m} = cell(nT,1);
    censoredPassages  = nan(nT,3);
    timeTotal = T*ones(nT,1);
    % compute for all times in the ave window, how long it takes to get to proximity.  
    for t=1:nT  % t=2 (loop over traces)
        isA = traj(t,:) > d_high & traj(t,:) < d_high + d_int;
        isB = traj(t,:) > 0 & traj(t,:) < d_low;
        % figure(3); clf; imagesc([isA; isB]);
        as = find(isA);
        bs = find(isB);
        
        % handle censoring
        %   each trace has most one censored passage at the start and one
        %   at the end
        %    traces that stay between the average and the contact are also censored  
        if ~(isempty(bs) && isempty(as))  % if both have no data, just keep going        
            % first deal with the missing data cases
            if isempty(as) % if no start points
                censoredPassages(t,1) = bs(1);  % save first contact point
            end
            if isempty(bs) % no end points 
                censoredPassages(t,2) = T-as(end); % handle end-censoring
            end
            if ~isempty(as) && ~isempty(bs) % neither is empty    
                if bs(1) < as(1)
                    censoredPassages(t,1) = bs(1); 
                elseif (as(end) > bs(end)) && (as(end) > T )
                    censoredPassages(t,2) = T-as(end);  % handle end-censoring
                end
            end
        else
            fully_inbetween = sum( (traj(t,:)>d_high) )==0 &&  sum((traj(t,:)<d_low))==0;
            if fully_inbetween
                censoredPassages(t,3) = T;
            end
        end

        tt = 0;  % counter over time
        k=0;
        gaps = zeros(length(as),1);
        while tt<T && ~isempty(as)  && ~isempty(bs) && k<pars.maxSamples
            try
            tt = as(1); % start at the first distal state
            bs(bs<tt) = []; % contact prior to the distal state don't count
            if ~isempty(bs)
                bt = bs(1); % the first waiting time is just the first contact event after the first distal state

                curr_distal =   as(as<=bs(1));
                [~,idmax] = max(curr_distal);  % start at the first time we leave the distal state 
                at = curr_distal(randi(idmax));  % chose a random event, among the times when it was at start distance, to avoid resampling the same passage.
                k=k+1;
                gaps(k) = bt-at;
                as(as<bt) = []; % remove all distal times before the current contact time     
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
        censored_passages = censoredPassages(:);
        censored_passages(isnan(censored_passages)) =[];
        passageTimes{m,1} = cat(1,passageTimesPerTrace{m}{:});
        passageTimes_censored{m} = censored_passages; 
        censored = [zeros(size(passageTimes{m,1})); ones(size(passageTimes_censored{m}))];
        pt_all = [passageTimes{m,1}; passageTimes_censored{m}];
        if isempty(pt_all)
            continue
        end
        censored(pt_all<=1) = -1; % any searches of 1 frame are left censored
        if sum(censored) ~= size(pt_all)
            boots = randi(length(pt_all),length(pt_all),pars.iters);
            meds = nan(pars.iters,1);
            for b=1:pars.iters
                [f,x] = ecdf(pt_all(boots(:,b)),'Censoring',censored(boots(:,b)));  
                med_KM_idx = find(f>=0.5,1);
                if isempty(med_KM_idx)
                    med_KM = max(x);
                else
                    med_KM = x(med_KM_idx);
                end
                if ~isempty(med_KM)
                    meds(b) = med_KM;
                else
                    meds(b) = nan;
                end
            end
            [f,x,flo,fup] = ecdf(pt_all,'Censoring',censored);  
            med_KM_idx = find(f>=0.5,1);
            if isempty(med_KM_idx)
                med_KM = max(x);
            else
                med_KM = x(med_KM_idx);
            end
            med_KM_lo = quantile(meds,1-pars.cI/100);
            med_KM_up = quantile(meds,pars.cI/100);

            passageTimes{m,2} = x;
            passageTimes{m,3} = f;
            passageTimes{m,4} = flo;
            passageTimes{m,5} = fup;
            % handle truncation
            % pdf = diff(f); 
            % med_KM =  pdf'*x(1:end-1); % average 
            % med_KM_lo=  diff(flo)'*x(1:end-1); % average 
            % med_KM_up=  diff(fup)'*x(1:end-1); % average 

            
            % med_lo_idx = find(flo>=0.5,1);
            % med_up_idx = find(fup>=0.5,1);     
            % if isempty(med_lo_idx)
            %     med_KM_lo = max(x);
            % else
            %     med_KM_lo = x(med_lo_idx);
            % end  
            % if isempty(med_up_idx)
            %     med_KM_up = max(x);
            % else
            %     med_KM_up = x(med_up_idx);
            % end  
        else
            if pars.boot_CI
                [med_KM,med_KM_ci] =  MedianWithCI(pt_all,'iters',pars.iters);
            else
                med_KM = nanmedian(pt_all);
                med_KM_ci = [nan,nan];
            end
        end
        passageTime_median(m,1) = nanmedian(passageTimes{m}); 
        passageTime_median(m,2) = nanmedian(pt_all); 
        passageTime_median(m,3) = med_KM; 

        if pars.boot_CI % slow
            [~,ci1] =  MedianWithCI(passageTimes{m},'iters',pars.iters);
            [~,ci2] =  MedianWithCI(pt_all,'iters',pars.iters);
        else
            ci1 = nan;
            ci2 = nan;
        end

        passageTime_CI(m,1:2,1) = ci1;
        passageTime_CI(m,1:2,2) = ci2;
        passageTime_CI(m,1:2,3) = [med_KM_lo,med_KM_up];
    catch  er
        warning(er.getReport);
        disp('debug here');
    end
end


catch er
    warning(er.getReport);
    disp('debug here');
end
