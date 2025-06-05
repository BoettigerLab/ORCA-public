function typicalPassage = ExamplePassages(traj_list,varargin)
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
defaults(end+1,:) = {'maxStep', 'positive', 40}; % max average step for examples 
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
typicalPassage = cell(nExp,1);
passageTimesPerTrace= cell(nExp,1);
passageExamples= cell(nExp,1);

for m=1:nExp
    traj = traj_list{m};
    d_high = pars.d_start*nanmedian(traj(:)); % 400; % nm
    nT = size(traj,1);
    passageTimesPerTrace{m} = cell(nT,1);
    passageExamples{m} = cell(nT,1);
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
        passageTrajs = cell(length(as),1);
        while tt<T && ~isempty(as)  && ~isempty(bs) && k<pars.maxSamples
            try
            tt = as(1); %
            bs(bs<tt) = [];
            if ~isempty(bs)
                bt = bs(1);
                at = max(as(as<=bs(1)));
                k=k+1;
                try
                    passageTrajs{k} = traj(t,at-10:bt+10);
                    gaps(k) = bt-at;
                catch
                end
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
        passageTrajs(gaps==0) = [];
        gaps(gaps==0) = [];
        if tt == T || isempty(gaps)
            % gaps = T;  % 
            if pars.verbose
                disp('no passage events')
            end
        end
        if ~isempty(gaps)
            passageTimesPerTrace{m}{t} = gaps;
            passageExamples{m}{t} =passageTrajs;
        end
    end


    try
        pt = cat(1,passageTimesPerTrace{m}{:});
        N = length(pt);
        pes = cell(N,1);
        N1 = length(passageExamples{m});
        k=0;
        for n=1:N1
            N2 = length(passageExamples{m}{n});
            for nn=1:N2
                k=k+1;
                pes{k} = passageExamples{m}{n}{nn};
            end
        end
        trace_length = cellfun(@length,pes)-20;
        dist_med = abs(trace_length - nanmedian(trace_length))/nanmedian(trace_length);
        step_mean = cellfun(@(x) nanmean(abs(diff(x))),pes);
        step_max  = cellfun(@(x) quantile(abs(diff(x(~isnan(x)))),.95),pes);
        frac_observed = cellfun(@(x) sum(~isnan(x)),pes)./trace_length;

        if length(trace_length) >= 2
            keep = frac_observed > 0.75 & step_mean < 50 & step_max < 250;
            typPes = pes(keep); 
            dist_keep = dist_med(keep);
            [dk,idx] = sort(dist_keep,'ascend');
            typicalPassage{m} = typPes(idx);
        else
            typicalPassage{m} =pes;
        end
    catch
    end
end



%%


% 
% 
% figure(10); clf; 
% ys = typicalPassage{m}{2};
% y1 = ys(11:end-10);
% ts = 1:length(ys);
% ts(isnan(ys)) = [];
% ys(isnan(ys)) = [];
% t1 = 10+(1:length(y1));
% t1(isnan(y1)) = [];
% y1(isnan(y1)) = [];
% plot(ts,ys,'.-'); hold on;
% plot(t1,y1,'.-'); box off;






% figure(10); clf; plot(passageExamples{1}{1}
% passLength  = cellfun(@length,passageExamples{1})
% passLength(passLength==0) = [];
% nanmedian(passLength)



