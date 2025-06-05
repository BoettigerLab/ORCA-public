function [passageAligned,traceNum,traceStart,traceLength] = AlignPassage(traj,varargin)



defaults = cell(0,3); % 
% shared parameters for moving average
defaults(end+1,:) = {'d_contact','positive',1}; % separation distance considered contact. also used for window size at d_start
defaults(end+1,:) = {'d_int','positive',0}; % separation distance considered contact. also used for window size at d_start
defaults(end+1,:) = {'d_start','positive',1}; % Distance used for start separation;  
defaults(end+1,:) = {'minPoints','positive',30}; % min observations in a passage event
defaults(end+1,:) = {'maxPoints','positive',5e3}; % max points in trace to record
defaults(end+1,:) = {'maxSamples','positive',6e3}; % max number of samples to analyze
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);

if pars.d_int == 0 
    pars.d_int = pars.d_contact;
end

d_low = pars.d_contact;
d_high = pars.d_start; % 30; % nanmedian(traj(:)) % 30; % 
d_int = pars.d_int; % d_low; 
nP = pars.maxPoints; % 5000; % round(40+dis_kb1(e)^(2/3)); % max number of points in passage
[nT,T] = size(traj);

% numNans = sum(isnan(traj),2);
%       [~,idx] = sort(numNans);
%       traj = traj(idx,:);



ptPerTrace = cell(nT,1);
passageAligned = nan(pars.maxSamples,nP+1);
traceNum = nan(pars.maxSamples,1);
traceStart = nan(pars.maxSamples,1);
c=0;
for t=1:nT  % t=2 (loop over traces)  t=1
    % d_high = nanmedian(traj(t,:));
    isA = traj(t,:) > d_high & traj(t,:) < d_high + d_int;
    isB = traj(t,:) > 0 & traj(t,:) < d_low;
    % figure(3); clf; imagesc([isA; isB]);
    as = find(isA);  % high
    bs = find(isB); % low
    tt = 0;  % counter over time
    k=0;
    gaps = zeros(length(as),1);
    while tt<T && ~isempty(as)  && ~isempty(bs) && c<pars.maxSamples
        % try
        tt = as(1); %
        bs(bs<tt) = [];
        if ~isempty(bs)
            bt = bs(1);
            at = max(as(as<=bs(1)));
            k=k+1;
            c=c+1;
            gaps(k) = bt-at;
            as(as<bt) = [];
            gt =  min([gaps(k),nP]);
            traj(t,bt); % this should be less than min
            t0 = max([bt-gt,1]);
            passA = fliplr(traj(t,t0:bt));
            passageAligned(c,1:length(passA)) = passA;
            traceNum(c) = t;
            traceStart(c) = t0;
        end
    end
    gaps(gaps==0) = [];
    if ~isempty(gaps)
        ptPerTrace{t} = gaps;
    end
end
c1=c;
traceLength = cat(1,ptPerTrace{:});
passageAligned = passageAligned(1:c1,:);
traceNum = traceNum(1:c1);
traceStart = traceStart(1:c1); 
%%


% [pa,ismiss] = fillmissing(passageAligned,"linear",2,"EndValues","none");
% nMiss = sum(ismiss,2);
% 
% a = diff(pa,[],2);    
% an = a; an(a>0) = nan; % back steps
% ap = a; ap(a<0) = nan; % closer
% b = abs(nansum(an,2));
% 
% [~,idx] = sort(b); % sort on least amount of backstep
% pa2 = passageAligned(idx,:);
% 
% %   best =1:c1;
% best =find( nansum(ap,2) > pars.d_high & abs(nansum(an,2)) < .2*pars.d_high & nMiss./traceLength < .1 & traceLength > pars.minPoints );
% % figure(3); clf; plot(passageAligned(1:6,:)','.-'); 
% box off; 
% %         title(['processive motion, dis=',num2str(dis_kb1(e))]); ylabel('3D dist. (nm)');
% xlabel('frame (0.5s)')
% proc_trace =500-1+traceNum(best);
% proc_t_start= traceStart(best);
% proc_t_end = traceStart(best)+traceLength(best);
