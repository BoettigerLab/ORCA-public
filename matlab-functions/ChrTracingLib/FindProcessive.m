function [traces,traceIdStartEnd,traceStats] = FindProcessive(dis3,varargin)
% 

defaults = cell(0,3); % 
% shared parameters for moving average
defaults(end+1,:) = {'minD','positive',[]}; % miniumum separation for a passage event, empty for auto-compute 2x median
defaults(end+1,:) = {'maxD','positive',inf}; % require the proximal state be within this distance
defaults(end+1,:) = {'minProx','positive',inf}; % require the proximal state be within this distance
defaults(end+1,:) = {'maxLen','positive',60}; % max number of observations for the passage event
defaults(end+1,:) = {'minLen','positive',20}; % miniumum number of observations for the passage event
defaults(end+1,:) = {'maxMiss','positive',0.1}; % maximum fraction missing data tolerated 
defaults(end+1,:) = {'maxSingle','positive',0.24}; % maximum fraction of total path in single step
defaults(end+1,:) = {'maxBackStep','positive',0.2}; % maximum fraction of back-step events tolerated
defaults(end+1,:) = {'maxBackDist','positive',0.5}; % maximum fraction of the passage distances that can be backwards
defaults(end+1,:) = {'showPlots','integer',0}; % maximum fraction of the passage distances that can be backwards
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);




% dis3 = dis3D{2}{7};
% pars.minD = 1000;
% pars.minLen = 40;
% pars.maxBack = 0.1;
% pars.maxMiss  = 0.1;


[nC,T] = size(dis3);
if isempty(pars.minD)
    pars.mind = 2*nanmedian(dis3(:));
end

% currTrace = dis3(c,:);  figure(3); clf; plot(currTrace,'.-');
[pa,ismiss] = fillmissing(dis3,"linear",2,"EndValues","none"); % dim=2  
a = diff(pa,[],2);    
an = a; an(a>0) = nan; % back steps
ap = a; ap(a<0) = nan; % closer
good = false(nC,T);

% good1= false(nC,T);
% good2= false(nC,T);
% good3= false(nC,T);
% good4= false(nC,T);

traceIdStartEnd = zeros(1000,3);
traceStats = struct();
traceStats(1000).nMiss = 0;
k=0;

if pars.showPlots
    figure(pars.showPlots); clf;
end

c0 = 0; t00 = 0; t10 = 0;
for c=1:nC
    for t=1:T-pars.maxLen
        % passageSize = max(pa(c,ts))-min(pa(c,ts));
        tt = t:t+pars.maxLen-1;
        [d1,t1] = max(pa(c,tt));
        [d0,t0] = min(pa(c,tt));
        passageSize = d1-d0; % max to min distance
        if t0-t1 > pars.minLen && passageSize > pars.minD  %  (doesn't truly need to be min to max)
            ts = t+(t1:t0)-1;
             
            nMiss = sum(ismiss(c,ts));
            closingSteps = -nansum(an(c,ts),2);
            openingSteps = nansum(ap(c,ts),2);
            numClosingSteps = sum(~isnan(an(c,ts)));
            numOpeningSteps = sum(~isnan(ap(c,ts)));  
            totSteps = numOpeningSteps + numClosingSteps;
            tr = dis3(c,ts+1);
            obsSteps = tr(~isnan(tr));
            stepSizes =  abs(diff(obsSteps));

            good(c,t+t1) = (passageSize > pars.minD) && ...
                (openingSteps < pars.maxBackDist*passageSize) && ...
                (numOpeningSteps < pars.maxBackStep*totSteps) && ...
                (nMiss/pars.minLen < pars.maxMiss) && ...
                max(stepSizes) < pars.maxSingle*passageSize && ...    %  single steps no more than x% of path (hardly 'processive' if it does it mostly in 1 jump) 
                d0 < pars.minProx && ...
                d1 < pars.maxD;
            if good(c,t+t1) && ~(c==c0 && t+t1-1==t00 && t+t0-1==t01 )
                k=k+1;
                traceIdStartEnd(k,1) = c;
                traceIdStartEnd(k,2) = t+t1-1;
                traceIdStartEnd(k,3) = t+t0-1;
                traceStats(k).passageSize = passageSize;
                traceStats(k).openingSteps = openingSteps;
                traceStats(k).closingSteps = closingSteps;
                traceStats(k).numClosingSteps = numClosingSteps;
                traceStats(k).numOpeningSteps = numOpeningSteps;
                traceStats(k).totSteps = totSteps;
                traceStats(k).stepSizes = stepSizes;
                traceStats(k).nMiss = nMiss;
                traceStats(k).d0 = d0;
                traceStats(k).d1 = d1;
                c0 =c ; t00 = t+t1-1; t01 = t+t0-1; 
                if pars.showPlots
                    figure(pars.showPlots); hold on; plot(dis3(c,t+t1-1:t+t0-1),'.-'); pause(.01);
                end
            end
            % good1(c,t) = (passageSize > pars.minD);
            % good2(c,t) = (openingSteps < pars.maxBackDist*passageSize);
            % good3(c,t) = (numOpeningSteps < pars.maxBackStep*numClosingSteps);
            % good4(c,t) = (nMiss/pars.minLen < pars.maxMiss);
        end
    end
end

traceIdStartEnd_full = traceIdStartEnd;
traceStats_full = traceStats; 

drop = traceIdStartEnd(:,1) == 0;
traceIdStartEnd(drop,:) = [];
traceStats(drop) = [];

% remove double counts
nT = size(traceIdStartEnd,1);

remove_id = false(1,nT); 
for t=1:nT
    sameCell = traceIdStartEnd(:,1)==traceIdStartEnd(t,1); %
    overlapTime = traceIdStartEnd(:,2) < traceIdStartEnd(t,2)+pars.maxLen ...
        & traceIdStartEnd(:,2) > traceIdStartEnd(t,2)-pars.maxLen;
    keep_first = find(sameCell & overlapTime,1);
    remove_id(sameCell & overlapTime) = true;
    remove_id(keep_first) = false;
end

traceIdStartEnd(remove_id,:) = [];
traceStats(remove_id) = [];


% figure(4); clf;
nT = size(traceIdStartEnd,1);
nTs = max(traceIdStartEnd(:,3)-traceIdStartEnd(:,2));
traces = nan(nT,nTs);
for t=1:nT
    currTrace = dis3(traceIdStartEnd(t,1),traceIdStartEnd(t,2):traceIdStartEnd(t,3));
    traces(t,1:length(currTrace)) = currTrace;
    % figure(4); plot(currTrace,'.-'); hold on;
end

