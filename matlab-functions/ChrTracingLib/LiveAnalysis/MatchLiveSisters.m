function [spotStacks,swapCount] = MatchLiveSisters(spotStacks,varargin)
%% Decide if we should swap sisters
% Rationale: 
% Assign spots in both channels to sister 1 or sister 2
% Sister1 is the original shortest-path trace through the dots
% Sister2 is any 2nd spot within a cut-off distance of the traced sister.  
% Spots that split in Z are handled in a second pass, if their x/y trajectories 
% are spatially divergent in the separate planes, they are likely two separate, 
% freely moving sisters just separated in z.  The shorter divergent trace data is 
% merged into sister2 and removed from sister1.  
% Note: G1 cells will have mostly empty sister2 traces (some occasional background
% spots filter in when we have 18,000 frames of observations), and stray background
% shouldn't correlate with data or other background in its 2d trajectory.  
% 
% Note: this version also removes jumps from the output data -- the sister
% trajectories haven't had any max step filtering yet.
% 
%% Inputs
% spotStacks nCells x 4, cell array of data for: chn1-sister1,
%           chn1-sister2, chn2-sister1, and chn2-sister2
%           each array is T-obs x 8-dim data (x,y,z,t,h,b,sigma,signficance)

defaults = cell(0,3);
defaults(end+1,:) = {'w','positive',30}; %  interpolation window for aligning 2-color data.  Must be at least as large as zDepth
defaults(end+1,:) = {'minFracLen','fraction',0.6}; % second sister must be at least this fraction of the the length of the first sister to be an acceptable alterantive match -- a good correlation for only a short window isn't good enough
defaults(end+1,:) = {'outlierFoldChange','positive',5}; % fold-change relative to median step-size for a point to be considered an outlier and ignored when computing correlation.  
defaults(end+1,:) = {'maxStep','positive',3}; % max single-step in xy pixels (relative to last observation) 
defaults(end+1,:) = {'npp','positive',108}; % nm per pixel (just for plotting/ troubleshooting)
defaults(end+1,:) = {'minLength','integer',300}; % traces must have at least this many points in common to consider switching
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'showFig','integer',0}; % note, 
defaults(end+1,:) = {'showSwaps','integer',0}; % note, 
pars = ParseVariableArguments(varargin,defaults,mfilename);

% s=88  % 30 88
nSpots = size(spotStacks,1);
swapCount = zeros(2,2);
for s=1:nSpots


    % Due to differences in z-position, the two channels are not always
    % in the same frame, but we can still use that data to measure
    % distances.  To avoid this problem, for the purposes of estimated
    % data alignment, we interpolate small gaps with a moving average
    % (this is much faster than regression). 
    % if we only interp one channel we also avoid inflating the data
    % volume: after filling in the gaps in chn2, we then remove all the
    % data where chn1 has gaps, so we still in general will keep only 1
    % observation per z
        w = pars.w; % 30; 
        d = 1; % dimension
        r1 = spotStacks{s,1}(:,d); % chn1 sisters 1
        r2 = spotStacks{s,3}(:,d); % chn 1 sister 2
        % r1 = fillmissing(spotStacks{s,1}(:,d),'movmean',w,'maxGap',w); % chn1 
        s1 = fillmissing(spotStacks{s,2}(:,d),'movmean',w,'maxGap',w); % chn2 sis1 
        s2 = fillmissing(spotStacks{s,4}(:,d),'movmean',w,'maxGap',w);  % chn2 sis2 
        
        % remove large single step outliers

        % REWRITE ME: remove nans before taking diff in order to remove
        % large steps. be careful to not skrew up computing overlap windows
        % either. 

        T = length(r1);
        ds2 = diff(s2);  mds2 = nonzeros(ds2); mds2 = min(pars.outlierFoldChange*nanmedian(abs(mds2)),pars.maxStep); %  5*median(abs(ds2(nonzeros(ds2))));
        % large step out         or       large step in 
        o2 = [0; ds2] > mds2 | [ds2; 0] > mds2; 
        s2(o2) = nan;
        c2s2 = 1:T; c2s2=c2s2(o2);% flag points to remove

        % same for s1 and r1
        
        ds1 = diff(s1);  mds1 = nonzeros(ds1); mds1 = min(pars.outlierFoldChange*nanmedian(abs(mds1)),pars.maxStep); % 
        
        o1 = [0; ds1] > mds1 | [ds1; 0] > mds1;
        s1(o1) = nan;
        c2s1 = 1:T; c2s1=c2s1(o1);% flag points to remove
        dr1 = diff(r1);  mdr1 = nonzeros(dr1); mdr1 = min(pars.outlierFoldChange*nanmedian(abs(mdr1)),pars.maxStep); %  
        or1 = [0; dr1] > mdr1 | [dr1; 0] > mdr1;
        r1(or1) = nan;
        c1s1 = 1:T; c1s1=c1s1(or1);% flag points to remove
        dr2 = diff(r2);  mdr2 = nonzeros(dr2); mdr2 = min(pars.outlierFoldChange*nanmedian(abs(mdr2)),pars.maxStep); % 
        or2 = [0; dr2] > mdr2 | [dr2; 0] > mdr2;
        r2(or2) = nan;    
        c1s2 = 1:T; c1s2=c1s2(or2);% flag points to remove
        keep1 = ~isnan(r1) & ~isnan(s1) & ~isnan(s2);
        keep2 = ~isnan(r2) & ~isnan(s1) & ~isnan(s2);
       
        % compute correlation
        corr_rs = zeros(2,2);
        leng_rs = zeros(2,2);
        corr_rs(1,1) = CorCoef(r1(keep1) - mean(r1(keep1)),s1(keep1) - mean(s1(keep1)));
        corr_rs(1,2)  = CorCoef(r1(keep1) - mean(r1(keep1)),s2(keep1) - mean(s2(keep1)));
        corr_rs(2,1)  = CorCoef(r2(keep2) - mean(r2(keep2)),s1(keep2) - mean(s1(keep2)));
        corr_rs(2,2) = CorCoef(r2(keep2) - mean(r2(keep2)),s2(keep2) - mean(s2(keep2)));
        % also compute length of overlapped data
        leng_rs(1,1) = sum( ~isnan(r1) & ~isnan(s1));
        leng_rs(1,2) = sum( ~isnan(r1) & ~isnan(s2));
        leng_rs(2,1) = sum( ~isnan(r2) & ~isnan(s1));
        leng_rs(2,2) = sum( ~isnan(r2) & ~isnan(s2));
        % corr_rs
        % leng_rs

        if pars.showFig > 0
            figure(pars.showFig); clf;
            subplot(2,1,1);
            plot(r1(keep1) - mean(r1(keep1))); hold on;
            plot(s1(keep1) - mean(s1(keep1))); 
            plot(s2(keep1) - mean(s2(keep1)));
            legend('chn1-s1','chn2-s1','chn2-s2')
            title('mean subtracted traces');
            xlabel('time (only co-observed windows)')
            subplot(2,1,2);
            plot(r2(keep2) - mean(r2(keep2))); hold on;
            plot(s1(keep2) - mean(s1(keep2))); 
            plot(s2(keep2) - mean(s2(keep2))); 
            legend('chn1-s2','chn2-s1','chn2-s2');
            xlabel('time (only co-observed windows)')
            pause(.1);
        end
        % plot(find(o2),s2(o2),'ro'); % show outliers;

        spotStacks{s,1}(c1s1,:) = nan;
        spotStacks{s,2}(c2s1,:) = nan;
        spotStacks{s,3}(c1s2,:) = nan;
        spotStacks{s,4}(c2s2,:) = nan;
        valid =   max(leng_rs(:)) > pars.minLength &&  max(corr_rs(:)) > 0.2;
        % choose pair to keep
        %  Highest correlation and within x% of the maximum length.
        %    x = 75% or 60% -- we don't want a super short trace that is
        %    better correlated.  
        % 1,1
        if corr_rs(1,1) >= max(corr_rs(:)) && leng_rs(1,1) >= pars.minFracLen*max(leng_rs(:)) && valid
            chn1a = spotStacks{s,1};
            chn2a = spotStacks{s,2};
            chn1b = spotStacks{s,3};
            chn2b = spotStacks{s,4};
            swapCount(1,1) = swapCount(1,1) + 1;
            swapped = false;
        % 1,2
        elseif all(corr_rs(1,2) >= corr_rs(:)) && leng_rs(1,2) >= pars.minFracLen*max(leng_rs(:))  && valid
            chn1a = spotStacks{s,1};
            chn2a = spotStacks{s,4};
            chn1b = spotStacks{s,3};
            chn2b = spotStacks{s,2};
            swapCount(1,2) = swapCount(1,2) + 1; 
            swapped = true; 
        % 2,1
        elseif all(corr_rs(2,1) >= corr_rs(:)) && leng_rs(2,1) >= pars.minFracLen*max(leng_rs(:))  && valid
            chn1a = spotStacks{s,3};
            chn2a = spotStacks{s,2};
            chn1b = spotStacks{s,1};
            chn2b = spotStacks{s,4};
            swapCount(2,1) = swapCount(2,1) + 1;
            swapped = true; 
        % 2,2
        elseif all(corr_rs(2,2) >= corr_rs(:)) && leng_rs(2,2) >= pars.minFracLen*max(leng_rs(:))  && valid
            chn1a = spotStacks{s,3};
            chn2a = spotStacks{s,4};
            chn1b = spotStacks{s,1};
            chn2b = spotStacks{s,2};
            swapCount(2,2) = swapCount(2,2) + 1;
            swapped = true; 
        else
            chn1a = spotStacks{s,1};
            chn2a = spotStacks{s,2};
            chn1b = spotStacks{s,3};
            chn2b = spotStacks{s,4};
            swapCount(1,1) = swapCount(1,1) + 1;
            swapped = false;            
        end
        
        % apply the swap if needed. The longest/best matched data is in 1,2
        spotStacks{s,1} = chn1a;
        spotStacks{s,2} = chn2a;
        spotStacks{s,3} = chn1b;
        spotStacks{s,4} = chn2b;


        if pars.showSwaps > 0 && swapped
            figure(pars.showSwaps); clf;
            subplot(2,1,1);
            plot(r1(keep1) - mean(r1(keep1))); hold on;
            plot(s1(keep1) - mean(s1(keep1))); 
            plot(s2(keep1) - mean(s2(keep1)));
            legend('chn1-s1','chn2-s1','chn2-s2')
            title(['s=',num2str(s),'  mean subtracted traces']);
            xlabel('time (only co-observed windows)')
            subplot(2,1,2);
            plot(r2(keep2) - mean(r2(keep2))); hold on;
            plot(s1(keep2) - mean(s1(keep2))); 
            plot(s2(keep2) - mean(s2(keep2))); 
            legend('chn1-s2','chn2-s1','chn2-s2');
            xlabel('time (only co-observed windows)');
            disp(corr_rs)
            disp(leng_rs)
            pause();
            % pause(.1);
        end

end
