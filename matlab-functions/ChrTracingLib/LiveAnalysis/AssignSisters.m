function [sisterRef,sisterConf] = AssignSisters(spotTrace_s,varargin)
% use all data from a single channel (both sisters, all z) to assign the
% sister identies, using the minimum jumpiness / continuum hypothesis.
% 
% Approach


defaults = cell(0,3); % 
defaults(end+1,:) = {'minSisObs','nonnegative',0.2}; % 
defaults(end+1,:) = {'minSisWindows','nonnegative',6}; % 
defaults(end+1,:) = {'movAveWin','positive',40}; % 
defaults(end+1,:) = {'minObsPerAve','positive',2}; % min number of observations in moving average window to count (per chn, per z)
defaults(end+1,:) = {'subSample','positive',2}; % sampling interval per moving average window (values must be from 1 to size moveAveWin)
defaults(end+1,:) = {'maxDistToSeed','positive',30}; % 
defaults(end+1,:) = {'maxDistToPrev','positive',10}; %
defaults(end+1,:) = {'maxSisterDist','positive',10}; %
defaults(end+1,:) = {'figSis','nonnegative',0};
defaults(end+1,:) = {'figMovAve','nonnegative',0};
defaults(end+1,:) = {'verbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);


try % for error handling

% some parameter short-hand
tw = pars.movAveWin;% 40; % window size for moving average
maxDistToSeed = pars.maxDistToSeed; 
maxDistToPrev = pars.maxDistToPrev;
zDepth = max(spotTrace_s(:,3,1));
tObs = size(spotTrace_s,1);
tzObs = tObs/zDepth;

% re-organizing data into z-stacks for easier processing
currStack1_s1 = nan(tzObs,8,zDepth);
currStack1_s2 = nan(tzObs,8,zDepth);
for z =1:zDepth
    currStack1_s1(:,:,z) = squeeze(spotTrace_s(z:zDepth:end,:,1));
    currStack1_s2(:,:,z) = squeeze(spotTrace_s(z:zDepth:end,:,2));
end
d=1;
obsPerZ_1 = squeeze(sum(~isnan(currStack1_s1(:,d,:)),1));
obsPerZ_2 = squeeze(sum(~isnan(currStack1_s2(:,d,:)),1));

tot_obs_s1 = sum(obsPerZ_1) ;
tot_obs_s2 = sum(obsPerZ_2);
% 
% % here we look for data in non-consecutive z-frames as evidence of a
% % sister doublet
% zsplits = zeros(tzObs,zDepth-2);
% for z = 1:zDepth -2; 
% zsplits(:,z) = ~isnan( currStack1_s1(:,1,z) ) ...
%             & isnan( currStack1_s1(:,1,z+1) ) ...
%             & ~isnan( currStack1_s1(:,1,z+2) );
% end
% tot_zsplits = sum(sum(zsplits,1));

% maybe divergent xy positions across a z-stack would be better
zsplits = zeros(tzObs,zDepth-1);
for z = 1:zDepth-1
zsplits(:,z) =  sqrt(sum( (currStack1_s1(:,1:2,z) - currStack1_s1(:,1:2,z+1)).^2,2))     >5;
end
tot_zsplits = sum(zsplits(:));

tot_sister = tot_obs_s2 + tot_zsplits; 

noSister = tot_sister < pars.minSisObs*tot_obs_s1;

% --- Pt 1
% moving average window for each z and each sister, with a minObs filter
%   the ideal sister paths will go through the moving average data, as a
%   way of avoiding getting confused by outlier noise events / detection of
%   more distal spots that drift into frame.
% Key parameters here: 
%    The size (in frames) of the moving average window tw (pars.movAveWin). 
%    The subsampling of the window (pars.subSample), 
%    The number of observations to be considered a legit median 
ts = 1:round(tw/pars.subSample):tzObs-tw;
tStps = length(ts);
t1s1 = nan(tStps,zDepth,2);
t1s2 = nan(tStps,zDepth,2);
for z=1:zDepth
    k=0;
    for t = ts
        b1 = t:t+tw;
        k=k+1;
        for d=1:2          
            l1 = sum(~isnan(currStack1_s1(b1,d,z)));
            if l1 >= pars.minObsPerAve
                t1s1(k,z,d) = nanmedian(currStack1_s1(b1,d,z)); 
            end
            l2 = sum(~isnan(currStack1_s2(b1,d,z)));
            if l2 >= pars.minObsPerAve
                t1s2(k,z,d) = nanmedian(currStack1_s2(b1,d,z));
            end
        end
    end
end

% % plotting the averages (just for troubleshooting)
% figure(1); clf;
% cmapz1 = hsv(zDepth)*.5;
% cmapz2 = cmapz1+.5;
% cmapz3 = jet(zDepth)*.5;
% for d=1:2
%     subplot(2,1,d);
%     for z=1:zDepth
%         plot(squeeze(currStack1_s1(:,d,z)),'.-','color',cmapz1(z,:)); hold on;
%         plot(squeeze(currStack1_s2(:,d,z)),'.-','color',cmapz2(z,:)); hold on;
%         plot(ts,t1s1(:,z,d),'+','color',cmapz3(z,:)); hold on;
%         plot(ts,t1s2(:,z,d),'x','color',cmapz3(z,:));
%     end
% end

% ---- Pt 2 -  % remove the average motion (largely cellular)
%   we'll add back exactly what we remove, so this doesn't need to be very
%   precise, but it does need to be smooth. 
%   Getting rid of the extreme and correlated step-to-step motion makes it
%   much easier / more robust to match the steps.  Otherwise in fast,
%   coherent motion, (e.g. two parallel spots both zooming to +x) in the
%   next frame it is easy for the other spot to be closer to you.
% 


tr1 = cat(2,t1s1,t1s2);  % t x z x d

xy = nan(tStps,2);
for d=1:2
xy(:,d) = smooth( nanmedian(tr1(:,:,d),2),.9,'rloess');
end


% tr1a = tr1; % preserving an original for troubleshooting
% figure(11); clf; % show the smooth data before subtraction
% for d=1:2
%     subplot(2,1,d); plot(xy(:,d)); hold on;
%     plot(tr1a(:,:,d),'.');
% end

for d=1:2
    tr1(:,:,d) = tr1(:,:,d) -  repmat(xy(:,d),1,2*zDepth);
    t1s1(:,:,d) = t1s1(:,:,d)-  repmat(xy(:,d),1,zDepth);
    t1s2(:,:,d) = t1s2(:,:,d)-  repmat(xy(:,d),1,zDepth);
end

% figure(10); clf;  % show the motion-corrected paths
% for d=1:2
%     subplot(2,1,d); % plot(xy(:,d));
%     hold on;  plot(tr1(:,:,d),'.');
% end


%  %   4/1/24 -- replaced this with moving average subtraction
% just use the first 1/4th of non-empty data to find the starting point
%   We need a starting point. Using median of all the data is problematic,
%   as we aren't guarenteed that the one spot doesn't move to be closer to
%   the second spots starting point by the time the second /sister spot
%   appears, in the case where both were not visible from the beginning. 
hasData1 = find(~isnan(nanmedian(t1s1(:,:,1),2)));
first_quarter_with_data = hasData1(1):hasData1(round(length(hasData1)/4));
startXY_1 = squeeze(nanmedian(nanmedian(t1s1(first_quarter_with_data,:,:),2),1))';
hasData2 = find(~isnan(nanmedian(t1s2(:,:,1),2)));
if length(hasData2) < pars.minSisWindows || noSister % no convincing sister
    refPts = startXY_1;
    startXY_2 = [NaN,NaN];
else
    first_quarter_with_data = hasData2(1):hasData2(round(length(hasData2)/4));
    startXY_2 = squeeze(nanmedian(nanmedian(t1s2(first_quarter_with_data,:,:),2),1))';
    refPts = [startXY_1; startXY_2];
end


% ----- Pt3 now we construct our outlier resistent shortest trace for both
% sisters.
% 
% % -- 3.1. Identify some starting points for the trace
% %    with motion correction complete, we can just use the medians
% %      don't use means here as they can get distracted by other data that
% %      drift into frame for a short period, and the path will lie in the
% %      distant middle far from either dataset. 
% startXY_1 = squeeze(nanmedian(nanmedian(t1s1(:,:,:),2),1))';
% startXY_2 = squeeze(nanmedian(nanmedian(t1s2(:,:,:),2),1))';
% if isnan(startXY_2(1))
%     refPts = startXY_1;
% elseif isnan(startXY_1(1))
%     refPts = startXY_2;
% else
%     refPts = [startXY_1; startXY_2];
% end

% -- 3.2  slide along window and construct trace
movAve = nan(2,tStps,2);
movAveRef = nan(2,tStps,2);
maxDistToLast = repmat(maxDistToSeed,2,1);
maxDistToPrev = repmat(maxDistToPrev,2,1);
% figure(10); clf; xlim([1,tStps]);
 for t=1:tStps % t=20
    currPts = squeeze(tr1(t,:,:));   % Nx2
    currPts = currPts(~isnan(currPts(:,1)),:);
    if size(currPts,1) >= 3
        [~,dis1] = knnsearch(currPts,currPts,'K',2);
        reject = dis1(:,2)>pars.maxSisterDist; % remove outliers before matching     % SHOULD consider adding this to LinkXYdoublets as well, though it slows things down 
        currPts(reject,:) = [];
        if size(currPts,1) >= 3
            [~,newCnts] = kmeans(currPts,2);
            currPts = newCnts;
        end
    end
    if ~isempty(currPts)
        % match both points with the least deplacement step-to-step
        %   (really want a criteria to position two traces that minimize
        %   the distance of all localization to one of the two traces, with
        %   an ability to ignore large jumps
        maxDist_t = max(maxDistToLast); 
        [matched,cost] = MatchPoints(currPts,refPts,'maxDist',maxDist_t);  % for speed, won't attempt to compute distance greater than maxDist
        % % if your match is bigger than a cut-off distance, you drop it
        if size(matched,1)==2
            m =matched(cost<maxDistToLast(matched(:,2)),:); 
        else
            m =matched(cost<maxDist_t,:); 
        end

        if ~isempty(m)
            movAve(m(:,2),t,1:2) = currPts(m(:,1),1:2);
            refPts(m(:,2),1:2) = currPts(m(:,1),1:2);
            maxDistToLast(m(:,2)) = maxDistToPrev(m(:,1));  
        end     
    end
    movAveRef(1:size(refPts,1),t,:) = refPts;
    % % % % for debugging - real time trace reconstruction
    % figure(10); 
    % plot(t,refPts(1,1),'b.'); hold on;
    % plot(t,refPts(2,1),'r.') ; hold on;
    % pause(.01);
 end

if pars.figMovAve > 0
    figure(pars.figMovAve); clf;
    plot(ts,movAveRef(1,:,d),'b.'); hold on;
    plot(ts,movAveRef(2,:,d),'r.') ; hold on;
    pause(.01);
end

 % 3.3  add back the moving/cell / chromosome frame
movAve (1,:,:) = squeeze(movAve (1,:,:)) + xy;
movAve (2,:,:) = squeeze(movAve (2,:,:)) + xy;
for d=1:2
t1s1(:,:,d) = t1s1(:,:,d)+  repmat(xy(:,d),1,zDepth);
t1s2(:,:,d) = t1s2(:,:,d)+  repmat(xy(:,d),1,zDepth);
end

% --- 3.4, Optional display results
if pars.figSis > 0
    figure(pars.figSis); clf;
    for d=1:2
        subplot(2,1,d);
        plot(ts,squeeze(movAve(1,:,d))','.-','color',[.8 .8 1],'linewidth',12,'MarkerSize',25); hold on;
        plot(ts,squeeze(movAve(2,:,d))','.-','color',[1 .8 .8],'linewidth',12,'MarkerSize',25); hold on;
        plot(squeeze(currStack1_s1(:,d,:)),'.','MarkerSize',4); hold on;
        plot(squeeze(currStack1_s2(:,d,:)),'.','MarkerSize',4); hold on;
        plot(ts,t1s1(:,:,d),'+','MarkerSize',4); hold on;
        plot(ts,t1s2(:,:,d),'x','MarkerSize',4);
    end
end


% ------ Pt 5: compute confidence (amount of data) supporting each step
%   This is returned as an optional output
%   Not currently used, though this is not very comp. expensive. 
stpZ = tw/pars.subSample;
movConf = nan(2,tStps);
for t = 1:tStps
    b1 = (t-1)*stpZ+1:t*stpZ;
    pts_x = cat(1,currStack1_s1(b1,1,:),currStack1_s2(b1,1,:));
    pts_y = cat(1,currStack1_s1(b1,2,:),currStack1_s2(b1,2,:));
    pts_xy = [pts_x(~isnan(pts_x)),pts_y(~isnan(pts_y))];
    np = size(pts_xy,1);
    dis1 = sqrt(sum( (pts_xy - repmat(squeeze(movAve(1,t,:))',np,1)).^2, 2));
    dis2 = sqrt(sum( (pts_xy - repmat(squeeze(movAve(2,t,:))',np,1)).^2, 2));
    movConf(1,t) = sum(dis1<3);
    movConf(2,t) = sum(dis2<3);
end

% --- Pt 6, Convert back to image frames / absolute time
% (rather than indices of frame averages). This makes it easier to plot
sisterRef = nan(2,tObs,2);
sisterConf = nan(2,tObs);
stp = zDepth*tw/pars.subSample;
to = 1:stp:(tObs-tw*zDepth);
sisterRef(:,to,:) = movAve;
sisterConf(:,to) = movConf;
for s=1:2
    sisterConf(s,:) = fillmissing(sisterConf(s,:),'previous','maxGap',stp);
    for d=1:2  % 
        sisterRef(s,:,d) = fillmissing(sisterRef(s,:,d),'linear','EndValues','nearest');  % fillmissing(sisterRef(s,:,d),'linear','maxGap',5*tObs);  % 
        % sisterRef(s,:,d) = fillmissing(sisterRef(s,:,d),'movmean',2*stp,'maxGap',4*stp);  % fillmissing(sisterRef(s,:,d),'linear','maxGap',5*tObs);  % 
    end
end

% 
% figure(3); clf;
% for d=1:2
%     subplot(2,1,d);
%     plot(spotTrace_s(:,d,1),'o','color',[.5 .5 1 .5], 'MarkerSize',2); hold on;
%     plot(spotTrace_s(:,d,2),'o','color',[1 .5 1 .5],'MarkerSize',2);
%     plot(sisterRef(1,:,d),'-','color',[0 0 .5]); hold on;
%     plot(sisterRef(2,:,d),'-','color',[0 .5 .5]);
% end


if pars.figSis > 0
    figure(pars.figSis); 
    for d=1:2
        subplot(2,1,d);
        plot(sisterRef(1,1:zDepth:end,d),'-','color',[0 0 .5]); hold on;
        plot(sisterRef(2,1:zDepth:end,d),'-','color',[0 .5 .5]);
    end
end

catch er
    warning(er.getReport);
    disp('debug');
end
