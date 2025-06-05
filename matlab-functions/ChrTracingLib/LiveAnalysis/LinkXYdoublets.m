function [spotTrace1,spotTrace2,spotsPerStep] = LinkXYdoublets(fits1,fits2c,xy1p,xy2p,zDepth,varargin)
%% Method
% Uses a moving average window to determine the approximate trajectory of
% the particle relative to noisy background localizations. This track 
% defines an spatial window in which additional localizations can be added 
% to the path.  This reference prevents the trace from getting 'off track' 
% by following a noisy interlude burst.
% 
% The moving average is constructed in 2 steps. First, we consider the full
% field of spots, and enforce unique linking of the clusters (defined by
% simple binning using hist3) between all frames.  As in LinkLiveStruct, a
% maximum step is enforced.
% 
% Perhaps it makes sense that closest point to the last step that is also 
% within the acceptable distance from the time average is the ONLY hard
% filter condition?  We don't also need to cut-off on the size of that step
% It doesn't hurt to leave these other parameters in the code.
% 
%% Inputs
%    fits1 - fit table (DaoSTORM data is much easier in matlab as a table)
%            corresponding to channel1
%   fits2c - registered fit table for channel2 (already mapped onto chn1)
%   xy1p   - x,y centers for fits1, already paired with channel 2 data.
%   xy2p   - x,y center for fits2 already paired with channel 1 data.
%   zDepth - number of steps per z
% 
% 
%% Optional Inputs
%    see below
% 
%% Outputs
%    spotTrace1 - nSpots x nFrames x d-dims 
%                 8 data dims:  x,y,z,frame,height,bkd,sigma,significance
%                 this is the channel 1 data. 
%    spotTrace2 - same as above, but for channel 2 data. 
%          (maybe we should combines these?)
%% Method
%   
% 
%% History / change log
%  Evolved from LinkLiveStruct / LinkLiveTable
%  * Enriched the concept of 'seed points' - we don't use only 1 initial
%  x,y as a seed, now we time average batch windows to get an idea of the
%  true trace.
%  * Does not enforce that each spot added to a chain is added to only 1
%   unique chain. (if the spots overlap it's reasonable they can share the
%   fit, and the frequency of overlap with 1 locus per cell is very low
%   anyway).
%  * Does not use MatchPoints algorithm to uniquely assign points
%      this will be slightly less robust to large drift jumps, which the
%      match points could reliably dissambiguate. But we find it more
%      helpful filtering to avoid large jumps, and again, these events are
%      rare (though 1 of the last 5 datasets clearly did have such a jump).
%  Challenge:
%     sometimes two spots are seen in a crop-box, and the time windows in
%     which they are detected have substantial non-overlap.  
%     This was previously handled by the unique match
% 
%% To consider for future
%  * should we consider data from both channels explicitly in deciding
%  which points to accept?
%  * should we add reporting on how many doublet frames were actually
%  detected? 


%% updates
%  evolved from LinkZ.m, this version uses both traces 


%% optional pars
defaults = cell(0,3); % 
% shared parameters for moving average
defaults(end+1,:) = {'movAveSteps','positive',10}; % Number of windows in which to split the trace in when computing the moving average. Note, multiplied by z-depth    
defaults(end+1,:) = {'movAveStepMaxFold','positive',3}; % max fold change in step size relative to local average, used by Remove Jumps
defaults(end+1,:) = {'moveAveMaxPixStep','positive',7}; % max absolute step size (in pixels) for moving avearage trace, used by Remove Jumps
defaults(end+1,:) = {'movAveLocalStep','positive',6}; % window size to determine local jumps  (both rounds use this filter)
% parameters for coarse linking of whole field
defaults(end+1,:) = {'movAveMaxGap','positive',4}; % Downsampling scale to determine spot clusters (spots more than this number of pixels apart will be considered as separate clusters)    
defaults(end+1,:) = {'minObsPerAve','positive',10}; % minimum number of observations in the moving average window to be counted as a valid obs in the coarse downsampling
defaults(end+1,:) = {'maxDistToSeed','positive',35}; % maximum distance to seed point to be counted in the coarse linking
% parameters for refined measurement of moving average
defaults(end+1,:) = {'windowR','positive',35}; % 'radius' of window around the original seed point in which valid localizations may occur
defaults(end+1,:) = {'maxDistFromRoughAve','positive',5}; % Distance in pixels from the rough moving average trace in constructing the refined trace (no binning this time, just a restricted xy window)
defaults(end+1,:) = {'minObsPerAve2','positive',6}; % minimum number of observations in the moving average window to be counted as a valid obs
defaults(end+1,:) = {'coarseInterpWindow','positive',20}; % window size of points used for interpolation filter (allows recovery in fine trace of time blocks not detected in coarse trace)
% parameters for trace assembly
defaults(end+1,:) = {'maxDistFromMovingAve','positive',4}; % Max distance an acceptable spot can be from its local-time-moving average centroid reported as fold-change relative to moving-average stepsize
defaults(end+1,:) = {'maxStep','positive',1.5}; % distance in pixels, the max distance from the last localization to be added to the trace (permisive)   
defaults(end+1,:) = {'maxStepFrame1','positive',30}; % distance in pixels, the max distance from the starting point to start counting (starting point is first non-nan in the time average)
% other parameters
defaults(end+1,:) = {'tLapseMov','array',0}; % raw image movie for overlay plotting (optional, 0 = don't use)
defaults(end+1,:) = {'figCoarseMovAve','integer',0}; % overlay of image data with the coarse moving  average plot.  Mostly good for troubleshooting. Displays once per FOV
defaults(end+1,:) = {'movAveFig','integer',1}; % per trace
defaults(end+1,:) = {'spotMapFig','integer',2}; % per trace
defaults(end+1,:) = {'traceFig','integer',3}; %  per trace
defaults(end+1,:) = {'dynFig','integer',0}; %% per step per trace
defaults(end+1,:) = {'nFrames','integer',0};
defaults(end+1,:) = {'npp','positive',108}; % nm per pixel (just for plotting)
defaults(end+1,:) = {'selSpots','integer',0}; % 
defaults(end+1,:) = {'parallel','boolean',false}; % 
defaults(end+1,:) = {'pause','boolean',false}; % 

defaults(end+1,:) = {'saveFigs','boolean',false}; % 
defaults(end+1,:) = {'overwrite','boolean',true}; % 

pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);
tLapseMov = pars.tLapseMov;

%% setup 
% moving average
tStps = zDepth*pars.movAveSteps; % needs to be an even multiple of zDepth
% cmap = jet(tStps);
if pars.nFrames==0
    nFrames = max(fits1.frame);
else
    nFrames = pars.nFrames;
end
tBlock = ceil(nFrames/tStps);

% initializing some empty variables and adding the frame data
nSpots = size(xy1p,1);
spotTrace1 = nan(nSpots,nFrames,8); % (s,tm,1)  % spot x time x dim (dims = x,y,z,frame,h,b,sigma,signficance)
spotTrace2 = nan(nSpots,nFrames,8); % (s,tm,1)  % spot x time x dim (dims = x,y,z,frame,h,b,sigma,signficance)
spotTrace1(:,:,4) = repmat( (1:nFrames),nSpots,1);
z = mod(spotTrace1(:,:,4),zDepth); z(z==0) = zDepth;
spotTrace1(:,:,3) = z;
spotTrace2(:,:,4) = repmat( (1:nFrames),nSpots,1);
z = mod(spotTrace2(:,:,4),zDepth); z(z==0) = zDepth;
spotTrace2(:,:,3) = z;

spotTrace1 = cat(4,spotTrace1,spotTrace1);
spotTrace2 = cat(4,spotTrace2,spotTrace2);

%% whole FOV moving average
% collapse blocks of time into single pile ups and fit centroids
%   Uniquely assign new points to seed points. 

if tLapseMov~=0
    [H,~,~]=size(tLapseMov);
    imW = H/pars.movAveMaxGap;
else
    imW = 2304/pars.movAveMaxGap;
end

movAveAll1 = nan(nSpots,tStps,2);
movAveAll2 = nan(nSpots,tStps,2);
x0 = 0; y0=0; % offset for zoom in
% cmap = jet(tStps);
startPoints = xy1p;
 for t=1:tStps % t=20
     % chn1
        frame0 = (t-1)*tBlock; % indexed from 0
        frame1 = min(t*tBlock,nFrames); % 1-block more or last-frame,  
        selFrames1 = fits1.frame >= frame0 & fits1.frame < frame1;
        tFits1 = fits1(selFrames1,:);
        xs1 = round((tFits1.x-x0)/pars.movAveMaxGap);
        ys1 = round((tFits1.y-y0)/pars.movAveMaxGap);
        obs = hist3([ys1,xs1],'Ctrs',{1:imW,1:imW});
        % figure(17); clf; imagesc(obs); colorbar; caxis([0,20]); hold on; 
        regs = regionprops(obs>pars.minObsPerAve,'Centroid'); % 'PixelIdxList',
        currPts = cat(1,regs.Centroid)*pars.movAveMaxGap;
        [matched,cost] = MatchPoints(currPts,startPoints,'maxDist',pars.maxDistToSeed);
        m =matched(cost<pars.maxDistToSeed,:); 
        if ~isempty(m)
            movAveAll1(m(:,2),t,1:2) = currPts(m(:,1),1:2);
        end

        % chn2
        frame0 = (t-1)*tBlock; % indexed from 0
        frame1 = min(t*tBlock,nFrames); % 1-block more or last-frame,  
        selFrames2 = fits2c.frame >= frame0 & fits2c.frame < frame1;
        tFits2 = fits2c(selFrames2,:);
        xs1 = round((tFits2.x-x0)/pars.movAveMaxGap);
        ys1 = round((tFits2.y-y0)/pars.movAveMaxGap);
        obs = hist3([ys1,xs1],'Ctrs',{1:imW,1:imW});
        regs = regionprops(obs>pars.minObsPerAve,'Centroid'); %
        currPts = cat(1,regs.Centroid)*pars.movAveMaxGap;
        [matched,cost] = MatchPoints(currPts,startPoints,'maxDist',pars.maxDistToSeed);
        m =matched(cost<pars.maxDistToSeed,:); 
        if ~isempty(m)
            movAveAll2(m(:,2),t,1:2) = currPts(m(:,1),1:2);
        end
 end

 pars.movAveCoarseMaxPixStep= 20;
 % even if spots werent 
w=pars.coarseInterpWindow; % 20;
moveAveAllInt1 =  movAveAll1;
moveAveAllInt2 =  movAveAll2;
for n=1:nSpots
    moveAveAllInt1(n,:,1:2) = RemoveJumps(squeeze(moveAveAllInt1(n,:,1:2)),'maxDiff',pars.movAveStepMaxFold,'localRegion',pars.movAveLocalStep,'maxAbsStep',pars.movAveCoarseMaxPixStep);
    moveAveAllInt2(n,:,1:2) = RemoveJumps(squeeze(moveAveAllInt2(n,:,1:2)),'maxDiff',pars.movAveStepMaxFold,'localRegion',pars.movAveLocalStep,'maxAbsStep',pars.movAveCoarseMaxPixStep);
    moveAveAllInt1(n,:,1) = fillmissing(squeeze(moveAveAllInt1(n,:,1)),'movmean',w,'EndValues','nearest');
    moveAveAllInt1(n,:,2) = fillmissing(squeeze(moveAveAllInt1(n,:,2)),'movmean',w,'EndValues','nearest');
    moveAveAllInt2(n,:,1) = fillmissing(squeeze(moveAveAllInt2(n,:,1)),'movmean',w,'EndValues','nearest');
    moveAveAllInt2(n,:,2) = fillmissing(squeeze(moveAveAllInt2(n,:,2)),'movmean',w,'EndValues','nearest');
end

if pars.figCoarseMovAve > 0
   
    figure(pars.figCoarseMovAve); clf;
    if tLapseMov ~=0
        Ncolor(tLapseMov); hold on; 
    end
    ts = 1:tStps;
    for n=1:nSpots % n=75
        scatter(squeeze(moveAveAllInt1(n,:,1)),squeeze(moveAveAllInt1(n,:,2)),[],ts','filled'); hold on; 
        xs = squeeze(moveAveAllInt1(n,:,1));
        ys = squeeze(moveAveAllInt1(n,:,2));
        plot(xs(~isnan(xs)),ys(~isnan(xs)),'.-','color',[.5 0 0],'lineWidth',2);
        scatter(squeeze(moveAveAllInt2(n,:,1)),squeeze(moveAveAllInt2(n,:,2)),[],ts','filled');
        xs = squeeze(moveAveAllInt2(n,:,1));
        ys = squeeze(moveAveAllInt2(n,:,2));
        plot(xs(~isnan(xs)),ys(~isnan(xs)),'.-','color',[0 0 1]);
        % text(nanmean(xs)+2,nanmean(ys),num2str(n),'color','y'); hold on;
    end
    % hold on; plot(xy1p(:,1),xy1p(:,2),'yo','MarkerSize',20); % display the seed points for sanity check    
    % hold on; plot(xy2p(:,1),xy2p(:,2),'gs','MarkerSize',20); % display the seed points for sanity check  
    sNum = cellstr(num2str((1:size(xy1p,1))'));
    if tLapseMov ~=0
    text(xy1p(:,1)+5,xy1p(:,2),sNum,'color','w'); hold on;
    else
        text(xy1p(:,1)+5,xy1p(:,2),sNum,'color','k'); hold on;
    end
end


% figure(17); clf; 
% subplot(1,2,1); imagesc(moveAveAllInt1(:,:,1)); title('interpolated moving average paths');
% subplot(1,2,2); imagesc(movAveAll1(:,:,1)); title('original moving average paths')



%% compute moving average
% 52 is good.  50 is pretty bad
if pars.selSpots == 0
    pars.selSpots =1:nSpots;
end
spotsPerStep = zeros(nSpots,nFrames,2);
%% linear  (good for plots and error checking)
if ~pars.parallel
    for s = pars.selSpots % s=1
        %%
        x0 = xy1p(s,1) - pars.windowR;
        x1 = xy1p(s,1) + pars.windowR;
        y0 = xy1p(s,2) - pars.windowR;
        y1 = xy1p(s,2) + pars.windowR;
        selPts = fits1.x > x0 & fits1.x < x1 & fits1.y > y0 & fits1.y < y1;
        sFit1 = fits1(selPts,:);
        selPts = fits2c.x > x0 & fits2c.x < x1 & fits2c.y > y0 & fits2c.y < y1;
        sFit2 = fits2c(selPts,:);
        
        % just plotting
        if pars.spotMapFig > 0
            figure(pars.spotMapFig); clf;
            if tLapseMov~= 0
                Ncolor(tLapseMov); hold on; axis image;
            end
            plot(xy1p(:,1),xy1p(:,2),'yo','MarkerSize',20); hold on;
            plot(xy2p(:,1),xy2p(:,2),'gs','MarkerSize',20);
            sNum = cellstr(num2str((1:size(xy1p,1))'));
            text(xy1p(:,1)+25,xy1p(:,2),sNum,'color','w'); hold on;
            xlim([x0,x1]); ylim([y0,y1]);
            title(['spot ',num2str(s)]);
            scatter(sFit1.x,sFit1.y,[],sFit1.frame,"filled"); colormap(jet); hold on;
            plot(sFit1.x,sFit1.y,'m.','MarkerSize',10); colormap(jet); hold on;
            scatter(sFit2.x,sFit2.y,[],sFit2.frame,"filled"); colormap(jet); hold on;
            plot(sFit2.x,sFit2.y,'c.','MarkerSize',10); colormap(jet); hold on;
        end
        
    
    
        movAve1 = nan(tStps,2);
        movAve2 = nan(tStps,2); 
        for t=1:tStps % t=9
            % to avoid the occassional stray-background point in an early frame
            % seeding the trace in an wrong position (which typically runs a short
            % and terminal path), we average the movie in blocks 
            % 
            % block averaging
            frame0 = (t-1)*tBlock; % indexed from 0
            frame1 = min(t*tBlock,nFrames); % 1-block more or last-frame,  
            sx = moveAveAllInt1(s,t,1);
            sy = moveAveAllInt1(s,t,2);
            w = pars.maxDistFromRoughAve;
    
            % get frames in time window and within distance from rough moving average  
            frameT = sFit1.frame >= frame0 & sFit1.frame < frame1 ;
            frameXY = (sFit1.x > (sx-w)) & (sFit1.x < (sx+w)) & (sFit1.y > (sy-w)) & (sFit1.y < (sy+w));
            selFrames1 = frameT & frameXY;
            if sum(selFrames1) > pars.minObsPerAve2
                fitsB1 = sFit1(selFrames1,:);
                movAve1(t,:) = [mean(fitsB1.x),mean(fitsB1.y)]; % time local mean 
            end
            % chn2 
            sx = moveAveAllInt2(s,t,1);
            sy = moveAveAllInt2(s,t,2);
            w = pars.maxDistFromRoughAve;
            selFrames2 = sFit2.frame >= frame0 & sFit2.frame < frame1  & ...
                 sFit2.x > sx - w & sFit2.x < sx + w & sFit2.y > sy - w & sFit2.y < sy+ w;
            if sum(selFrames2)  > pars.minObsPerAve
                fitsB2 = sFit2(selFrames2,:);
                movAve2(t,:) = [mean(fitsB2.x),mean(fitsB2.y)];  % local mean  
            end
        end
    
        % clean up moving average
        % use moving average to define the expected stepsize 
        % remove outliers in the moving average
        if pars.movAveFig > 0
             figure(pars.movAveFig);  % figure(15); 
             clf; Ncolor(tLapseMov);  axis image; hold on;
             scatter(movAve1(:,1),movAve1(:,2),[],(1:tStps)','filled'); colormap('jet'); hold on;
             plot(movAve1(:,1),movAve1(:,2),'.','color',[.5 0 0],'MarkerSize',10);
             hold on; scatter(movAve2(:,1),movAve2(:,2),[],(1:tStps)','filled'); colormap('jet')
             plot(movAve2(:,1),movAve2(:,2),'.','color',[0 0 1],'MarkerSize',10);
                    xlim([x0,x1]); ylim([y0,y1]);
                title(['all ave. spot ',num2str(s)]); pause(.1);
        end
    
        moveAveClean1 = RemoveJumps(movAve1,'maxDiff',pars.movAveStepMaxFold,'localRegion',pars.movAveLocalStep);
        hasDat1 = ~isnan(moveAveClean1(:,1));
        aveStep1 = sqrt(sum(diff(moveAveClean1(hasDat1,:)).^2,2));
        medAveStep1 = nanmedian(aveStep1); %#ok<*NANMEDIAN>
        moveAveClean2 = RemoveJumps(movAve2,'maxDiff',pars.movAveStepMaxFold,'localRegion',pars.movAveLocalStep);
        hasDat2 = ~isnan(moveAveClean2(:,1));
        aveStep2 = sqrt(sum(diff(moveAveClean2(hasDat2,:)).^2,2));
        medAveStep2 = nanmedian(aveStep2); 
       
       noData1 = sum(isnan(moveAveClean1(:,1))) == size(moveAveClean1,1);
       noData2 = sum(isnan(moveAveClean2(:,2))) == size(moveAveClean2,1);
       if noData1 || noData2 
           if pars.verbose
               disp(['unable to find trace for spot s=',num2str(s)]);
           end
           continue;
       end

        if pars.movAveFig > 0
             figure(pars.movAveFig); clf; Ncolor(tLapseMov);  axis image; hold on;
             x = moveAveClean1(:,1); y = moveAveClean1(:,2);
             plot(x(~isnan(x)),y(~isnan(y)),'-','color',[.5 0 0]); hold on;
             x = moveAveClean2(:,1); y = moveAveClean2(:,2);
             plot(x(~isnan(x)),y(~isnan(y)),'-','color',[0 0 1]); hold on;
             scatter(moveAveClean1(:,1),moveAveClean1(:,2),[],(1:tStps)','filled'); colormap('jet'); hold on;
             plot(moveAveClean1(:,1),moveAveClean1(:,2),'.','color',[.5 0 0],'MarkerSize',10);
             hold on; scatter(moveAveClean2(:,1),moveAveClean2(:,2),[],(1:tStps)','filled'); colormap('jet')
             plot(moveAveClean2(:,1),moveAveClean2(:,2),'.','color',[0 0 1],'MarkerSize',10);
                    xlim([x0,x1]); ylim([y0,y1]);
                title(['clean ave. spot ',num2str(s)]);
        end
    
    
        %% compute traces
        if pars.dynFig > 0
            im = tLapseMov(y0:y1,x0:x1,:);
        end
        
        try
            valid1 = moveAveClean1(~isnan(moveAveClean1(:,1)),:);
            valid2 = moveAveClean2(~isnan(moveAveClean2(:,1)),:);
            refPt1 = valid1(1,:); % start with first non-nan
            refPt2 = valid2(1,:);

            maxStep1 = pars.maxStepFrame1;
            maxStep2 = pars.maxStepFrame1;
            maxDistFromMovAve1 = pars.maxDistFromMovingAve*medAveStep1;
            maxDistFromMovAve2 = pars.maxDistFromMovingAve*medAveStep2;
            for z=1:zDepth
                refTime1 = find(~isnan(moveAveClean1(:,1)),1);
                refTime2 = find(~isnan(moveAveClean2(:,1)),1);

                for t=z:zDepth:nFrames 
                    % accept the closest spot to the current reference that is also
                    % within the required distance from the reference centroid. 
                    tt = floor((t-1)/tBlock)+1;
                    xy1 = movAve1(tt,:);
                    xy2 = movAve2(tt,:);
                    if ~isnan(xy1)
                        currZstack = sFit1.frame == t-1; 
                        fr1 = sFit1(currZstack,:);
                        if ~isempty(fr1)
                            newTime1 = t;
                            newPts = [fr1.x,fr1.y];
                            curDat1 = [fr1.height,fr1.background,fr1.xsigma,fr1.significance];
                            distLast = sqrt(sum((newPts - repmat(refPt1,size(newPts,1),1)).^2,2));
                            distMov = sqrt(sum((newPts - repmat(xy1,size(newPts,1),1)).^2,2));
                            [minDistLast,idMinDist] = min(distLast);  % only tracing a single path at this point so min distance is okay, don't need MatchPairs 
                            % sort by ascend and keep the top two
                            [~,idxS] = sort(distLast,'ascend');
                            if length(idxS) >= 2
                            spotTrace1(s,t,1:2,2) = newPts(idxS(2),:)';
                            spotTrace1(s,t,5:end,2) = curDat1(idxS(2),:)';
                            end
                            % may want to toss second doublet if first doublet is nan  
                            dt = sqrt((newTime1-refTime1)/zDepth);  % allow larger steps for larger time lags  
                            if minDistLast < maxStep1*dt && distMov(idMinDist) < maxDistFromMovAve1
                                if pars.dynFig>0 % dynamically and sequentially add spots  (really just for troubleshooting)
                                    figure(pars.dynFig); clf;
                                    Ncolor(im);  axis image; hold on;
                                    plot(refPt1(1)-x0,refPt1(2)-y0,'r+');
                                    plot(newPts(:,1)-x0,newPts(:,2)-y0,'w.');
                                    plot(newPts(idMinDist,1)-x0,newPts(idMinDist,1)-y0,'ro'); 
                                    title('chn1'); pause(.1); 
                                end
                                spotsPerStep(s,t,1) = size(newPts,1);
                                refPt1 =  newPts(idMinDist,:);
                                refTime1 = newTime1;
                                spotTrace1(s,t,1:2,1) = refPt1;
                                spotTrace1(s,t,5:end,1) = curDat1(idMinDist,:);
                                maxStep1 = pars.maxStep;
                            end
                        end   
                    end
                    if ~isnan(xy2)
                        currZstack = sFit2.frame == t-1; % 
                        fr2 = sFit2(currZstack,:);
                        if ~isempty(fr2)
                            newTime2 = t;
                            newPts2 = [fr2.x,fr2.y]; % z,t, h,bkd,sigma
                            curDat2 = [fr2.height,fr2.background,fr2.xsigma,fr2.significance];
                            distLast = sqrt(sum((newPts2 - repmat(refPt2,size(newPts2,1),1)).^2,2));
                            distMov = sqrt(sum((newPts2 - repmat(xy2,size(newPts2,1),1)).^2,2))  ;
                            % select
                            %  currently: closest point that is within the step size
                            %  tolerance and that is within the appropriate distance from
                            %  the moving average.  
                            [minDistLast,idMinDist] = min(distLast);% sort by ascend and keep the top two
                            [~,idxS] = sort(distLast,'ascend');
                            if length(idxS) >= 2
                            spotTrace2(s,t,1:2,2) = newPts2(idxS(2),:)';
                            spotTrace2(s,t,5:end,2) = curDat2(idxS(2),:)';
                            end
                            % may want to toss second doublet if first doublet is nan
                            dt = sqrt((newTime2-refTime2)/zDepth);
                            if minDistLast < maxStep2*dt && distMov(idMinDist) < maxDistFromMovAve2
                                spotsPerStep(s,t,2) = size(newPts2,1);
                                refPt2 =  newPts2(idMinDist,:);
                                refTime2 = newTime2;
                                spotTrace2(s,t,1:2,1) = refPt2;
                                spotTrace2(s,t,5:end,1) = curDat2(idMinDist,:);
                                maxStep2 = pars.maxStep;
                            end
                        end
                    end
                end
            end
        catch er
            warning(er.getReport);
            disp('place debug here');
        end
        
        
        if pars.movAveFig > 0
        fMovAve = figure(pars.movAveFig); hold on;
        fMovAve.Position = [0 0 1200 1200];
            for z=1:zDepth
                xy = squeeze(spotTrace1(s,z:zDepth:end,1:2,1));
                plot(xy(:,1),xy(:,2),'.','color',[.5 0 0],'MarkerSize',2);
                xy = squeeze(spotTrace2(s,z:zDepth:end,1:2,1));
                plot(xy(:,1),xy(:,2),'b.','MarkerSize',2);
            end
            if pars.saveFigs
                SaveFigure(fMovAve,'name', ['s',num2str(s),'MovingAve'],'formats',{'png'} ,'overwrite' ,pars.overwrite);
            end
        end
        
        npp=pars.npp;
        if pars.traceFig > 0
            fTrace = figure(pars.traceFig); clf;
            fTrace.Position = [0,0,600,1200];
            for z=1:zDepth
                subplot(zDepth,1,z);
                tim = z:zDepth:nFrames;
                d12 = abs(spotTrace1(s,z:zDepth:end,1,1) - spotTrace2(s,z:zDepth:end,1,1));
                dd12 = nanmedian(abs(diff(d12(~isnan(d12)))));
                x = squeeze(spotTrace1(s,z:zDepth:end,1,1));
                plot(tim(~isnan(x)),x(~isnan(x)),'.-'); hold on;
                dx1 = npp*nanmedian(abs(diff(x(~isnan(x)))));
                nx1 = sum(~isnan(x));
                x = squeeze(spotTrace2(s,z:zDepth:end,1,1));
                plot(tim(~isnan(x)),x(~isnan(x)),'.-'); hold on;
                dx2 = npp*nanmedian(abs(diff(x(~isnan(x)))));
                nx2 = sum(~isnan(x));
                title(['s=',num2str(s), '  d: ',num2str(npp*nanmedian(d12),3), 'nm.  \delta d:',num2str(npp*dd12,3),...
                    'nm.    x-step: ',num2str(dx1,3),'nm.   ',num2str(dx2,3),'nm.',...
                    '  n1=',num2str(nx1),'  n2=',num2str(nx2)]);
                xlim([0,nFrames]);
            end
            if pars.saveFigs
                SaveFigure(fTrace,'name', ['s',num2str(s),'TraceFig'],'formats',{'png'},'overwrite' ,pars.overwrite);
            end


            if pars.pause
                 disp('change pars if needed')
            end
        end
    end

else
%% parallel

    % split into a cell array for more efficient parallel processing
    pars.tLapseMov = 0; % avoid replicating a large array in passing pars
    
    sFits1 = cell(nSpots,1);
    sFits2 = cell(nSpots,1);
    spotTraces1 = cell(nSpots,1);
    spotTraces2 = cell(nSpots,1);
    spotsPerSteps = cell(nSpots,1);
    moveAveAllInts1 = cell(nSpots,1);
    moveAveAllInts2 = cell(nSpots,1);
    for s = pars.selSpots % s=8
        %%
        x0 = xy1p(s,1) - pars.windowR;
        x1 = xy1p(s,1) + pars.windowR;
        y0 = xy1p(s,2) - pars.windowR;
        y1 = xy1p(s,2) + pars.windowR;
        selPts = fits1.x > x0 & fits1.x < x1 & fits1.y > y0 & fits1.y < y1;
        sFits1{s} = fits1(selPts,:);
        selPts = fits2c.x > x0 & fits2c.x < x1 & fits2c.y > y0 & fits2c.y < y1;
        sFits2{s} = fits2c(selPts,:);
        moveAveAllInts1{s} = squeeze(moveAveAllInt1(s,:,:));
        moveAveAllInts2{s} = squeeze(moveAveAllInt2(s,:,:));
    end
    
     for s = pars.selSpots % s=8  parfor this
        %%
        x0 = xy1p(s,1) - pars.windowR;
        x1 = xy1p(s,1) + pars.windowR;
        y0 = xy1p(s,2) - pars.windowR;
        y1 = xy1p(s,2) + pars.windowR;
        selPts = fits1.x > x0 & fits1.x < x1 & fits1.y > y0 & fits1.y < y1;
        sFit1 = fits1(selPts,:);
        selPts = fits2c.x > x0 & fits2c.x < x1 & fits2c.y > y0 & fits2c.y < y1;
        sFit2 = fits2c(selPts,:);
        
        % just plotting
        if pars.spotMapFig > 0
            figure(pars.spotMapFig); clf; Ncolor(tLapseMov);  axis image;
            hold on; plot(xy1p(:,1),xy1p(:,2),'yo','MarkerSize',20);
            hold on; plot(xy2p(:,1),xy2p(:,2),'gs','MarkerSize',20);
            sNum = cellstr(num2str((1:size(xy1p,1))'));
            text(xy1p(:,1)+25,xy1p(:,2),sNum,'color','w'); hold on;
            xlim([x0,x1]); ylim([y0,y1]);
            title(['spot ',num2str(s)]);
            scatter(sFit1.x,sFit1.y,[],sFit1.frame,"filled"); colormap(jet); hold on;
            plot(sFit1.x,sFit1.y,'m.','MarkerSize',10); colormap(jet); hold on;
            scatter(sFit2.x,sFit2.y,[],sFit2.frame,"filled"); colormap(jet); hold on;
            plot(sFit2.x,sFit2.y,'c.','MarkerSize',10); colormap(jet); hold on;
        end
        
    
    
        movAve1 = nan(tStps,2);
        movAve2 = nan(tStps,2); 
        for t=1:tStps % t=9
            % to avoid the occassional stray-background point in an early frame
            % seeding the trace in an wrong position (which typically runs a short
            % and terminal path), we average the movie in blocks 
            % 
            % block averaging
            frame0 = (t-1)*tBlock; % indexed from 0
            frame1 = min(t*tBlock,nFrames); % 1-block more or last-frame,  
            sx = moveAveAllInts1{s}(t,1);
            sy = moveAveAllInts1{s}(t,2);
            w = pars.maxDistFromRoughAve;
    
            % get frames in time window and within distance from rough moving average  
            frameT = sFit1.frame >= frame0 & sFit1.frame < frame1 ;
            frameXY = (sFit1.x > (sx-w)) & (sFit1.x < (sx+w)) & (sFit1.y > (sy-w)) & (sFit1.y < (sy+w));
            selFrames1 = frameT & frameXY;
            if sum(selFrames1) > pars.minObsPerAve2
                fitsB1 = sFit1(selFrames1,:);
                movAve1(t,:) = [mean(fitsB1.x),mean(fitsB1.y)]; % time local mean 
            end
            % chn2 
            sx = moveAveAllInts2{s}(t,1);
            sy = moveAveAllInts2{s}(t,2);
            w = pars.maxDistFromRoughAve;
            selFrames2 = sFit2.frame >= frame0 & sFit2.frame < frame1  & ...
                 sFit2.x > sx - w & sFit2.x < sx + w & sFit2.y > sy - w & sFit2.y < sy+ w;
            if sum(selFrames2)  > pars.minObsPerAve
                fitsB2 = sFit2(selFrames2,:);
                movAve2(t,:) = [mean(fitsB2.x),mean(fitsB2.y)];  % local mean  
            end
        end
    
        % clean up moving average
        % use moving average to define the expected stepsize 
        % remove outliers in the moving average
        if pars.movAveFig > 0
             figure(pars.movAveFig);  % figure(15); 
             clf; Ncolor(tLapseMov);  axis image; hold on;
             scatter(movAve1(:,1),movAve1(:,2),[],(1:tStps)','filled'); colormap('jet'); hold on;
             plot(movAve1(:,1),movAve1(:,2),'.','color',[.5 0 0],'MarkerSize',10);
             hold on; scatter(movAve2(:,1),movAve2(:,2),[],(1:tStps)','filled'); colormap('jet')
             plot(movAve2(:,1),movAve2(:,2),'.','color',[0 0 1],'MarkerSize',10);
                    xlim([x0,x1]); ylim([y0,y1]);
                title(['all ave. spot ',num2str(s)]); pause(.1);
        end
    
        moveAveClean1 = RemoveJumps(movAve1,'maxDiff',pars.movAveStepMaxFold,'localRegion',pars.movAveLocalStep);
        hasDat1 = ~isnan(moveAveClean1(:,1));
        aveStep1 = sqrt(sum(diff(moveAveClean1(hasDat1,:)).^2,2));
        medAveStep1 = nanmedian(aveStep1); %#ok<*NANMEDIAN>
        moveAveClean2 = RemoveJumps(movAve2,'maxDiff',pars.movAveStepMaxFold,'localRegion',pars.movAveLocalStep);
        hasDat2 = ~isnan(moveAveClean2(:,1));
        aveStep2 = sqrt(sum(diff(moveAveClean2(hasDat2,:)).^2,2));
        medAveStep2 = nanmedian(aveStep2); 
       
        if pars.movAveFig > 0
             figure(pars.movAveFig); clf; Ncolor(tLapseMov);  axis image; hold on;
             x = moveAveClean1(:,1); y = moveAveClean1(:,2);
             plot(x(~isnan(x)),y(~isnan(y)),'-','color',[.5 0 0]); hold on;
             x = moveAveClean2(:,1); y = moveAveClean2(:,2);
             plot(x(~isnan(x)),y(~isnan(y)),'-','color',[0 0 1]); hold on;
             scatter(moveAveClean1(:,1),moveAveClean1(:,2),[],(1:tStps)','filled'); colormap('jet'); hold on;
             plot(moveAveClean1(:,1),moveAveClean1(:,2),'.','color',[.5 0 0],'MarkerSize',10);
             hold on; scatter(moveAveClean2(:,1),moveAveClean2(:,2),[],(1:tStps)','filled'); colormap('jet')
             plot(moveAveClean2(:,1),moveAveClean2(:,2),'.','color',[0 0 1],'MarkerSize',10);
                    xlim([x0,x1]); ylim([y0,y1]);
                title(['clean ave. spot ',num2str(s)]);
        end
    
    
        %% compute traces
        if pars.dynFig > 0
            im = tLapseMov(y0:y1,x0:x1,:);
        end
        
        try
            valid1 = moveAveClean1(~isnan(moveAveClean1(:,1)),:);
            valid2 = moveAveClean2(~isnan(moveAveClean2(:,1)),:);
            refPt1 = valid1(1,:); % start with first non-nan
            refPt2 = valid2(1,:);
            spotTraces1{s} = nan(nFrames,8,2); % (t,1:2,2
            spotTraces2{s} = nan(nFrames,8,2);% (t,1:2,2)
            maxStep1 = pars.maxStepFrame1;
            maxStep2 = pars.maxStepFrame1;
            maxDistFromMovAve1 = pars.maxDistFromMovingAve*medAveStep1;
            maxDistFromMovAve2 = pars.maxDistFromMovingAve*medAveStep2;
            for z=1:zDepth
                refTime1 = find(~isnan(moveAveClean1(:,1)),1);
                refTime2 = find(~isnan(moveAveClean2(:,1)),1);

                for t=z:zDepth:nFrames 
                    % accept the closest spot to the current reference that is also
                    % within the required distance from the reference centroid. 
                    tt = floor((t-1)/tBlock)+1;
                    xy1 = movAve1(tt,:);
                    xy2 = movAve2(tt,:);
                    if ~isnan(xy1)
                        currZstack = sFit1.frame == t-1; 
                        fr1 = sFit1(currZstack,:);
                        if ~isempty(fr1)
                            newTime1 = t;
                            newPts = [fr1.x,fr1.y];
                            curDat1 = [fr1.height,fr1.background,fr1.xsigma,fr1.significance];
                            distLast = sqrt(sum((newPts - repmat(refPt1,size(newPts,1),1)).^2,2));
                            distMov = sqrt(sum((newPts - repmat(xy1,size(newPts,1),1)).^2,2));
                            [minDistLast,idMinDist] = min(distLast);  % only tracing a single path at this point so min distance is okay, don't need MatchPairs 
                            % sort by ascend and keep the top two
                            [~,idxS] = sort(distLast,'ascend');
                            if length(idxS) >= 2
                            spotTraces1{s}(t,1:2,2) = newPts(idxS(2),:)';
                            spotTraces1{s}(t,5:end,2) = curDat1(idxS(2),:)';
                            end
                            % may want to toss second doublet if first doublet is nan  
                            dt = sqrt((newTime1-refTime1)/zDepth);  % allow larger steps for larger time lags  
                            if minDistLast < maxStep1*dt && distMov(idMinDist) < maxDistFromMovAve1
                                if pars.dynFig>0 % dynamically and sequentially add spots  (really just for troubleshooting)
                                    figure(pars.dynFig); clf;
                                    Ncolor(im);  axis image; hold on;
                                    plot(refPt1(1)-x0,refPt1(2)-y0,'r+');
                                    plot(newPts(:,1)-x0,newPts(:,2)-y0,'w.');
                                    plot(newPts(idMinDist,1)-x0,newPts(idMinDist,1)-y0,'ro'); 
                                    title('chn1'); pause(.1); 
                                end
                                spotsPerSteps{s}(t,1) = size(newPts,1);
                                refPt1 =  newPts(idMinDist,:);
                                refTime1 = newTime1;
                                spotTraces1{s}(t,1:2,1) = refPt1;
                                spotTraces1{s}(t,5:end,1) = curDat1(idMinDist,:);
                                maxStep1 = pars.maxStep;
                            end
                        end   
                    end
                    if ~isnan(xy2)
                        currZstack = sFit2.frame == t-1; % 
                        fr2 = sFit2(currZstack,:);
                        if ~isempty(fr2)
                            newTime2 = t;
                            newPts2 = [fr2.x,fr2.y]; % z,t, h,bkd,sigma
                            curDat2 = [fr2.height,fr2.background,fr2.xsigma,fr2.significance];
                            distLast = sqrt(sum((newPts2 - repmat(refPt2,size(newPts2,1),1)).^2,2));
                            distMov = sqrt(sum((newPts2 - repmat(xy2,size(newPts2,1),1)).^2,2))  ;
                            % select
                            %  currently: closest point that is within the step size
                            %  tolerance and that is within the appropriate distance from
                            %  the moving average.  
                            [minDistLast,idMinDist] = min(distLast);% sort by ascend and keep the top two
                            [~,idxS] = sort(distLast,'ascend');
                            if length(idxS) >= 2
                            spotTraces2{s}(t,1:2,2) = newPts2(idxS(2),:)';
                            spotTraces2{s}(t,5:end,2) = curDat2(idxS(2),:)';
                            end
                            % may want to toss second doublet if first doublet is nan
                            dt = sqrt((newTime2-refTime2)/zDepth);
                            if minDistLast < maxStep2*dt && distMov(idMinDist) < maxDistFromMovAve2
                                spotsPerSteps{s}(t,2) = size(newPts2,1);
                                refPt2 =  newPts2(idMinDist,:);
                                refTime2 = newTime2;
                                spotTraces2{s}(t,1:2,1) = refPt2;
                                spotTraces2{s}(t,5:end,1) = curDat2(idMinDist,:);
                                maxStep2 = pars.maxStep;
                            end
                        end
                    end
                end
            end
        catch er
            warning(er.getReport);
            disp('place debug here');
        end
        
        
        if pars.movAveFig > 0
        fMovAve = figure(pars.movAveFig); hold on;
        fMovAve.Position = [0 0 1200 1200];
            for z=1:zDepth
                xy = squeeze(spotTraces1{s}(z:zDepth:end,1:2,1));
                plot(xy(:,1),xy(:,2),'.','color',[.5 0 0],'MarkerSize',2);
                xy = squeeze(spotTraces2{s}(z:zDepth:end,1:2,1));
                plot(xy(:,1),xy(:,2),'b.','MarkerSize',2);
            end
            if pars.saveFigs
                SaveFigure(fMovAve,'name', ['s',num2str(s),'MovingAve'],'formats',{'png'} ,'overwrite' ,pars.overwrite);
            end
        end
        
        npp=pars.npp;
        if pars.traceFig > 0
            fTrace = figure(pars.traceFig); clf;
            fTrace.Position = [0,0,600,1200];
            for z=1:zDepth
                subplot(zDepth,1,z);
                tim = z:zDepth:nFrames;
                d12 = abs(spotTraces1{s}(z:zDepth:end,1,1) - spotTraces2{s}(z:zDepth:end,1,1));
                dd12 = nanmedian(abs(diff(d12(~isnan(d12)))));
                x = squeeze(spotTraces1{s}(z:zDepth:end,1,1));
                plot(tim(~isnan(x)),x(~isnan(x)),'.-'); hold on;
                dx1 = npp*nanmedian(abs(diff(x(~isnan(x)))));
                nx1 = sum(~isnan(x));
                x = squeeze(spotTraces2{s}(z:zDepth:end,1,1));
                plot(tim(~isnan(x)),x(~isnan(x)),'.-'); hold on;
                dx2 = npp*nanmedian(abs(diff(x(~isnan(x)))));
                nx2 = sum(~isnan(x));
                title(['s=',num2str(s), '  d: ',num2str(npp*nanmedian(d12),3), 'nm.  \delta d:',num2str(npp*dd12,3),...
                    'nm.    x-step: ',num2str(dx1,3),'nm.   ',num2str(dx2,3),'nm.',...
                    '  n1=',num2str(nx1),'  n2=',num2str(nx2)]);
                xlim([0,nFrames]);
            end
            if pars.saveFigs
                SaveFigure(fTrace,'name', ['s',num2str(s),'TraceFig'],'formats',{'png'},'overwrite' ,pars.overwrite);
            end


            if pars.pause
                 disp('change pars if needed')
            end
        end
    end



    try
        spotTrace1 = cat(4,spotTraces1{:});
        spotTrace1 = permute(spotTrace1,[4,1,2,3]); % swap dims 
        spotTrace2 = cat(4,spotTraces2{:});
        spotTrace2 = permute(spotTrace2,[4,1,2,3]); % swap dims 
        spotsPerStep = cat(3,spotsPerSteps{:});
        spotsPerStep = permute(spotsPerStep,[3,1,2]); % swap dims
    catch er
        warning(er.getReport);
        disp('place debug here');
    end

end
