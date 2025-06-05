function [pntArray,idx,pntChain] = LinkLiveZstack(fitStruct,varargin)
% 
% fitStruct.  has elements
% x,y,height,background,width,frame
% 
%% outputs
%  (should make pntArray the primary output)
% pntChain n-chains x 2-dim x t-frames matrix of traces
%     has nans for coordinates where a given chain was not observed. 
% 
%% update history
% developed from LinkLiveStruct()
%  forked to process frames based on their position in the z-scan cycle

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'zDepth','integer',5};
defaults(end+1,:) = {'selFrames','float',0}; % select frames.  0 for all
defaults(end+1,:) = {'seedPoints','float',[0,0]}; % pass seed points to use in searching for data (if zero, it will use the first frame's data as the seed. any point not visible in frame 0 will be lost) 
defaults(end+1,:) = {'maxStep','positive',1}; % 
defaults(end+1,:) = {'maxDistToSeed','positive',15}; % 
defaults(end+1,:) = {'uniqueMatch','boolean','false'}; 
defaults(end+1,:) = {'symbol','string','.'}; 
defaults(end+1,:) = {'figHandle','freeType',[]}; % 
pars = ParseVariableArguments(varargin,defaults,mfilename);

% in case there are any frames with no data, we still want to make sure we
% don't accidently drop the empty frames, for the sake of caution we pass
% around an empty cell for any empty frames. 
if istable(fitStruct)
    xs =cell(max(fitStruct.frame),1);
    xs(fitStruct.frame) = fitStruct.x;
    ys =cell(max(fitStruct.frame),1);
    ys(fitStruct.frame) = fitStruct.y;
elseif isstruct(fitStruct)
    xs =cell(max([fitStruct.frame]),1);
    xs([fitStruct.frame]) = {fitStruct.x};
    ys =cell(max([fitStruct.frame]),1);
    ys([fitStruct.frame]) = {fitStruct.y};
else
    error('input must be either a table or a structure')
end


if pars.selFrames <= 0 
    nFrames = length(xs);
    selFrames = 1:nFrames;
else
    selFrames = pars.selFrames;
end
nFrames = length(selFrames);

% show frames if not to many
if nFrames < 100
    showPlot = true;
    cmap = hsv(nFrames); 
    if isempty(pars.figHandle)
        pars.figHandle = figure(3); 
    end
else
    showPlot = false;
end

if sum(pars.seedPoints) == 0
    xy0 = [xs{1},ys{1}];
else
    xy0 = pars.seedPoints;
end
% 
% xya = [xs{1},ys{1}];
% figure(10); clf; plot(xya(:,1),xya(:,2),'r+'); hold on;
% plot(pars.seedPoints(:,1),pars.seedPoints(:,2),'bo');
if pars.verbose
    disp('linking data...')
end

tic
nPts = size(xy0,1); % length(xs{1});
idx = cell(nFrames,1); 
pntChain = nan(nPts,2,nFrames+1);
for f=1:nFrames % -1
    try
    fn = selFrames(f);
    if pars.verbose
        if rem(f,1000)==0
            disp(['linking ',num2str(f/nFrames*100,3),'% complete'])
        end
    end
    if f==1    
        xy1 = xy0;   % seed point 
        xy2 = [xs{fn},ys{fn}];   % first data frame  
        [matched,cost] = MatchPoints(xy2,xy1,'maxDist',pars.maxDistToSeed); %
        m =matched(cost<pars.maxDistToSeed,:);  
        % this maximizes unique mapping of frame1 pts to the seed points.
        %   given the potentially large gaps between a spots
        %   motion-centroid and its first frame, this is probably a good
        %   idea.  The unique mapping makes sure the sum of the
        %   distance between all matches is minimum, rather than assigning
        %   each nearest neighbors.
        %     
        nMatched = size(m,1);
        pntChain(m(:,2),1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point
    else
        % then we take these points from frame 1 and start matching to
        % frame 2.  The original approach is to again use MatchPoints
        % unique mapping.  
        xy1 = pntChain(:,1:2,f-1); % this is 1 frame ahead of the data, because we kept the seed
        xy2 = [xs{fn},ys{fn}];
        % if missed in last frame, use previous frame
        missed = isnan(xy1(:,1));
        ff = f-2;
        while any(missed) && ff >=1
            xy1(missed,:) = pntChain(missed,1:2,ff);
            missed = isnan(xy1(:,1));
            ff=ff-1;
        end
        xy1n = xy1;  % this is the new xy1, which has missing values back filled from the last-spot seen in the trace.  
        stillMissing = isnan(xy1(:,1));
        xy1n(stillMissing,:) = xy0(stillMissing,:); % if still missing, replace with seed.  
       

        % if still missing, keep using the dist
        [matched,cost] = MatchPoints(xy2, xy1n(stillMissing,:),'maxDist',pars.maxDistToSeed);
        m =matched(cost<pars.maxDistToSeed,:);  
        idxMissing = find(stillMissing);
        xy1n(idxMissing(m(:,2)),:) = xy2(m(:,1),:); % update still missing and matched
       

        if pars.uniqueMatch
            [matched,cost] = MatchPoints(xy2,xy1n,'maxDist',pars.maxStep);
            m =matched(cost<pars.maxStep,:);  
            % this uniquely assigns points
            %   maybe we should just assign ne
        else
            [matched,cost] = knnsearch(xy2,xy1n); % knnsearch(X,Y) - nearest in X for each query point in Y   
            m2 = matched(cost<pars.maxStep); % values in xy2
            m1 = (1:size(xy1n,1))';
            m1(cost>=pars.maxStep) = [];
            m=[m2,m1];  % (m1 is original)  %  This structure already implies a 1-1 mapping. 
            %   I suppose we can add duplicates to m2 and m1 to allow non 1-1 but 
     
           %  xy2n(:,1:2) = xy2(m2,1:2); % points in 2
        end
        
        if ~isempty(m)
            % pntChain(m(:,2),1:2,f-1) = xy1(m(:,2),1:2);
           % shouldn't be necessary? (doesn't seem to effect outcome either)  
            pntChain(m(:,2),1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point
        end
    end
    idx{f} = m;
        % data(m(:,2),all_data_types,frameNum) = dataInFrame_f(m(:,1),all_data_types) 
    if showPlot && f>2
        figure(pars.figHandle); clf;
        imagesc(zeros(2048,2048,'uint8')); hold on;
        x1 = squeeze(pntChain(:,1,f-1))';
        x2 = squeeze(pntChain(:,1,f))';
        y1 = squeeze(pntChain(:,2,f-1))';
        y2 = squeeze(pntChain(:,2,f))';
        plot([x1;x2],[y1;y2],'color',cmap(f,:)); hold on;
        plot([x1;x2],[y1;y2],pars.symbol,'color',cmap(f,:));
        plot(xy2(:,1),xy2(:,2),'.','color',cmap(f,:));
        plot(xy0(:,1),xy0(:,2),'k+');
    end
    catch er
        disp(['error on ',num2str(f)])
        disp(er.getReport);
        disp('stop here')
    end
end
toc
disp('all points linked')

% should check that there is extra data columns first?

% Pull the rest of the data into the point chains
pntArray = nan(size(pntChain,1),nFrames,7);  % x,y,z,t,h,b 
% (maybe this should be a structure for each spot?)  
% (I kinda like matrices though
% data(m(:,2),frameNum,all_data_types) = dataInFrame_f(m(:,1),all_data_types) 
for f=1:nFrames  % f = 5
    m = idx{f};
    if ~isempty(m)
        pntArray(m(:,2),f,1) = fitStruct(f).x(m(:,1)); % that's compact but not super legible;  
        pntArray(m(:,2),f,2) = fitStruct(f).y(m(:,1)); % the pattern is helping with the reading
        pntArray(m(:,2),f,5) = fitStruct(f).height(m(:,1)); % spot brightness  (will need this later for calculating z)
        pntArray(m(:,2),f,4) = fitStruct(f).frame*ones(length(m(:,1)),1); % time (by frame)
        pntArray(m(:,2),f,6) = fitStruct(f).background(m(:,1)); % spot brightness  (will need this later for calculating z)
        pntArray(m(:,2),f,7) = fitStruct(f).sigma(m(:,1)); % spot brightness  (can use this later for calculating z)
    end
end

