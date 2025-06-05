function [pntChain,idx] = LinkLiveTable(fitStruct,varargin)
% pntTable (or structure)
% x,y,height,background,width,frame
% 
% outputs
% pntChain n-chains x 2-dim x t-frames matrix of traces
%     has nans for coordinates where a given chain was not observed. 
% 

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'selFrames','float',0}; % select frames.  0 for all
defaults(end+1,:) = {'seedPoints','float',[0,0]}; % pass seed points to use in searching for data (if zero, it will use the first frame's data as the seed. any point not visible in frame 0 will be lost) 
defaults(end+1,:) = {'maxStep','positive',1}; % 
defaults(end+1,:) = {'maxDistToSeed','positive',15}; % 
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

tic
nPts = size(xy0,1); % length(xs{1});
idx = cell(nFrames,1); 
pntChain = nan(nPts,2,nFrames);
for f=2:nFrames % -1
    try

    fn = selFrames(f);
    if pars.verbose
        if rem(f,1000)==0
            disp([num2str(f/nFrames*100,3),'% complete'])
        end
    end
    if f==2    
        xy1 = xy0; 
        xy2 = [xs{fn},ys{fn}];      
        [matched,cost] = MatchPoints(xy2,xy1,'maxDist',pars.maxDistToSeed); % s 
        m =matched(cost<pars.maxDistToSeed,:);      
        nMatched = size(m,1);
        pntChain(1:nMatched,1:2,f-1) = xy1(m(:,2),1:2); % set point 1 to the start point (might want to toss this out in the end).
        pntChain(1:nMatched,1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point
        % pntChain(nMatched+1:end,:,:) = []; % anyone who's first point didn't get matched gets tossed.  
    else
        xy1 = pntChain(:,1:2,f-1);
        xy2 = [xs{fn},ys{fn}];
        % if missed in last frame, use previous frame
        missed = isnan(xy1(:,1));
        ff = f-2;
        while any(missed) && ff >=1
            xy1(missed,:) = pntChain(missed,1:2,ff);
            missed = isnan(xy1(:,1));
            ff=ff-1;
        end
        xy1n = xy1;  % this is the new xy1, which hassing missing values back filled from the last-spot seen in the trace.  
        %  xy1n(isnan(xy1(:,1)),:) = 0; % if any are still missing, replace with 0
        xy1n(isnan(xy1(:,1)),:) = xy0(isnan(xy1(:,1)),:); % if still missing, replace with seed.  
        [matched,cost] = MatchPoints(xy2,xy1n,'maxDist',pars.maxStep);
        m =matched(cost<pars.maxStep,:);      
        pntChain(m(:,2),1:2,f-1) = xy1(m(:,2),1:2);
        pntChain(m(:,2),1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point
       
        idx{f} = m;
        % data(m(:,2),all_data_types,frameNum) = dataInFrame_f(m(:,1),all_data_types) 
    end
    if showPlot && f>2
        figure(pars.figHandle);
        x1 = squeeze(pntChain(:,1,f-1))';
        x2 = squeeze(pntChain(:,1,f))';
        y1 = squeeze(pntChain(:,2,f-1))';
        y2 = squeeze(pntChain(:,2,f))';
        plot([x1;x2],[y1;y2],'color',cmap(f,:)); hold on;
        plot([x1;x2],[y1;y2],pars.symbol,'color',cmap(f,:));
        plot(xy2(:,1),xy2(:,2),'.','color',cmap(f,:));
    end
    catch er
        disp(['error on ',num2str(f)])
        disp(er.getReport);
        disp('stop here')
    end
end
toc
disp('all points linked')

% Pull the rest of the data into the point chains
pntArray = nan(size(pntChain,1),nFrames,6);  % x,y,z,t,h,b
% data(m(:,2),frameNum,all_data_types) = dataInFrame_f(m(:,1),all_data_types) 
for f=2:nFrames  % f = 5
    m = idx{f};
    if ~isempty(m)
        pntArray(m(:,2),f,1) = fitStruct(f).x(m(:,1)); % that's compact but not super legible;  
        pntArray(m(:,2),f,2) = fitStruct(f).y(m(:,1)); % the pattern is helping with the reading
        pntArray(m(:,2),f,5) = fitStruct(f).height(m(:,1)); % spot brightness  (will need this later for calculating z)
        pntArray(m(:,2),f,4) = fitStruct(f).frame*ones(length(m(:,1)),1); % time (by frame)
        pntArray(m(:,2),f,6) = fitStruct(f).background(m(:,1)); % spot brightness  (will need this later for calculating z)
    end
end

