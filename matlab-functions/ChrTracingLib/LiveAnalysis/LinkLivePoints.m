function pntChain = LinkLivePoints(xs1,ys1,varargin)
% xs1  cell array of length nFrames of x coords
%   UPDATE this to use a table or a structure for more elegance. 

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'selFrames','float',0}; % select frames.  0 for all
defaults(end+1,:) = {'seedPoints','float',[0,0]}; % 
defaults(end+1,:) = {'maxStep','positive',1}; % 
defaults(end+1,:) = {'maxDistToSeed','positive',15}; % 
defaults(end+1,:) = {'symbol','string','.'}; 
defaults(end+1,:) = {'figHandle','freeType',[]}; % 
pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.selFrames <= 0 
    nFrames = length(xs1);
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
    xy0 = [xs1{1},ys1{1}];
else
    xy0 = pars.seedPoints;
end
% 
% xya = [xs1{1},ys1{1}];
% figure(10); clf; plot(xya(:,1),xya(:,2),'r+'); hold on;
% plot(pars.seedPoints(:,1),pars.seedPoints(:,2),'bo');

tic
nPts = size(xy0,1); % length(xs1{1});
pntChain = nan(nPts,2,nFrames);
for f=2:nFrames
    fn = selFrames(f);
    if pars.verbose
        if rem(f,1000)==0
            disp([num2str(f/nFrames*100,3),'% complete'])
        end
    end
    if f==2    
        xy1 = xy0; 
        xy2 = [xs1{fn},ys1{fn}];      
        [matched,cost] = MatchPoints(xy1,xy2);
        m =matched(cost<pars.maxDistToSeed,:);      
        nMatched = size(m,1);
        pntChain(1:nMatched,1:2,f-1) = xy1(m(:,1),1:2); % set point 1 to the start point (might want to toss this out in the end).
        pntChain(1:nMatched,1:2,f) = xy2(m(:,2),1:2); % indexed in order of starting point
        % pntChain(nMatched+1:end,:,:) = []; % anyone who's first point didn't get matched gets tossed.  
    else
        xy1 = pntChain(:,1:2,f-1);
        xy2 = [xs1{fn},ys1{fn}];
        % if missed in last frame, use previous frame
        missed = isnan(xy1(:,1));
        ff = f-2;
        while any(missed) && ff >=1
            xy1(missed,:) = pntChain(missed,1:2,ff);
            missed = isnan(xy1(:,1));
            ff=ff-1;
        end
        % if still missing, replace with seed.  
        xy1n = xy1; 
        %  xy1n(isnan(xy1(:,1)),:) = 0;
        xy1n(isnan(xy1(:,1)),:) = xy0(isnan(xy1(:,1)),:);  % 0
        [matched,cost] = MatchPoints(xy1n,xy2);
        m =matched(cost<pars.maxStep,:);      
        pntChain(m(:,1),1:2,f-1) = xy1(m(:,1),1:2);
        pntChain(m(:,1),1:2,f) = xy2(m(:,2),1:2); % indexed in order of starting point

       
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
end
toc
disp('all points linked')
