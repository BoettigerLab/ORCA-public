function [outTable,pntChain] = SimpleSpotLinking(inTable,xy0,varargin)
% takes data from a table with entries x, y, and frame, and a list of seed
% loci, xy0.  For each row in the xy0, we assign up to one (x,y) pair per
% frame, thus creating a chain of spots.
% 
% Inputs 
% - inTable.  a matlab table.
%     Must have Columns labeled "x", "y" and "frame". May have
%     may have arbitrary amounts of additional columns if desired.
% - xy0 - an N x 2 column matrix
%     Lists the seed locations to use in constructing linked trajectories
% 
% Outputs - outTable
%    Identical column labels as the input table, but now all rows contain a
%    traceID column, which uniquely links the spots through the frames.
% 
% Optional Inputs  (see below)
% 
% 
% 
% Optional Output - pntChain  
%      - data tensor, nPts by 2 by nFrames+1; where nPts is the number of
% rows in xy0, nFrames is the number of frames in the datatable fits1, and
% the second dimension is the x,y coordinates.
% example
%  pntChain(3,:,1:nFrames) is a 2xnFrames matrix of recording the xy 
%  position of spot number 3 for all frames
% 
% 
%  Notes: 
%    should we start tracing from frame 1?  Often this has more noise.
%    Maybe we should start tracing backwards from last frame? Sometimes
%    this has less signal, or some spots have bleached or drifted away. 


defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'selFrames','float',0}; % select frames.  0 for all
defaults(end+1,:) = {'showSeed','boolean',false}; % autoselect threshold downsampling for grouping localizations in z. 
defaults(end+1,:) = {'maxStep','positive',3}; % max 2D step distance allowed between consecutive frames (last detected frame?)  [if you didn't detect it, should you allow the distance to grow?]   
defaults(end+1,:) = {'maxDistToSeed','positive',10}; %  max distance allowed relative to the seed point for the initial trace 
pars = ParseVariableArguments(varargin,defaults,mfilename);


showPlot = false; 
  
%%
if pars.selFrames == 0
    maxFrame = max(inTable.frame);
    minFrame = min(inTable.frame);
    selFrames = 1:(maxFrame-minFrame+1);
else
    selFrames = pars.selFrames;
end
nFrames = length(selFrames);
cmap = hsv(nFrames);

tic
nPts = size(xy0,1); % length(xs{1});
idx = cell(nFrames,1); 
pntChain = nan(nPts,2,nFrames+1);
nData = length(inTable.Properties.VariableNames);
dataChain = nan(nPts,nData,nFrames+1);
noChains = true;
for f=1:nFrames % -1
    % disp(f)
    try
        fn = selFrames(f);
        if pars.verbose
            if rem(f,1000)==0
                disp(['linking ',num2str(f/nFrames*100,3),'% complete'])
            end
        end

        isF = inTable.frame==fn;
        xy2 = [inTable.x(isF),inTable.y(isF)];
        % isColumnNotXY = ~(strcmp(fits1.Properties.VariableNames,'x') | strcmp(fits1.Properties.VariableNames,'y'));
        tableData = inTable(isF,:);

        if noChains    
            xy1 = xy0;   % seed point 
            

            if ~isempty(xy2)
                [matched,cost] = MatchPoints(xy2,xy1,'maxDist',pars.maxDistToSeed); % s 
                m =matched(cost<pars.maxDistToSeed,:);  
            else
                m = [];
            end
            if ~isempty(m)
                nMatched = size(m,1);
                pntChain(m(:,2),1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point
                dataChain(m(:,2),:,f)= tableData{m(:,1),:};
                noChains = false;
            else
                noChains = true;
            end
        else
            xy1 = pntChain(:,1:2,f-1); % this is 1 frame ahead of the data, because we kept the seed
            % if missed in last frame, use previous frame
            missed = isnan(xy1(:,1));
            ff = f-2;
            while any(missed) && (ff >=1)
                xy1(missed,:) = pntChain(missed,1:2,ff);
                missed = isnan(xy1(:,1));
                ff=ff-1;
            end
            xy1n = xy1;  % this is the new xy1, which hassing missing values back filled from the last-spot seen in the trace.  
            %  xy1n(isnan(xy1(:,1)),:) = 0; % if any are still missing, replace with 0
            xy1n(isnan(xy1(:,1)),:) = xy0(isnan(xy1(:,1)),:); % if still missing, replace with seed.  
            [matched,cost] = MatchPoints(xy2,xy1n,'maxDist',pars.maxStep);
            m =matched(cost<pars.maxStep,:);      
            if ~isempty(m) 
                pntChain(m(:,2),1:2,f-1) = xy1(m(:,2),1:2); %  --- this is just needed for propigating the links  
                pntChain(m(:,2),1:2,f) = xy2(m(:,1),1:2); % indexed in order of starting point

                dataChain(m(:,2),:,f)= tableData{m(:,1),:};
            end
        end
        idx{f} = m;
            % data(m(:,2),all_data_types,frameNum) = dataInFrame_f(m(:,1),all_data_types) 
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

if pars.verbose
    toc
    disp('all points linked')
end
% should check that there is extra data columns first?

%%

% figure(11); clf;
% imagesc( squeeze(pntChain(:,1,:)) - squeeze(nanmean(pntChain(:,1,:),3))  ) ; colorbar; clim([-4 4]);
% GetColorMap('RedWhiteBlueSat');

oldLabels = inTable.Properties.VariableNames;
newLabels = [oldLabels,'traceID'];
[nPts,nCols,nFrames] = size(dataChain);
datTable = cell(nPts*nFrames,1);
k=0;
for s=1:nPts
    for f=1:nFrames
        k=k+1;
        rowDat = [dataChain(s,:,f),s];
        if ~isnan(rowDat(1))
            datTable{k} = rowDat;
        end
    end
end
datTable = cat(1,datTable{:});
outTable = array2table(datTable); % ,'VariableNames',newLabels);
if ~isempty(outTable)
    outTable.Properties.VariableNames = newLabels;
end





