function [finalOrder,overlapMatrix] = SortFOVsByOverlap(mosaicSize, frameSize, stageXY, varargin)
% Uses BFS to determine an order for image insertion into mosaic that
% ensures correctness. Prioritizes FOVs based on number of overlapping
% other FOVs

%% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'resize','nonnegative',0.1}; % scale down images to calculate overlap; loses some accuracy, but much faster
defaults(end+1,:) = {'data','cell',{}}; % determine nonempty pixels in data; 
defaults(end+1,:) = {'dataThreshold','fraction',.5};  % display intermediate figures
defaults(end+1,:) = {'debug','boolean',false};  % display intermediate figures
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Calculate Image Processing Order
numIms = length(stageXY); 
% rewrote to use UL not center position. 
%   this saves stupid padding to avoid negative indexing

% Generate boolean images from position data
boolIms = cell(numIms,1);
w = frameSize(1) * pars.resize;
h = frameSize(1) * pars.resize;
resizedMosaicSize = ceil(mosaicSize * pars.resize);
totalImage = false(resizedMosaicSize);
for f = 1:numIms
    curX = stageXY(f,1) * pars.resize;
    curY = stageXY(f,2) * pars.resize;
%     rowStart = round(curY - h/2);
%     rowEnd   = round(curY + h/2);
%     colStart = round(curX - w/2);
%     colEnd   = round(curX + w/2);
    rowStart = ceil(curY);
    rowEnd   = ceil(curY + h);
    colStart = ceil(curX);
    colEnd   = ceil(curX + w);
    
    curImage = false(resizedMosaicSize);   % can probably do this without making every image represented by something the size of the whole mosaic. 
    if isempty(pars.data)
        curImage(rowStart:rowEnd, colStart:colEnd) = true;  
    else
        hasData = pars.data{f};
        h_i = rowEnd-rowStart+1;
        w_i = colEnd-colStart+1;
        hasData = double(imresize(hasData,[h_i,w_i]));
        cut = quantile(hasData(:),pars.dataThreshold);
        hasData(hasData < cut) = 0;
        hasData(hasData >= cut)= 1;
        curImage(rowStart:rowEnd, colStart:colEnd) = boolean(hasData);  
    end
    boolIms{f} = curImage;
    
    totalImage = totalImage | curImage;
end




% Calculate overlap
overlapMatrix = zeros(numIms,numIms);
for i = 1:numIms
    for j = 1:i
        if (i == j)
            continue
        end
        overlap = sum(sum(boolIms{i} & boolIms{j}));
        overlapMatrix(i,j) = overlap;
        overlapMatrix(j,i) = overlap;
    end
end

% Create graph representation; run breadth first search 
G = graph(overlapMatrix);
[~, nodeOrder] = sort(degree(G), 'descend');
G = reordernodes(G, nodeOrder);  % Start with image that overlaps with most other images (maybe start with most overlapping pixels instead?)
nodes = 1:numIms;
sortedNodes = nodes(nodeOrder);

v = bfsearch(G, 1, 'Restart', true);
finalOrder = sortedNodes(v);

if (pars.debug)
    figure(12); clf; plot(G);

    figure(13); clf; imagesc(totalImage); hold on;
    for i = 1:numIms
        f = sortedNodes(v(i));  % final orderering (accounts for previous sorting of nodes)
        curX = stageXY(f,1) * pars.resize;
        curY = stageXY(f,2) * pars.resize;
        rectangle('Position',[curX curY w h],'EdgeColor','y');
        text(round(curX), round(curY - h/4), num2str(f) + "->" ,'color','r','FontSize',14);
        text(round(curX), round(curY - h/4), "        " + num2str(i),'color','g','FontSize',14);
%         rectangle('Position',[curX-w/2 curY-h/2 w h],'EdgeColor','y');
%         text(round(curX - w/2), round(curY - h/4), num2str(f) + "->" ,'color','r','FontSize',14);
%         text(round(curX - w/2), round(curY - h/4), "        " + num2str(i),'color','g','FontSize',14);
    end
end
