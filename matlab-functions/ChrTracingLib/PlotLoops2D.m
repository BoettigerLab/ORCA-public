function PlotLoops2D(mapin,varargin)


defaults = cell(0,3);
defaults(end+1,:) = {'theta','nonnegative',100};
defaults(end+1,:) = {'bondColor','colormap',.25*ones(1,3)};
pars = ParseVariableArguments(varargin, defaults, mfilename);

    nB = size(mapin,1);
    loopMap = mapin< pars.theta;
    loopMap = triu(loopMap,4);
   %  figure(1); clf; imagesc(loopMap)
    [loopx,loopy] = ind2sub([nB,nB],find(loopMap>0));
    
    currLoops = [loopx,loopy];
    G = digraph(currLoops(1,1),currLoops(1,2));
    nL=size(currLoops,1);
    for n=2:nL
        G = addedge(G,currLoops(n,1),currLoops(n,2));
    end
    for n=1:nB-1
        G = addedge(G,n,n+1);
    end
     %  subplot(K,2,2*(k-1)+1); cla;
    p = plot(G,'Layout','force','NodeLabel',{},'Marker','none','ArrowSize',0,'EdgeColor',pars.bondColor);
    x = p.XData; y = p.YData;
    nB = size(loopMap,1);
    nP = length(x);
    z = zeros(size(x)); col = 1:nP;
    surface([x;x],[y;y],[z;z],[col;col],'facecolor','no','edgecolor','interp','LineWidth',3);  
    hold on; scatter(x,y,[],col,'SizeData',3);
    GetColorMap('hsvCut',nB);
 