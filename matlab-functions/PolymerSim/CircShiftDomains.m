function [newNodeType,newDomainCoords,domainShift,nodeShift] = CircShiftDomains(nodeType,domainCoords,varargin) 
% newNodeType = CircShiftDomains(nodeType)
% computes a random circular shift of domains and returns the nodeType for
% every node.  The shift ensures that contiguous domains remain intact.  
% 
% newNodeType = CircShiftDomains(nodeType,n) 
% wraps n domains from the left back to the right end of the polymer.
% 


% % some test code 
% % startup;
% 
% % Number of nodes
% numNodes =8E3; % 45E3; 
% 
% % Number of domains
% numDomains =10;%  % total small colored domains to study
% maxDomainLength = 500; 
% minDomain = 10; 
% fracYvsK = .3;  % fraction of yellow vs black
% fracSticky =.25; % fraction of blue nodes thar are sticky
% numBlue = 2; 
% 
% [nodeType,domainCoords,domainColors,domainLengths,showDomains] = ...
%     GetRandomColorDomains(numDomains,numNodes,'numBlue',numBlue,...
%     'minDomain',minDomain,'minGap',30,'fracYvsK',fracYvsK,'fracSticky',fracSticky,'lengthenDomains',0,'maxDomainLength',maxDomainLength);
% nodeColor = nodeType;
% nodeColor(nodeColor==4) = 3; 
% 
% showLinear = 1;
% clrMap1 = .5*ones(numNodes,3);
% clrs = [.2 .2 .2;
%         .9 .8 0;
%          0  0 1];     
% for c=1:3
%     clrMap1(nodeColor==c,:) = repmat(clrs(c,:),sum(nodeColor==c),1);
% end
% if showLinear
%     figure(1); clf; 
%     data = repmat(1:numNodes,3,1)';
%     data = (data+.05*rand(size(data)))/400;
%     data(:,1) = sin(data(:,1));
%     data(:,3) = sin(data(:,3));
%     PlotTube(data,'r',1,'subDivisions',3,'interpPts',1); colormap(clrMap1); set(gcf,'color','w');
%     view(-90,-22); axis off;
% end

%% Main Function

% default paremters
plotNewDomains = false;
numNodes = length(nodeType); 

% dealing with multicolor blue
nodeColor = nodeType;
nodeColor(nodeType==4) = 3; % these are both types of blue.  

% figure(1); clf; plot(diff(nodeColor));
domainBreaks = find(diff(nodeColor)>0)-20;
numBreaks = length(domainBreaks);
numDomains = length(domainCoords); 

if numBreaks ~= numDomains
    warning('something super weird happened');
    numBreaks = numDomains;
end

if ~isempty(varargin)
    domainShift = varargin{1};
else
    domainShift = randi(numBreaks-1);
end
nodeShift = numNodes-domainBreaks(domainShift);
newNodeType = circshift(nodeType,nodeShift,2); 


newDomainCoords = domainCoords;
for d=1:numDomains
    newCoords = domainCoords{d} - (numNodes - nodeShift);
    if min(newCoords) < 0
        newCoords = numNodes  - (domainCoords{1+domainShift-1}(1)-20) + domainCoords{d}+1;
    end
    newDomainCoords{d} = newCoords;
end
        
    
%     test = zeros(1,numNodes);
%     for d=1:numDomains
%         test(newDomainCoords{d}) = domainColors(d);
%     end
%     figure(2); clf; imagesc([nodeIdx;newNodeIdx;test]); colorbar;
    
    

if plotNewDomains
    clrMap = .5*ones(numNodes,3);
    clrs = [.2 .2 .2;
            .9 .8 0;
             0  0 1];      
    for c=1:3
        clrMap(newNodeType==c,:) = repmat(clrs(c,:),sum(newNodeType==c),1);
    end

    if showLinear
        figure(2); clf; 
        data = repmat(1:numNodes,3,1)';
        data = (data+.05*rand(size(data)))/400;
        data(:,1) = sin(data(:,1));
        data(:,3) = sin(data(:,3));
        PlotTube(data,'r',1,'subDivisions',3,'interpPts',1); colormap(clrMap); set(gcf,'color','w');
        view(-90,-22); axis off;
    end
end


