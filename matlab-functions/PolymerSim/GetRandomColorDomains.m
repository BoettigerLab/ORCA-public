function [nodeType,domainCoords,domainColors,domainLengths,showDomains] = GetRandomColorDomains(numDomains,numNodes,varargin)
% specify domain structure

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'fracYvsK', 'fraction', .02};
defaults(end+1,:) = {'numBlue', 'nonnegative', 6};
defaults(end+1,:) = {'showLinear', 'boolean', true};
defaults(end+1,:) = {'showFraction', 'nonnegative', .2};
defaults(end+1,:) = {'maxShowLength', 'nonnegative', 300};
defaults(end+1,:) = {'minDomain', 'nonnegative', 0};
defaults(end+1,:) = {'minGap', 'nonnegative', 0};
defaults(end+1,:) = {'fracSticky', 'nonnegative', 0};
defaults(end+1,:) = {'maxDomainLength', 'nonnegative', 600};
defaults(end+1,:) = {'lengthenDomains', 'nonnegative', 0};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'requires numDomains,numNodes');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments('', defaults, mfilename);

%% Main Function

numBlue = parameters.numBlue; % 6;  
%%

% if parameters.minGap*numDomains > numNodes/numDomains/2
%     error(' minGap*numDomains must be smaller than numNodes/numDomains/2'); 
% end
% 
% if parameters.minDomain*numDomains > numNodes/numDomains/2
%     error(' minDomain*numDomains must be smaller than numNodes/numDomains/2'); 
% end


restart = true;
while restart == true
    e = 0;
    domainColors = zeros(numDomains,1);
    domainCoords = cell(numDomains,1);
    domainLengths = zeros(numDomains,1);
    nodeType = zeros(1,numNodes); 
    
    for d=1:numDomains
        restart = false;
        newDomainLength = inf;
        while newDomainLength >= parameters.maxDomainLength;
            newDomainLength = exprnd(numNodes/numDomains/2-parameters.minDomain + parameters.lengthenDomains);
        end
        s = round(e + parameters.minGap +  exprnd(numNodes/numDomains/2-parameters.minGap + parameters.lengthenDomains ));
        e = round(s + parameters.minDomain + newDomainLength);
        if e>numNodes
            e=numNodes;
        end
        if s>numNodes
            restart = true;
            warning('failed to partition domains, restarting search');
        end
        domainCoords{d} = s:e;
        clrChoice = 1 + (rand < parameters.fracYvsK);
        domainColors(d) = clrChoice;
        domainLengths(d) = e-s+1;
        nodeType(s:e) = domainColors(d);
    end
end
% evenly distribute the blue domains
blueDomainStep = ceil(numDomains/numBlue);
if numBlue>=1
    for b = 1:blueDomainStep:numDomains
        domainColors(b) = 3;
        if parameters.fracSticky == 1
            nodeType(domainCoords{b}) = 3; 
        else
           nodeFlavor = 4-(rand(1,domainLengths(b)) < parameters.fracSticky); % random distribution of sticky nodes  
           nodeType(domainCoords{b}) = nodeFlavor;  
        end
    end
end

% Build color maps 
k = [.2 .2 .2] ; % 
y = [1 .9 0];
b = [0 0 1];
bkd = [.8 .8 .8];
clrs = [k;y;b;bkd];
clrMap = repmat(bkd,numNodes,1);
for d=1:length(domainCoords);
    clrMap(domainCoords{d},:) = repmat(clrs(domainColors(d),:), length(domainCoords{d}),1);
end

if parameters.showLinear
    figure(1); clf; 
    data = repmat(1:numNodes,3,1)';
    data = (data+.05*rand(size(data)))/400;
    data(:,1) = sin(data(:,1));
    data(:,3) = sin(data(:,3));
    PlotTube(data,'r',1,'subDivisions',3,'interpPts',1); colormap(clrMap); set(gcf,'color','w');
    view(-90,-22); axis off;
end
 
 showDomains = (domainLengths < parameters.maxShowLength & rand(numDomains,1) < parameters.showFraction)   | domainColors == 3;

 % Make sure output makes sense
 
if min(cat(2,domainCoords{:})) <= 0 
   error('min less than 0!');
end
if max(cat(2,domainCoords{:})) > numNodes
   error('max greater than numNodes!')
end
 
      