function BoxPlot2D(x,data,varargin)
% BoxPlot2D(x,d) 
% generates a box plot for the data in each cell entry of d{i}, plotted at
% position x(i).
%
%--------------------------------------------------------------------------
% Optional Inputs
% 'width' / double / .03
% 'datanames' / cell / {}
% 'clrmap' / string or colormap
% 'showdots' / boolean / false
% 'MarkerSize' / double / 30
% 'dotsize' / double / 5

%% Main Function

%% Default parameters

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'width', 'boolean', .9};
defaults(end+1,:) = {'colormap', 'colormap', []};
defaults(end+1,:) = {'color', 'colormap', [0 0 1]};
defaults(end+1,:) = {'datanames', 'cell', {}};
defaults(end+1,:) = {'showdots', 'boolean', false};
defaults(end+1,:) = {'dotsize', 'nonnegative', 5};
defaults(end+1,:) = {'MarkerSize', 'nonnegative', 30};
defaults(end+1,:) = {'wisker', 'nonnegative', 0};
defaults(end+1,:) = {'wiskerSD', 'nonnegative', 0};
defaults(end+1,:) = {'showMean', 'boolean', false};
defaults(end+1,:) = {'showMedian', 'boolean', true};

% allow cell input with no x values
if iscell(x)
    data = x;
    x = 1:length(x);
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


width = parameters.width;
clrmap = parameters.colormap;
datanames = parameters.datanames;
showdots = parameters.showdots;
dotsize = parameters.dotsize;
MarkerSize = parameters.MarkerSize; 
wisker = parameters.wisker;
wiskerSD = parameters.wiskerSD;
showMean = parameters.showMean;
showMedian = parameters.showMedian;

%--------------------------------------------------------------------------
%% Actual Function
%--------------------------------------------------------------------------


xRange = max(x) - min(x);
numDataTypes = length(data);
w = xRange*width;

if isempty(clrmap)
    clrmap = repmat(parameters.color,numDataTypes,1);
end
if ischar(clrmap)
    try
    clrmap = eval([clrmap,'(numDataTypes)']);
    catch 
       disp([clrmap,' is not a valid colormap name']);  
    end
end

quartiles = zeros(numDataTypes,2);
medians = zeros(numDataTypes,1); 
means = zeros(numDataTypes,1); 
for i=1:numDataTypes
    quartiles(i,:) = quantile(data{i},[.25,.75]);
    quarts = quantile(data{i},[.25,.75]);
    medians(i) = nanmedian(data{i});
    means(i)= nanmean(data{i});
    h = max(1E-16,quarts(2)-quarts(1));
    boxes = [x(i)-w/2,quarts(1),w,h];  
    
    
    if showdots
        numPts = length(data{i}); 
        plot( x(i) + w*.9*(.5-rand(numPts,1)),data{i},'.',...
            'color',clrmap(i,:),'MarkerSize',dotsize);
    end
    
    hold on;
    rectangle('Position',boxes,'EdgeColor',clrmap(i,:)); 
    if showMedian
        plot(x(i),medians(i),'.','color',clrmap(i,:),'MarkerSize',MarkerSize);
    end
    if showMean
        plot(x(i),means(i),'+','color',clrmap(i,:),'MarkerSize',MarkerSize);
    end
        
    if ~isempty(datanames)
        text(x(i)+w,medians(i),datanames{i});
    end

    if wisker > 0;  
       wisk = quantile(data{i},[1-wisker,wisker]);
       plot([x(i),x(i)],[wisk(1) quarts(1)],'k-');
       plot([x(i)-w/2,x(i)+w/2],[wisk(1) wisk(1)],'k-');
       plot([x(i),x(i)],[quarts(2) wisk(2)],'k-');
       plot([x(i)-w/2,x(i)+w/2],[wisk(2) wisk(2)],'k-');
    end
     if wiskerSD > 0;  
       wisk(2) = min([max(data{i}), median(i)+wiskerSD*std(data{i})]) ;
       wisk(2) = data{i}(find(data{i} >= wisk(2),1,'first')); 
       wisk(1) =max([min(data{i}), median(i)-wiskerSD*std(data{i}) ]);
       wisk(1) = data{i}(find(data{i} <= wisk(1),1,'first'));
       if wisk(1) < quarts(1)
       plot([x(i),x(i)],[wisk(1) quarts(1)],'k-');
       plot([x(i)-w/2,x(i)+w/2],[wisk(1) wisk(1)],'k-');
       end
       if wisk(2) > quarts(2); 
        plot([x(i),x(i)],[quarts(2) wisk(2)],'k-');
        plot([x(i)-w/2,x(i)+w/2],[wisk(2) wisk(2)],'k-');
       end
    end
    
end

set(gcf,'color','w');

