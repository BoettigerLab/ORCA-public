
function fighand = normhistall(data,varargin)
%--------------------------------------------------------------------------
% normhistall(data)       -- histograms of the data in each cell of data
% normhistall(data,bins)  -- histograms of the data in each cell of data
% normhistall(data,x)     -- break data evenly into x bins
% normhistall(...,'groupnames',cell) -- cell of names for figure legend
% normhistall(...,'xlabel',string)   -- string for xlabel
% normhistall(...,'alpha',alphavalue) -- zero for alpha off
% normhistall(...,'colormap',colormap-name) -- colormap (e.g. 'hsv')
%--------------------------------------------------------------------------
% Alistair Boettiger
% CC BY 04/25/13
%--------------------------------------------------------------------------
% 

%--------------------------------------------------------------------------
%% Default Parameters
%--------------------------------------------------------------------------
groupnames = repmat({''},length(data),1); 
clrmap = 'hsv';
fighand = [];
xlab = '';
alphavalue = .15;
%--------------------------------------------------------------------------
%% Parse Variable Input Parameters
%--------------------------------------------------------------------------

% normhistall(data) 
% normhistall(data,'option',value);
if nargin == 1  || (nargin > 1 && ischar(varargin{1})) 
    [h,w] = size(data{1});
    if h>w
        catdim = 1;
    else
        catdim = 2;
    end
    alldat = cat(catdim,data{:});
    x = linspace(min(alldat),max(alldat),10);
    
% normhistall(data,bins)   
% normhistall(data,bins,'option',value)
elseif nargin > 1 && length(varargin{1}) == 1  
    [h,w] = size(data{1});
    if h>w
        catdim = 1;
    else
        catdim = 2;
    end
    alldat = cat(catdim,data{:});
    x = linspace(min(alldat),max(alldat),varargin{1}); 
    varargin(1) = [];
    
% normhistall(data,x)
% normhistall(data,x,'option',value); 
elseif nargin > 1 && ~ischar(varargin{1}) == 1  
    x = varargin{1}; 
    varargin(1) = [];
end

if nargin > 2 
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'groupnames'
                groupnames = CheckParameter(parameterValue, 'cell', 'groupnames');
            case 'xlabel'
                xlab = CheckParameter(parameterValue, 'string', 'xlabel');
            case 'colormap'
                clrmap = CheckParameter(parameterValue, 'colormap', 'colormap');
            case 'alpha'
                alphavalue = CheckParameter(parameterValue, 'nonnegative', 'alpha');
            otherwise
                error(['The parameter ''', parameterName,...
                    ''' is not recognized by the function, ''',...
                    mfilename '''.' '  See help ' mfilename]);
        end
    end
end

% more automatic parsing of input
D = length(data);
if ischar(clrmap)
    cmap =eval([clrmap,'(D+1)']);
else
    cmap = clrmap;
end

if alphavalue == 0
    face = 0;
else
    face = 1;
end

%--------------------------------------------------------------------------
%% Main function
%--------------------------------------------------------------------------

if iscell(x);
    for d=1:D
        h1 = hist(data{d},x{d});
        hnorm = sum(h1);
        bar(x{d},h1/hnorm,1);
        hold on;
    end
    try
        x = cat(2,x{:});
    catch er
        disp(er.message)
        error('multiple x-intervals must be passed as row vectors'); 
    end
else 
    for d=1:D
        h1 = hist(data{d},x);
        hnorm = sum(h1);
        bar(x,h1/hnorm,1);
        hold on;
    end
end
if ~isempty(groupnames{1})
 legend_hand = legend(groupnames);
 lh = findobj(legend_hand,'Type','patch');
else 
    lh = [];
end
 h = findobj(gca,'Type','patch');


for d=1:D
    if face ==1 
        set(h(D-d+1),'FaceColor',cmap(d,:),'EdgeColor',cmap(d,:),'LineWidth',2);
        if ~isempty(lh)
            set(lh(D-d+1),'FaceColor',cmap(d,:),'EdgeColor',cmap(d,:));
        end
    else 
    set(h(D-d+1),'FaceColor','none','EdgeColor',cmap(d,:),'LineWidth',2);
        if ~isempty(lh)
            set(lh(D-d+1),'FaceColor','none','EdgeColor',cmap(d,:));
        end
    end
end
xrange = max(x)-min(x);
xlim([min(x)-.1*xrange,max(x)+.1*xrange]);
xlabel(xlab,'FontSize',14);
% set(gca,'FontSize',14);
set(gcf,'color','w');
if alphavalue <= 1 && alphavalue > 0
    alpha(alphavalue);
end
% PresentationPlot('LineWidth',0);

