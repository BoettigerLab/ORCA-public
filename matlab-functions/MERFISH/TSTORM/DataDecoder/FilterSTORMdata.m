function [xf,yf] = FilterSTORMdata(ImLists,dxc,dyc,tform,xshift,yshift,varargin)

%--------------------------------------------------------------------------
%% Default parameters
%--------------------------------------------------------------------------
showPlots = true;

xmin = 0;  % Try analyze smaller region
xmax = 256;
ymin = 0;
ymax = 256;

minPhotons =  2E2; % [35E3,35E3,35E3]  ;  % min # photons per spot
maxD = .2; % max distance for density filter (in pixels)
minDensity = 5; % min number of localizations within maxD
startingFrame = 10; 

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 6
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'startframe'
                startingFrame = CheckParameter(parameterValue,'positive','startframe');
            case 'maxDistance'
                maxD = CheckParameter(parameterValue,'positive','maxDistance');
            case 'minDensity'
                minDensity = CheckParameter(parameterValue,'positive','minDensity');
            case 'minPhotons'
                minPhotons = CheckParameter(parameterValue,'positive','minPhotons');
            case 'showPlots'
                showPlots = CheckParameter(parameterValue,'boolean','showPlots');
            case 'xmin'
                xmin = CheckParameter(parameterValue,'nonNegative','xmin');
            case 'xmax'
                xmax = CheckParameter(parameterValue,'nonNegative','xmax');
            case 'ymin'
                ymin = CheckParameter(parameterValue,'nonNegative','ymin');
            case 'ymax'
                ymax = CheckParameter(parameterValue,'nonNegative','ymax');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
%% Main Function
%-------------------------------------------------------------------------
 
D = length(ImLists); 
xf = cell(1,D);
yf = cell(1,D); 
xr = cell(1,D);
yr = cell(1,D); 
xa = cell(1,D);
ya = cell(1,D); 
cMap = jet(D);
amin = minPhotons*ones(D,1);  

for d=1:D % d=2
    % Correct the drift within a movie
    maxframe = length(dxc{d}); % number of frames in bead movie;
    keepframes = ImLists{d}.frame < maxframe; 
    xo = ImLists{d}.x(keepframes) - dxc{d}(ImLists{d}.frame(keepframes)) + xshift(d);
    yo = ImLists{d}.y(keepframes) - dyc{d}(ImLists{d}.frame(keepframes)) + yshift(d); 
    
   % Align different stains
    [x,y] = tforminv(tform{d},double(xo),double(yo));
     
     % Filter on photons and frame
    filt1 = ImLists{d}.frame(keepframes)  >= startingFrame;     
    filt2 = ImLists{d}.a(keepframes)  > amin(d);
    % Filter on local density:
    [~,Dist] = knnsearch([x,y],[x,y],'K',minDensity);
    filt3 = max(Dist,[],2) < maxD;
    filt4 = x > xmin & x < xmax & y > ymin & y < ymax;
    
    % localizations without density filter
    xa{d} = x(filt1 & filt2 & filt4);
    ya{d} = y(filt1 & filt2 & filt4); 
    
    % localizations including density filter
    xf{d} = x(filt1 & filt2 & filt3 & filt4);
    yf{d} = y(filt1 & filt2 & filt3 & filt4);
    xr{d} = x( ~(filt1 & filt2 & filt3 & filt4));
    yr{d} = y( ~(filt1 & filt2 & filt3 & filt4));
    
    if showPlots
        colordef white; set(gcf,'color','w');
        plot(xf{d}, yf{d},'color',cMap(d,:),'Marker','.',...
            'LineStyle','none','MarkerSize',5); hold on; % plots filtered dots 
    end  
end
if showPlots
    legend(num2str([1:D]'))
    for d=1:D
            plot(xr{d}, yr{d},'color',cMap(d,:),'Marker','.',...
        'LineStyle','none','MarkerSize',1); hold on; % plots rejected dots 
    end
end