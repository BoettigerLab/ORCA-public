function [xf,yf,flists] = FilterConvData(imLists,tform,varargin)

%--------------------------------------------------------------------------
%% Default parameters
%--------------------------------------------------------------------------
showPlots = true;

xmin = 0;  % Try analyze smaller region
xmax = 256;
ymin = 0;
ymax = 256;

minPhotons =  2E2; % [35E3,35E3,35E3]  ;  % min # photons per spot
startingFrame = 1; 
endFrame = 1; 

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 2
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
            case 'endframe'
                endFrame = CheckParameter(parameterValue,'positive','endFrame');
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
 

numHybes = length(imLists); 
xf = cell(1,numHybes);
yf = cell(1,numHybes); 
flists = cell(1,numHybes);
cMap = jet(numHybes);
amin = minPhotons*ones(numHybes,1);  

for h=1:numHybes % d=2
    % Correct the drift within a movie
    xo = imLists{h}.x;
    yo = imLists{h}.y; 
    
   % Align different stains
   if ~isempty(tform)
    [x,y] = tforminv(tform{h},double(xo),double(yo));
   else
       x=xo;
       y=yo;
   end 
    
     % Filter on photons and frame
    filt1 = imLists{h}.frame  >= startingFrame & imLists{h}.frame <= endFrame;     
    filt2 = imLists{h}.a  > amin(h);
    filt4 = x > xmin & x < xmax & y > ymin & y < ymax;
    
    % localizations without density filter
    xf{h} = x(filt1 & filt2 & filt4);
    yf{h} = y(filt1 & filt2 & filt4); 
    
    flists{h} = IndexStructure(imLists{h},(filt1 & filt2 & filt4));
    flists{h}.xc = xf{h}; % need to apply warp filter
    flists{h}.yc = yf{h}; 
        
    if showPlots
        colordef white; set(gcf,'color','w');
        plot(xf{h}, yf{h},'color',cMap(h,:),'Marker','.',...
            'LineStyle','none','MarkerSize',5); hold on; % plots filtered dots 
    end  
end
if showPlots
    legend(num2str([1:numHybes]'));
end