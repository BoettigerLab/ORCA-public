function [fitHandle,legendEntry,allfit] = PlotLogLogFit(xs,volume,varargin)
% defaults(end+1,:) = {'color', 'colormap', [0,0,0]};
% defaults(end+1,:) = {'dataName', 'string', ''};
% defaults(end+1,:) = {'showLegend', 'boolean', true};
% defaults(end+1,:) = {'plotPoints','boolean',true};
% defaults(end+1,:) = {'linewidth','positive',2};
% defaults(end+1,:) = {'lineStyle','string','-'};
% defaults(end+1,:) = {'plotFit','boolean',true};
% defaults(end+1,:) = {'markerSize','nonnegative',30};
% defaults(end+1,:) = {'markerStyle','string','.'};
% defaults(end+1,:) = {'markerLineWidth','nonnegative',1};
% defaults(end+1,:) = {'markerFaceColor','freeType','none'};
% defaults(end+1,:) = {'legConf','boolean',true};

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'color', 'colormap', [0,0,0]};
defaults(end+1,:) = {'dataName', 'string', ''};
defaults(end+1,:) = {'showLegend', 'boolean', true};
defaults(end+1,:) = {'plotPoints','boolean',true};
defaults(end+1,:) = {'linewidth','positive',2};
defaults(end+1,:) = {'lineStyle','string','-'};
defaults(end+1,:) = {'lineColor','colormap',[]};
defaults(end+1,:) = {'plotFit','boolean',true};
defaults(end+1,:) = {'markerSize','nonnegative',30};
defaults(end+1,:) = {'markerStyle','string','.'};
defaults(end+1,:) = {'markerLineWidth','nonnegative',1};
defaults(end+1,:) = {'markerFaceColor','freeType','none'};
defaults(end+1,:) = {'legConf','boolean',false};
defaults(end+1,:) = {'legSE','boolean',true};
defaults(end+1,:) = {'sigFigs','positive',2};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a nx3 data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);

%% 
% parameters
fade = 1; 
color = parameters.color; 
dataName = parameters.dataName;
sigfigs = parameters.sigFigs; 

if isempty(parameters.lineColor)
    lineColor = color;
else
    lineColor = parameters.lineColor;
end

% total domain median and errors
% idx = colorData.idx == StringFind(colorData.subGroupNames,'Blue');
numTotData = length(volume);
vM  = zeros(numTotData,1);
vCI  = zeros(numTotData,2);
if iscell(volume)
for i=1:numTotData
    if length(volume{i}) > 1
        [vM(i),vCI(i,:)] = MedianWithCI(volume{i});
    else
        vM(i) = volume{i};
    end
end
else
    vM = volume;
end

if parameters.plotPoints 
    if iscell(volume)
        for i=1:numTotData
            loglog(xs(i),vM(i),parameters.markerStyle,'color',color,...
                'MarkerSize',parameters.markerSize,...
                'MarkerFaceColor',parameters.markerFaceColor,...
                'LineWidth',parameters.markerLineWidth); hold on;
            ploterr(xs(i),vM(i),[],{vCI(i,1),vCI(i,2)},'','logxy','color',color);
        end
    else
    loglog(xs,vM,parameters.markerStyle,'color',color,...
                'MarkerSize',parameters.markerSize,...
                'MarkerFaceColor',parameters.markerFaceColor,...
                'LineWidth',parameters.markerLineWidth); hold on;
    end
    axis tight;
else
        loglog(10,10,'w.'); hold on;
end

skipPts = isnan(vM) | vM==0 | isnan(xs) | xs==0;

if parameters.plotFit
    x =  linspace(min(xs),max(xs))';
    [y,yl,yu,yp,allfit,conf95,se] = LinLogFit(xs(~skipPts),vM(~skipPts),x); % figure(3); clf;
   
    fitHandle = plot(x,y,'color',lineColor,'linewidth',parameters.linewidth,'linestyle',parameters.lineStyle); hold on;
    if parameters.legConf
        legendEntry = [dataName,'{\it c} = ',num2str( allfit.c,sigfigs),' (',num2str(conf95(1,1),sigfigs),',',num2str(conf95(2,1),sigfigs),')'];
    elseif parameters.legSE
        legendEntry = [dataName,'{\it c} =',num2str( allfit.c,sigfigs),' \pm ',num2str(se(2),sigfigs)];
    else
        legendEntry = [dataName,'{\it c} =',num2str( allfit.c,sigfigs)];
    end
else
    parameters.showLegend = false;
end
if parameters.showLegend
    legend(legendEntry,'Location','Best');
end
set(gcf,'color','w');




% for i=1:numTotData;  
%     loglog(xs(i),vM(i),'o','color',color,'MarkerSize',8,'MarkerFaceColor',color); hold on;
%     ploterr(xs(i),vM(i),[],{vCI(i,1),vCI(i,2)},'','logxy','color',color);
% end


% patch( yp(:,1),yp(:,2),ones(size(yp,1),1),'FaceColor',(fade+colorData.colorClass(c,:))/(1+fade),'EdgeColor','none'); hold on;


% 
% 
% for i=1:numTotData;  
%     loglog(locuslengths(i),vM(i),'o','color',colorData.colors(i,:),'MarkerSize',8,'MarkerFaceColor',colorData.colors(i,:)); hold on;
%     ploterr(locuslengths(i),vM(i),[],{vCI(i,1),vCI(i,2)},'','logxy','color',colorData.colors(i,:));
% end