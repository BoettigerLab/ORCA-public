function pearsonCorr = PlotCorr(x,y,varargin)
% pearsonCorr = PlotCorr(x,y)
% pearsonCorr = PlotCorr(x,y)
% creates log-log plot of correlation betweent x and y
% stats = PlotCorr(x,y) returns stats.log10rho, stast.log10pvalue
%  stats.rho and stats.pvalue for the correlation in addition to the plot
%
%  non-zero points are removed from log-log correlation and correlation
%  plot, but not from the linear correlations computed.  

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'MarkerSize', 'positive', 10};
defaults(end+1,:) = {'FontSize', 'positive', 6};
defaults(end+1,:) = {'colorMap', 'colormap', 'jet'};
defaults(end+1,:) = {'color', 'colormap', [0,0,0]};
defaults(end+1,:) = {'showNames', 'boolean', false};
defaults(end+1,:) = {'names', 'array', {}};
defaults(end+1,:) = {'nameBuffer', 'positive', .1};
defaults(end+1,:) = {'remove1s', 'boolean', false};
defaults(end+1,:) = {'log', 'boolean', true};

if length(varargin) == 1
    pointNames = varargin{1};
    varin = {};
else
    varin = varargin;
    pointNames = {};
end

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires x,y');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varin, defaults, mfilename);

if ~isempty(pointNames)
   parameters.showNames = ~parameters.showNames;  % change default behavior
end

% required column vectors;
x = Column(x);
y = Column(y); 

if parameters.log
    non0 = x>0 & y>0;
else
   non0 = true(length(x),1); 
end
nonNaN = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y) );
x = x(non0 & nonNaN);
y = y(non0 & nonNaN);

if parameters.remove1s
    is1 = x==1 | y==1;
    x(is1) = [];
    y(is1) = [];
end


try
if isempty(x)
   warning('no nonzero data to plot'); 
   pearsonCorr.log10rho = [];
   pearsonCorr.log10pvalue = [];
   pearsonCorr.rho = [];
   pearsonCorr.pvalue = []; 
else
    if nargin>2 && parameters.showNames
        if parameters.log
            loglog(x,y,'.','color',parameters.color,'MarkerSize',parameters.MarkerSize); hold on;
        else
            plot(x,y,'.','color',parameters.color,'MarkerSize',parameters.MarkerSize); hold on;
        end
        text(x+parameters.nameBuffer*x,y,pointNames,'FontSize',parameters.FontSize); 
    else
        if parameters.log
            loglog(x,y,'.','color',parameters.color,'MarkerSize',parameters.MarkerSize); hold on;
        else
            plot(x,y,'.','color',parameters.color,'markerSize',parameters.MarkerSize); hold on;
        end
    end
    
    [c1,p1] = corr(x,y);
    c0 = 0; p0 =0;
    if parameters.log
        [c0,p0] = corr(log10(x),log10(y));
        title(['\rho_l_o_g_1_0 = ',num2str(c0,2),' (p=',num2str(p0,2),')',...
            '  \rho = ',num2str(c1,2),' (p=',num2str(p1,2),')']);
    else
      title(['  \rho = ',num2str(c1,2),' (p=',num2str(p1,2),')']);
    end
    pearsonCorr.log10rho = c0;
    pearsonCorr.log10pvalue = p0;
    pearsonCorr.rho = c1;
    pearsonCorr.pvalue = p1; 
end
set(gcf,'color','w');


catch er
    disp(er.getReport);
    error('debug here');
end