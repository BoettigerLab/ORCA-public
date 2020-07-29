function [y,yl,yu,errorPatch,linfit,conf95,se] = LinLogFit(xdata,ydata,varargin)


% -------------------------------------------------------------------------
% Optional Inputs
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'minPlotX', 'nonnegative', []};
defaults(end+1,:) = {'maxPlotX', 'nonnegative', []};


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a nx3 data array is required');
end

varinput = [];

if nargin > 3
    varinput = varargin(2:end);
elseif nargin > 2
    x = varargin{1};
else
    x = xdata; % linspace(min(xdata),max(xdata)); 
end



% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varinput, defaults, mfilename);


% fittype('c*x+b','coeff',{'c','b'},'ind','x');
% lower_c = conf95(1,1);
% upper_c = conf95(2,1);
% lower_b = conf95(1,2);
% upper_b = conf95(2,2);


%% Main Function

ftype = fittype('c*x+b','coeff',{'c','b'},'ind','x');
linfit = fit(log10(xdata),log10(ydata),ftype,'StartPoint',[1,0]);
y = x.^linfit.c*10^linfit.b;
try
    conf95 = confint(linfit);
catch
    conf95 = NaN(4);
end
lower_c = conf95(1,1);
upper_c = conf95(2,1);
lower_b = conf95(1,2);
upper_b = conf95(2,2);
% yl = x.^lower_c*10^upper_b;% linfit.b  ;
% yu = x.^upper_c*10^lower_b; % linfit.b  ;

y2u = (log10(x(end))*upper_c+lower_b);
y1u = (log10(x(1))*lower_c+upper_b);
y2l = (log10(x(end))*lower_c+upper_b);
y1l = (log10(x(1))*upper_c+lower_b);
x2 =  log10(x(end));
x1 =  log10(x(1));

slope_upper = (  y2u - y1u )/ ( x2 - x1  );
slope_lower = (  y2l - y1l )/ ( x2 - x1  );
b_upper = y1u - x1*slope_upper;
b_lower = y1l - x1*slope_lower;

yu =  x.^slope_upper*10^b_upper;
yl =  x.^slope_lower*10^b_lower;

[ymin,ymax] = LineConfidenceBounds(log10(xdata),log10(ydata),log10(x));

yl = 10.^ymin;
yu = 10.^ymax;

errorPatch = [[x; x],[ymin; ymax]];

stats = regstats(log10(ydata),log10(xdata),'linear','tstat');  
se = stats.tstat.se;

