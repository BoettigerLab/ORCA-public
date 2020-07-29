function [rho,p] = NanCorr(x,y,varargin)
% like corr but removes NaNs 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'remove1s', 'boolean', false};
defaults(end+1,:) = {'log', 'boolean', false};


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires x,y');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
pars= ParseVariableArguments(varargin, defaults, mfilename);


% required column vectors;
x = Column(x);
y = Column(y);
if pars.remove1s
    x(x==1) = nan;
    y(y==1) = nan;
end
if pars.log
    x = log10(x);
    y = log10(y);
end
nonNaN = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y) );

if isempty(x(nonNaN))
    if pars.verbose
        warning('no nonzero data'); 
    end
    rho = nan;
    p = nan;
else
    [rho,p] = corr(x(nonNaN),y(nonNaN));
end
