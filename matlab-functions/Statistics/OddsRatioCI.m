function [odds,odds_CI,odds_stdev,oddsM] = OddsRatioCI(condition,exposure,varargin)

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'cI', 'positive', .95};
defaults(end+1,:) = {'iters', 'positive', 10000};
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabFunctions:invalidArguments', 'data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
pars = ParseVariableArguments(varargin, defaults, mfilename);
% pars = ParseVariableArguments([], defaults, mfilename);

% condition = pab;
% exposure = pad;



a = sum(condition & exposure);
d = sum(~condition & ~exposure);
b = sum(~condition & exposure);
c = sum(condition & ~exposure);
odds = a*d/(b*c);

if nargout > 1 
    lowLim = (1 - pars.cI);
    highLim = pars.cI;
    n = length(condition); 
    randSample = randi(n,pars.iters,n);

    Cond = condition(randSample);
    Exp = exposure(randSample);
    A = sum(Cond & Exp ,2);
    D = sum(~Cond & ~Exp, 2);
    B = sum(~Cond & Exp, 2);
    C = sum(Cond & ~Exp, 2);

    oddsM = A.*D./(B.*C);
    oddsM = sort(oddsM);

    odds_CI =  oddsM(round([lowLim,highLim]*pars.iters))';
    odds_stdev = std(oddsM);
else
    odds_CI = [NaN,NaN];
end

