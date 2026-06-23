function [odds,odds_CI,contengencyTable,odds_stdev,oddsM] = OddsRatioCI(condition,exposure,varargin)

% change 05/21/2026 by Derek Le:
% (1) output contengencyTable to see counts in the contingency table
% (2) added Haldane-Anscombe correction to account for possibility of 0's 
% in the 2x2 contingency tables 
% (3) fixed lowLim and highLim for cI calculation


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'haldane', 'boolean', true};
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


a = sum(condition & exposure);
d = sum(~condition & ~exposure);
b = sum(~condition & exposure);
c = sum(condition & ~exposure);

contengencyTable = [a,b; c,d]; % aka confusion matrix 

% Haldane-Anscombe correction: add 0.5 to every cell if any cell is zero
if pars.haldane && (a==0 || b==0 || c==0 || d==0)
    a = a + 0.5;
    b = b + 0.5;
    c = c + 0.5;
    d = d + 0.5;
end

odds = a*d/(b*c);

if nargout > 1 
	lowLim = (1 - pars.cI)/2;
	highLim = (1 - pars.cI)/2 + pars.cI;
    n = length(condition); 
    randSample = randi(n,pars.iters,n);

    Cond = condition(randSample);
    Exp = exposure(randSample);
    A = sum(Cond & Exp ,2);
    D = sum(~Cond & ~Exp, 2);
    B = sum(~Cond & Exp, 2);
    C = sum(Cond & ~Exp, 2);

    % Haldane-Anscombe correction also for bootstrap replicates
    if pars.haldane
        zeroCell = (A==0) | (B==0) | (C==0) | (D==0);
        A(zeroCell) = A(zeroCell) + 0.5;
        B(zeroCell) = B(zeroCell) + 0.5;
        C(zeroCell) = C(zeroCell) + 0.5;
        D(zeroCell) = D(zeroCell) + 0.5;
    end

    oddsM = A.*D./(B.*C);
    oddsM = sort(oddsM);

    odds_CI =  oddsM(round([lowLim,highLim]*pars.iters))';
    odds_stdev = std(oddsM);
else
    odds_CI = [NaN,NaN];
end

