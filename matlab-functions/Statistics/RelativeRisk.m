function [rr,rr_CI,rr_stdev,RR] = RelativeRisk(condition,exposure,varargin)
%  [rr,rr_CI] = RelativeRisk(condition,exposure)
%  [rr,rr_CI] = RelativeRisk(condition,exposure,'cI',.95)
%  
% Matrix
%                          Exposure:
%                          Yes   No
%                   Yes   [a     b
%       Condition: 
%                   No     c    d]
%
%
%                  a/(a+c)
%            RR =  -------
%                 (b/(b+d)
%
%
%       std(log(rr))/sqrt(N) * 
%




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

if nargout < 2
    pars.iters = 1;  % no need to do iterations if not returning a CI
end

% condition = pab;
% exposure = pad;

lowLim = (1 - pars.cI);
highLim = pars.cI;
n = length(condition); 

rr = sum(condition & exposure)./sum(exposure) / (sum(condition & ~exposure)./sum(~exposure));

if pars.iters > 1
    randSample = randi(n,pars.iters,n);
    C = condition(randSample);
    E = exposure(randSample);
    RR = sum(C & E,2)./sum(E,2) ./ (sum(C & ~E,2)./sum(~E,2));
    RR = sort(RR); 
    rr_CI =  RR(round([lowLim,highLim]*pars.iters))';
    rr_stdev = std(RR);
else
    rr_CI = NaN;
end
