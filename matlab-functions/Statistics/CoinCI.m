function [obs,exp,obsCI,expCI,pValue] = CoinCI(data,varargin)


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

lowLim = (1 - pars.cI);
highLim = pars.cI;
[nObs,nTrials] = size(data);
randSample = randi(nObs,pars.iters,nObs);

% %% 2d
% p1 = data(:,1); p2 = data(:,2); 
% n = length(p1);
% obs = sort(sum(p1(randSample) & p2(randSample)  ,2)/n);
% exp = sort(sum(p1(randSample),2)/n .* sum(p2(randSample),2)/n);
% obsCI = obs(round([lowLim,highLim]*pars.iters))';
% expCI = exp(round([lowLim,highLim]*pars.iters))';
% obs = sum(p1 & p2)/n;
% exp = sum(p1)/n*sum(p2)/n;

%% multi-dimensional version
% compute CI

obs = ones(pars.iters,nObs);
exp = ones(pars.iters,1);
for i=1:nTrials
    p_i = data(:,i);
    obs = obs & p_i(randSample);
    exp = exp .* sum(p_i(randSample),2)/nObs;
end
obs = sum(obs,2)/nObs;
obs = sort(obs);
exp = sort(exp);
obsCI = obs(round([lowLim,highLim]*pars.iters))';
expCI = exp(round([lowLim,highLim]*pars.iters))';

% now compute actual values 
obs = ones(nObs,1);
exp = 1;
for i=1:nTrials
    p_i = data(:,i);
    obs = obs & p_i;
    exp = exp .* sum(p_i)/nObs;
end    
cnt = sum(obs);
obs = cnt/nObs;
exp;

% p-values
pValue = myBinomTest(cnt,nObs,exp,'one');
    