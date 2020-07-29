function [len, rG_w] = PolymerInternalScalingWindow(data,varargin)
% [len, rG_w] = PolymerInternalScalingWindow(data,varargin)
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
% Optional Inputs
% -------------------------------------------------------------------------
% maxSamples
% maxWindow
% minWindow
% windowStep

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'maxSamples', 'positive', 20};
defaults(end+1,:) = {'maxWindow', 'positive', 5E3};
defaults(end+1,:) = {'minWindow', 'positive', 10};
defaults(end+1,:) = {'numWindows', 'positive', 100};
defaults(end+1,:) = {'windowStep', 'positive', 10};
defaults(end+1,:) = {'scale', 'string', 'log'};


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

%% Internal scaling
 
N = length(data); 
maxWindow = min(parameters.maxWindow,N-1);

if strcmp(parameters.scale,'log')
    windows = round( logspace( log10( parameters.minWindow ), log10(maxWindow),parameters.numWindows) );
else
     windows = round(linspace(parameters.minWindow,maxWindow,parameters.numWindows) );
end
rG_w = zeros(length(windows),1);
len =  zeros(length(windows),1);
for k = 1:length(windows)
    window = windows(k);
    rG_window = NaN(N-window,1); 
    ws = 1:N-window;
    wi= randperm(length(ws),min(length(ws),parameters.maxSamples));
    wsTest = ws(wi);
    for w = wsTest;
        rG_window(w) = RadiusOfGyration( data(w:w+window-1,:) );
    end
    rG_w(k) = nanmedian(rG_window);
    len(k) = window;
end

