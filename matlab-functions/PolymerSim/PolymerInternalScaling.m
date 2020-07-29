function [locLength, rGint] = PolymerInternalScaling(data,varargin)
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'minLength', 'positive', 1};

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


%% Internal scaling
 
N = length(data); 
% sects = fliplr(1:N)';
rGint = zeros(N,1); 
locLength = zeros(N,1); 
for s = 1:N
%    sect = sects(s); 
%     nmin = max(floor(N/sect),1);
%     nmax =  min(N,ceil(2*N/sect));
    nmin = parameters.minLength;
    nmax =  s;
    locLength(s) = nmax-nmin; 
    rGint(s) = RadiusOfGyration(data(nmin:nmax,:));
end

%% 



