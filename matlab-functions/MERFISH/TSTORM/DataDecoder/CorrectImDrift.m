function [dxc,dyc,FedPos] = CorrectImDrift(beads,varargin)
% [dxc,dyc,FedPos] = CorrectImDrift(beads,'maxdrift',maxDrift,...
%                                   'showplots',showDriftPlots,...
%                                   'fmin',minFramesForFeducials);



%--------------------------------------------------------------------------
%% Default Parameters
%--------------------------------------------------------------------------
startframe = 1; % frame to use to find feducials
maxdrift = 2.5; % max distance a feducial can get from its starting position and still be considered the same molecule
integrateframes = 200; % number of frames to integrate
fmin = .5; 
showplots = true;

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'startframe'
                startframe = CheckParameter(parameterValue,'positive','startframe');
            case 'maxdrift'
                maxdrift = CheckParameter(parameterValue,'positive','maxdrift');
            case 'integrateframes'
                integrateframes = CheckParameter(parameterValue,'positive','integrateframes');
            case 'fmin'
                fmin = CheckParameter(parameterValue,'positive','fmin');
            case 'showplots'
                showplots = CheckParameter(parameterValue,'boolean','showplots');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
%% Main Function
%--------------------------------------------------------------------------
D = length(beads)
dxc = cell(D,1);
dyc = cell(D,1);
Fc = cell(D,1);
FedPos = cell(D,1);
for d=1:D % d=2
  [dxc{d},dyc{d},Fc{d}] = feducialDriftCorrection([beads{d},'_list.bin'],...
      'maxdrift',maxdrift,'showplots',showplots,...
      'fmin',fmin,'integrateframes',integrateframes,...
      'startframe',startframe);
%   [~,~,Fc{d}] = feducialDriftCorrection([beads{d},'_list.bin'],...
%       'maxdrift',maxDrift,'showplots',false,...
%       'fmin', 11/double(max(beadlists{d}.frame)));
%  % Fc is Tframes x Nbeads x 2-dimensions
  temp = nanmedian(Fc{d}(1:10,:,:),1);
  [~,Nfeducials,~] = size(temp);
  FedPos{d} = [temp(:,:,1)',temp(:,:,2)'];
end
