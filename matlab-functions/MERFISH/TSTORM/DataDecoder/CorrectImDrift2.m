function [dxc,dyc,fedPos,driftErrors] = CorrectImDrift(beads,ImLists,varargin)
% [dxc,dyc,FedPos] = CorrectImDrift(beads,'maxdrift',maxDrift,...
%                                   'showplots',showDriftPlots,...
%                                   'fmin',minFramesForFeducials);



%--------------------------------------------------------------------------
%% Default Parameters
%--------------------------------------------------------------------------
startframe = 1; % frame to use to find feducials
maxdrift = 2.5; % max distance a feducial can get from its starting position and still be considered the same molecule
integrateframes = 200; % number of frames to integrate
samplingrate = 1; 
fmin = .5; 
fighandle = [];

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
            case 'samplingrate'
                samplingrate = CheckParameter(parameterValue,'positive','samplingrate');
            case 'integrateframes'
                integrateframes = CheckParameter(parameterValue,'positive','integrateframes');
            case 'fmin'
                fmin = CheckParameter(parameterValue,'positive','fmin');
            case 'fighandle'
                fighandle = CheckParameter(parameterValue,'handle','fighandle');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
%% Main Function
%--------------------------------------------------------------------------
numHybes = length(beads); 
cT = cell(numHybes,1);
rT = cell(numHybes,1);
dxc = cell(numHybes,1);
dyc = cell(numHybes,1);
fedPos = cell(numHybes,1);
driftErrors = zeros(numHybes,1);

for i=1:numHybes  % i = 9  
    figure(fighandle);
        if numHybes < 9
            subfighandle = subplot(ceil(numHybes/2),2,i);
        else 
            subfighandle = subplot(ceil(numHybes/4),4,i);
        end
    try
    [dx,dy,cT{i},rT{i},drift_error] = FeducialDriftCorrection(beads{i},... % Compute drift
    'samplingrate',samplingrate,'integrateframes',integrateframes,...
    'maxdrift',maxdrift,'startframe',startframe,'fmin',fmin,...
     'fighandle',subfighandle);
    catch er
        disp(er.getReport);
        disp(er.message);
       cT{i} = cT{i-1};
       rT{i} = rT{i-1}; 
       drift_error = NaN;
    end
    fedPos{i} = [cT{i}(1,:,1);cT{i}(1,:,2)]';                              % Record starting position of beads in each movie
    driftErrors(i) = drift_error;
    
    dx(isnan(dx)) = [];
    dy(isnan(dy)) = [];
    targetmlist = ImLists{i};                                              % Convert beadframes to movie frames 
    targetFrames = double(max(targetmlist.frame)); 
    beadFrames = linspace(1,targetFrames,length(dx));    
    dxc{i} = interp1(beadFrames,dx,1:targetFrames,'linear' )'; 
    dyc{i} = interp1(beadFrames,dy,1:targetFrames,'linear' )';
end