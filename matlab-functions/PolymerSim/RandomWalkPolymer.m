function P = RandomWalkPolymer(N,varargin)
% returns a 3D random-walk of N links

showPlots = false;
dim = 3;
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
            case 'showPlots'
                showPlots = CheckParameter(parameterValue,'boolean','showPlots');
            case 'dim'
                dim = CheckParameter(parameterValue,'integer','dim');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end


%% Main Function


M = 2*N+2;
P = zeros(N,3); % not placed
P(1,:) = [M/2,M/2,M/2]; % center start;
n=1;
while n<N
    validlink = false;
    tries = 0; % number of attempts to find valid lattice step
    chainIdx = sub2indFast([M,M,M],P(1:n,1),P(1:n,2),P(1:n,3));
    while ~validlink && tries < 100;
        tries = tries + 1; % count attempts to find valid lattice step
        stepDir = randi([1,dim],1);
        if stepDir == 1
            newlink = round([P(n,1) + 2*(rand-.5),P(n,2), P(n,3)]);
        elseif stepDir == 2
            newlink = round([P(n,1), P(n,2) + 2*(rand-.5),P(n,3)]);
        elseif stepDir == 3
            newlink = round([P(n,1),P(n,2), P(n,3) + 2*(rand-.5)]);
        end
       
        % Use linear indexing to determine if object vertices overlap
        newlinkIdx = sub2indFast([M,M,M],newlink(1,1),newlink(1,2),newlink(1,3));
        validlink = sum(chainIdx==newlinkIdx) == 0;
    end
    if validlink
        P(n+1,:) = newlink;
        n = n+1; 
    else % if tries is exceeded, back up 10 links and start again; 
        disp('chain stuck, backing up...'); 
        backsteps = min(n-3,20);
        P(n-backsteps+1:n,:) = zeros(backsteps,3);
        n = n-backsteps;
    end
        
end

if showPlots
    figure(2); clf;
    plot3(P(:,1),P(:,2),P(:,3)); hold on;
    plot3(P(:,1),P(:,2),P(:,3),'k.');
    xlim([min(P(:,1)),max(P(:,1))]);
    ylim([min(P(:,2)),max(P(:,2))]);
    zlim([min(P(:,3)),max(P(:,3))]);
end