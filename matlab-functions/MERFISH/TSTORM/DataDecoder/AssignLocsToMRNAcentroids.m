function [Mcents,datPos] = AssignLocsToMRNAcentroids(in1,in2,mRNAcents,varargin)
% [Mcents,datPos] = AssignLocsToMRNAcentroids(xf,yf,mRNAcents,varargin)

%% DefaultParameters
minLocsPerStain = 3; % min number of localizations within MaxDtoCentroid to be called 'detected in stain d'
maxDtoCentroid = .5; % max distance you can be from the centroid of the mRNA
method = 'dots';

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 3
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'minLocsPerStain'
                minLocsPerStain = CheckParameter(parameterValue,'positive','minLocsPerStain');
            case 'maxDtoCentroid'
                maxDtoCentroid = CheckParameter(parameterValue,'positive','maxDtoCentroid');
            case 'method'
                method = CheckParameter(parameterValue,'string','method');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%% Main Function
numHybes = length(in1); 

if strcmp(method,'clusters')
    stainCents = in1;
    stainCounts = in2;
else
    xf = in1;
    yf = in2; 
end

% minLocsPerStain = [5 ,3];
if length(minLocsPerStain) < numHybes
    minLocsPerStain = [minLocsPerStain, repmat(minLocsPerStain(end),1, numHybes - length(minLocsPerStain)) ];
end

% Assign localizations to mRNA centroids
datPos = cell(numHybes,1); 
Mcents = zeros(length(mRNAcents),numHybes);
for d=1:numHybes
    if strcmp(method,'clusters')
        positionsInStain = stainCents{d}(stainCounts{d} > minLocsPerStain(d),:); 
        [idxstain,Dist2stain] = knnsearch(positionsInStain,mRNAcents,'K',1);  
        Mcents(:,d) = max(Dist2stain,[],2)<maxDtoCentroid;
        datPos{d} = positionsInStain(idxstain,:);
    end
    
%      figure(7); clf; 
%       plot(mRNAcents(:,1),mRNAcents(:,2),'k+'); hold on;
%      plot(positionsInStain(:,1),positionsInStain(:,2),'mo');
%      
%      figure(5); clf; 
%      plot(mRNAcents(:,1),mRNAcents(:,2),'k+');
%     hold on; plot(datPos{d}(:,1),datPos{d}(:,2),'mo');

    if strcmp(method,'dots')
    [idx,Dist] = knnsearch([xf{d},yf{d}],mRNAcents,'K',minLocsPerStain(d));
    xx = [];
    yy = [];
    for m=1:minLocsPerStain(d)
        xx = [xx, xf{d}(idx(:,m))];
        yy = [yy, yf{d}(idx(:,m))];
    end
    datPos{d} = [mean(xx,2),mean(yy,2)];
    Mcents(:,d) = max(Dist,[],2)<maxDtoCentroid;
    end
end





