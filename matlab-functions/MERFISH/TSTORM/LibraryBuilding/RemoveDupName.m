function ProbeData = RemoveDupName(ProbeData,varargin)
%--------------------------------------------------------------------------
% Remove probedata entry that has non-unique name.
%
%% Parse input
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
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end
    
%% Main Function
%-------------------------------------------------------------------------
[~,ia] = unique([ProbeData.CommonName]);
DupIdx = setdiff(1:length([ProbeData.CommonName]),ia);
NonUniqueGenesIdx =[];
for n = 1:length(DupIdx)
    NonUniqueGenesIdx = [NonUniqueGenesIdx find(strcmp([ProbeData.CommonName], ProbeData(DupIdx(n)).CommonName)==1)];
%     ProbeData(strcmp([ProbeData.CommonName], ProbeData(DupIdx(n)).CommonName)).CommonName
end
ProbeData(NonUniqueGenesIdx)=[];
end