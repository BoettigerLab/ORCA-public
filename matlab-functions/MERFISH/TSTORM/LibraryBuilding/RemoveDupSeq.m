function ProbeData = RemoveDupSeq(ProbeData,varargin)
%--------------------------------------------------------------------------
% Remove probedata entry that has non-unique sequence. Leave only one
% sequence if repeated sequences are found on the same gene. Remove all duplicates if repeated sequences are found on different gene.      
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
ProbeData([ProbeData.Nprobes]==0)=[];

AllSeqs = vertcat(ProbeData.Sequence);
AllGeneNumber = vertcat(ProbeData.GeneNumber);
[~,ia] = unique(AllSeqs);
DupIdx = setdiff(1:length(AllSeqs),ia);
rptSeqs = {AllSeqs{DupIdx(:)}};
SameGene = zeros(length(rptSeqs),1);
for r=1:length(rptSeqs)
    rptIndex = find(strcmp(AllSeqs,rptSeqs{r}));
    n = AllGeneNumber(rptIndex);
    SameGene(r) = sum(n == n(1))/length(n);
    if SameGene(r) == 1 %Leave only one if repeated sequences are found on the same gene
        rptIdx = find(strcmp(vertcat(ProbeData(AllGeneNumber(DupIdx(r))).Sequence), rptSeqs{r}));
        ProbeData(AllGeneNumber(DupIdx(r))).Sequence(rptIdx(2:end)) = [];
        ProbeData(AllGeneNumber(DupIdx(r))).FivePrimeEnd(rptIdx(2:end)) = [];
%          display(['Same gene ' ProbeData(AllGeneNumber(DupIdx(r))).CommonName AllGeneNumber(DupIdx(r))]);
        ProbeData(AllGeneNumber(DupIdx(r))).Nprobes = length(ProbeData(AllGeneNumber(DupIdx(r))).Sequence);
    else %Remove all duplicates if repeated sequences are found on different gene       
        for i = 1:length(n)
            ProbeData(n(i)).Sequence(strcmp(ProbeData(n(i)).Sequence, rptSeqs{r})) = [];
            ProbeData(n(i)).FivePrimeEnd(strcmp(ProbeData(n(i)).Sequence, rptSeqs{r})) = [];
%              display(['Different gene ' ProbeData(n(i)).CommonName n(i)]);
            ProbeData(n(i)).Nprobes = length(ProbeData(n(i)).Sequence);
        end    
    end
end
end