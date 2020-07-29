function COntSeq = CodonOptimize(varargin)
%--------------------------------------------------------------------------
% COntSeq = CodonOptimize(varargin)
% This file returns an E. coli nt sequence that has been codon optimized in 
% that its original codons have been replaced with synonymous codons with
% probability roughly equivalent to those used
%--------------------------------------------------------------------------
% Necessary Inputs
% ntSeq/char array: a string containing a nt sequence
%
%--------------------------------------------------------------------------
% Outputs
% COntSeq/char array: a string containg a nt sequence with the codons
% randomnly assigned based on E. coli usage
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% July 3, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
flags = {'verbose'};

%--------------------------------------------------------------------------
% E. coli codon usage
%--------------------------------------------------------------------------
% Taken from http://openwetware.org/wiki/Escherichia_coli/Codon_usage
aaCodes = {'G', 'E', 'D', 'V', 'A', 'R', 'K', 'N', 'M', 'I', 'T', 'W', 'C', ...
    '*', 'Y', 'F', 'S', 'Q', 'H', 'L', 'P'};

codons = { {'GGG', 'GGA', 'GGT', 'GGC'}, {'GAG', 'GAA'}, {'GAT', 'GAC'}, ...
    {'GTG', 'GTA', 'GTT', 'GTC'}, {'GCG', 'GCA', 'GCT', 'GCG'}, ...
    {'AGG', 'AGA', 'CGG', 'CGA', 'CGT', 'CGC'}, {'AAG', 'AAA'}, ...
    {'AAT', 'AAC'}, {'ATG'}, {'ATA', 'ATT', 'ATC'}, {'ACG', 'ACA', 'ACT', 'ACC'}, ...
    {'TGG'}, {'TGT', 'TGC'}, {'TAG', 'TAA', 'TGA'}, {'TAT', 'TAC'}, ...
    {'TTT', 'TTC'}, {'AGT', 'AGC', 'TCG', 'TCA', 'TCT', 'TCC'}, ...
    {'CAG', 'CAA'}, {'CAT', 'CAC'}, {'TTG', 'TTA', 'CTG', 'CTA', 'CTT', 'CTC'}, ...
    {'CCG', 'CCA', 'CCT', 'CCC'}};

freq = { [.15 .11 .34 .40], [.31 .69], [.63 .37], ...
    [.37 .15 .26 .22], [.36 .21 .16 .27], ...
    [.02 0.04 .1 0.06 .38 .40], [.23 .77], ...
    [.45 .55], [1], [.07 0.51 .42], [.27 .13 .16 .44], ...
    [1], [.45 .55], [0.07 0.64 0.29], [.57 .43], ...
    [.57 .43], [.15 .28 .14 .13 .15 .15], ...
    [.65 .35], [.57 .43], [.13 .13 .50 .04 .10 .10], ...
    [.52 .20 .16 .12]};


%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = true;

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin < 1
    error('A nt sequence is required')
end

ntSeq = varargin{1};
varargin = varargin(2:end);

if mod(length(ntSeq), 3)
    ntSeq = ntSeq(1:3*floor(length(ntSeq)/3));
end

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'verbose'
                verbose = CheckParameter(parameterValue,'boolean','verbose');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Generate New nt Sequence
%--------------------------------------------------------------------------
COntSeq = ntSeq;

for i=1:(length(ntSeq)/3)
    aa = nt2aa(ntSeq( (3*(i-1)+1) : (3*i) ) , 'AlternativeStartCodons', false); % Translate
    
    ind = find(strcmp(aa, aaCodes)); % Get frequencies
    aaFreq = freq{ind};
    
    prob = cumsum([0 aaFreq]);
    
    p = rand(1);
    newCodonInd = find(prob(1:(end-1)) <= p & prob(2:end) > p);
    possibleCodons = codons{ind};
    COntSeq((3*(i-1)+1) : (3*i)) = possibleCodons{newCodonInd};
end
    
% for i=1:length(aaCodes)
%     display('------------------------------------------------------------');
%     display(aaCodes{i});
%     specCodons = codons{i};
%     display(nt2aa([specCodons{:}], 'AlternativeStartCodons', false));
%     display(num2str(sum(freq{i})));
% end     
