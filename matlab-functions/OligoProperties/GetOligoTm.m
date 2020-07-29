function oligoTm = GetOligoTm(sequence, varargin)
%OLIGOPROP calculates properties of DNA oligonucleotide.
%
%   OLIGOPROP(SeqNT) returns oligonucleotide properties for a DNA sequence in
%   a structure with the following fields:
%
%       Tm: vector with melting temperature values, in degrees Celsius,
%       calculated by six different methods, listed in the following order:
%            Basic (Marmur et al., 1962)
%            Salt adjusted (Howley et al., 1979)
%            Nearest-neighbor (Breslauer et al., 1986)
%            Nearest-neighbor (SantaLucia Jr. et al., 1996)
%            Nearest-neighbor (SantaLucia Jr., 1998)
%            Nearest-neighbor (Sugimoto et al., 1996)
%       Ambiguous N characters in SeqNT are considered to potentially be
%       any nucleotide. If SeqNT contains ambiguous N characters, Tm is the
%       midpoint value, and its uncertainty is expressed by Tmdelta.
%
%
%   OLIGOPROP(...,'SALT',SALT) specifies a salt concentration in moles per
%   liter for melting temperature calculations. The default is 0.05M.
%
%   OLIGOPROP(...,'TEMP',TEMP) specifies the temperature for nearest
%   neighbor calculations of free energy. The default is 25 degrees Celsius.
%
%   OLIGOPROP(...,'PRIMERCONC',PRIMERCONC) specifies what concentration
%   will be used for melting temperature calculations. The default is 50e-6M.
%
%   OLIGOPROP(...,'HPBASE',HAIRPINBASE) specifies the minimum number of
%   paired bases forming the neck of the hairpin. The default is 4.
%
%   OLIGOPROP(...,'HPLOOP',HAIRPINLOOP) specifies the minimum number of
%   bases forming a hairpin. The default is 2.
%
%   OLIGOPROP(...,'DIMERLENGTH',LENGTH) the minimum number of aligned bases
%   between the sequence and its reverse. The default is 4.
%
%   Examples:
%
%       S1 = oligoprop(randseq(25))
%       S2 = oligoprop('ACGTAGAGGACGTN')
%
%   See also NTDENSITY, PALINDROMES, PRIMERDEMO, RANDSEQ, ISOELECTRIC,
%   MOLWEIGHT.

%   Copyright 2005-2012 The MathWorks, Inc.

%   References:
%   [1] Panjkovich A., Melo F., "Comparison of different melting
%       temperature calculation methods for short DNA sequences",
%       Bioinformatics 21(6):711-722 (2004).
%   [2] Breslauer K.J., Frank R., Blocker H., Marky L.A., "Predicting DNA
%       duplex stability from the base sequence" PNAS 83:3746-3750 (1986).
%   [3] Sugimoto N., Nakano S., Yoneyama M., Honda K., "Improved
%       thermodynamic parameters and helix initiation factor to predict
%       stability of DNA duplexes", Nucleic Acids Research
%       24(22):4501-4505 (1996).
%   [4] SantaLucia J., "A unified view of polymer, dumbbell, and
%       oligonucleotide DNA nearest-neighbor thermodynamics", PNAS
%       95:1460-1465 (1998).
%   [5] SantaLucia, J., Allawi, H.T., Seneviratne P.A., "Improved
%       Nearest-Neighbor Parameters for Predicting DNA Duplex Stability",
%       Biochemistry 35:3555-3562 (1996).
%   [6] Howley P.M., Israel M.A., Law M.F., Martin M.A., "A rapid method
%       for detecting and mapping homology between heterologous DNAs.
%       Evaluation of polyomavirus genomes", Journal of Biological
%       Chemistry, 254(11):4876-4883 (1979).
%   [7] Marmur J., Doty P., "Determination of the base composition of
%       deoxyribonucleic acid from its thermal denaturation temperature",
%       JMB 5:109-118 (1962).
%   [8] Chen S.H., Lin C.Y., Lo C.Z., Cho C.S., Hsiung C.A., "Primer Design
%       Assistant (PDA): a Web-based Primer Design Tool" Nucleic Acids
%       Research 31:3751-3754 (2003).
%   [9] http://www.basic.northwestern.edu/biotools/oligocalc.html for
%       weight calculations

% determine the format of the sequence
if isstruct(sequence)
    sequence = bioinfoprivate.seqfromstruct(sequence);
end
seq = upper(sequence);


if(isnumeric(seq))
    seq = upper(int2nt(seq));
else
    seq = upper(seq);
end

if (~all(seq=='A'|seq=='C'|seq=='G'|seq=='T'|seq=='N'))
    error(message('bioinfo:oligoprop:IncorrectSequenceType'));
elseif(any(seq=='N'))
    nFlag=1;
else
    nFlag=0;
end

% defaults parameters
minHpinLoop = 2;      % minimum number of bases in the hairpin loop
minHpinBases = 4;     % minimum number of bases in the hairpin stem
temp = 25;            % temperature in Celsius
salt = 0.05;          % salt concentration in moles per liter (M)
primerConc= 50e-6;    % concentration of primers in mole per liter (M)
dimerLength = 4;      % minimum number of bases for dimers
weight = [313.21 289.18 329.21 304.2]; % molecular weight A,C,G,T respectively

% check arguments
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:oligoprop:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'SALT','PRIMERCONC','HPBASE','HPLOOP','DIMERLENGTH','TEMP'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strcmpi(upper(pname), okargs)); %#ok
        if isempty(k)
            error(message('bioinfo:oligoprop:UnknownParameterName', pname));
        else
            switch(k)
                case 1  % SALT
                    if (pval > -1 && isnumeric(pval) && isreal(pval))
                        salt =  pval;
                    else
                        error(message('bioinfo:oligoprop:BadSaltConcentrationParam'));
                    end
                case 2  % PRIMERCONC
                    if (pval > -1 && isnumeric(pval) && isreal(pval))
                        primerConc =  pval;
                    else
                        error(message('bioinfo:oligoprop:BadConcentrationParam'));
                    end
                case 3  % HPBASE
                    if (isnumeric(pval) && isreal(pval) && pval > 1 && pval < length(seq))
                        minHpinBases =  pval;
                    else
                        error(message('bioinfo:oligoprop:BadHairpinBaseParam'));
                    end
                case 4  % HPLOOP
                    if (isnumeric(pval) && isreal(pval) && pval > 1 && pval < length(seq))
                        minHpinLoop =  pval;
                    else
                        error(message('bioinfo:oligoprop:BadHairpinLoopParam'));
                    end
                case 5  % DIMERlength
                    if (isnumeric(pval) && isreal(pval) && pval < length(seq) && pval > 1)
                        dimerLength =  pval;
                    else
                        error(message('bioinfo:oligoprop:BadDimerLengthParam'));
                    end
                case 6  % Temp
                    if (isnumeric(pval) && isreal(pval))
                        temp =  pval;
                    else
                        error(message('bioinfo:oligoprop:BadtempParam'));
                    end
            end
        end
    end
end

% compute sequence properties
numSeq = double(nt2int(seq));
baseNum = [sum(numSeq == 1) sum(numSeq == 2) sum(numSeq == 3) sum(numSeq == 4) sum(numSeq == 15)];

if(numel(numSeq)<8)
    warning(message('bioinfo:oligoprop:SeqLengthTooShort'));
end

if (~nFlag) % no ambiguous symbols 'N'

  
    [tm tmdelta NN NNdelta]  = get_tm_NN(numSeq, baseNum, salt, primerConc, nFlag);

else % occurrences of ambiguous N symbols
    % warning('bioinfo:oligoprop:ambiguousCheck', 'Ambiguous symbols N in input sequence.');

      [tm tmdelta NN NNdelta]  = get_tm_NN(numSeq, baseNum, salt, primerConc, nFlag);
end

% build output structure
oligoTm = tm;



% FUNCTIONS


% calculate melting temperature and thermo values (37 degrees C)
function [tm tmdelta NN NNdelta]  = get_tm_NN(numSeq, baseNum, salt, primerConc, nFlag)
% melting temperatures and thermodynamic values are returned as average +/- delta level.
% If no ambiguous symbols are present, the delta level is zero.

selfCompFlag = all(5-numSeq == numSeq(end:-1:1));

if(selfCompFlag) % self complementary sequence
    b = 1; % correction value for nearest neighbor melting temperature calculation
else
    b = 4;
end

tmdelta = zeros(1,6);

if (~nFlag) % no ambiguous symbols
    if (sum(baseNum)<14)
        basic = 2 * (baseNum(1) + baseNum(4)) + 4 * (baseNum(2) + baseNum(3)); % TM BASIC [9]
        saltadj = basic - 16.6 * log10(0.05) + 16.6 * log10(salt); % TM SALT ADJUSTED [9]
    else
        basic = 64.9 + (41 * ((baseNum(3) + baseNum(2) - 16.4) / sum(baseNum))); %TM BASIC [1],[9]
        saltadj = 100.5 + (41 * ((baseNum(3) + baseNum(2))/ sum(baseNum))) - (820/sum(baseNum)) + (16.6 * log10(salt)); %TM SALT ADJUSTED [9]
    end

    [NN, NNdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag);
    tm = (NN(:,1) * 1000 ./ (NN(:,2) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15; %TM NEAREST NEIGHBOR

else % occurrences of 'N'
    if(sum(baseNum)<14)
        basic = 2 * (baseNum(1) + baseNum(4)) + 4 * (baseNum(2) + baseNum(3)) + 3 * baseNum(5); % TM BASIC [9]
        tmdelta(1) = baseNum(5);
        saltadj = basic - 16.6 * log10(0.05) + 16.6 * log10(salt); % TM SALT ADJUSTED [9]
        tmdelta(2) = 3 * baseNum(5);
    else
        basic = 64.9 + (41 * (baseNum(3) + baseNum(2) - 16.4 + baseNum(5)/2)/sum(baseNum)); % avg TM BASIC
        tmdelta(1)= 1/2 * 41 * (baseNum(5)/sum(baseNum));
        saltadj = 100.5 - (820/sum(baseNum)) + (16.6 * log10(salt)) +  41 * ((baseNum(2) + baseNum(3)+ baseNum(5)/2)/sum(baseNum));% avg TM SALT ADJUSTED
        tmdelta(2)= 41 * 1/2 * (baseNum(5) / sum(baseNum));
    end
    [NN, NNdelta] = near_neigh(numSeq, length(numSeq), selfCompFlag, nFlag);
    tm = (((NN(:,1)+ NNdelta(:,1)) * 1000 ./ ((NN(:,2)+ NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 + ...
        ((NN(:,1)- NNdelta(:,1)) * 1000 ./ ((NN(:,2)- NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)* 1/2  ; % NEAREST NEIGHBOR
    tmdelta(3:6)=(((NN(:,1)+ NNdelta(:,1)) * 1000 ./ ((NN(:,2)+ NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15 - ...
        (((NN(:,1)- NNdelta(:,1)) * 1000 ./ ((NN(:,2)- NNdelta(:,2)) + (1.9872 * log(primerConc./b)))) + (16.6 * log10(salt)) - 273.15)) * 1/2 ;
end
tm = [basic; saltadj; tm]';

% compute thermo values using Nearest Neighbor methods
function [NN NNdelta] = near_neigh(seq, seq_length, selfCompFlag, nFlag)

if (~nFlag)
    % nearest neighbor parameters from: Panjkovich and Melo, Bioinformatics  Vol 21 no 6 pp 711-722 2004 [1]
    % rows corresponds to A,C,G,T respectively; columns correspond to A,C,G,T respectively
    Bres86_H = [-9.1,-6.5,-7.8,-8.6,;-5.8,-11,-11.9,-7.8,;-5.6,-11.1,-11,-6.5,;-6,-5.6,-5.8,-9.1,];
    Bres86_S = [-24,-17.3,-20.8,-23.9,;-12.9,-26.6,-27.8,-20.8,;-13.5,-26.7,-26.6,-17.3,;-16.9,-13.5,-12.9,-24,];
   
    ind = sub2indFast([4 4],seq(1:seq_length-1),seq(2:seq_length));
else
    % nearest neighbor parameters as in [1] with added average values for
    % all possible combinations involving 'N'
    % rows corresponds to A,C,G,T,N respectively; columns correspond to
    % A,C,G,T,N respectively
    Bres86_H = [-9.1,-6.5,-7.8,-8.6,-8;-5.8,-11,-11.9,-7.8,-9.125; -5.6,-11.1,-11,-6.5,-8.55;-6,-5.6,-5.8,-9.1,-6.625;-6.625,-8.55,-9.125,-8,-8.075];
    Bres86_S = [-24,-17.3,-20.8,-23.9,-21.5;-12.9,-26.6,-27.8,-20.8,-22.025;-13.5,-26.7,-26.6,-17.3,-21.025;-16.9,-13.5,-12.9,-24,-16.825;-16.825,-21.025,-22.025,-21.5,-20.3438];
  
    seq(seq==15)=5; % substitute numeric value of 'N' with 5
    ind = sub2ind([5 5],seq(1:seq_length-1),seq(2:seq_length));
end

% NN is 4x2 matrix. Columns are DeltaH and DeltaS. Rows correspond to
% methods by Bres86, SantaLucia96, SantaLucia98 and Sugimoto96.
NN = [sum(Bres86_H(ind)),sum(Bres86_S(ind)); ...
    ];

% Corrections: all AT pairs, any GC pairs, symmetry, initiation
if(~nFlag) % ambiguous symbols 'N' not present

    % only AT pairs or any GC pairs?
    if(all((seq == 1)|(seq == 4)))
        NN =  NN + [0 -20.13;];
    else
        NN = NN +  [0 -16.77; ];
    end

    % symmetry
    if(selfCompFlag)
        NN = NN + [0 -1.34 ;  ];
    end

      NNdelta=zeros(1,2);

else % 'N' symbols are present

   error('Ns not allowed in Temp calculation'); 
end

% function optical_density to be added
% function complexity to be added
