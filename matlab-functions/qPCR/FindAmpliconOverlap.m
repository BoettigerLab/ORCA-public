function [overlap,errorPairs] = FindAmpliconOverlap(primerPairs,varargin)
%
% currently assumes primers are unique in genome
% assumes forward and reverse primers are on the same chromosome

genomeDefault =  'D:\Data\Genomics\DmelGenome\BLASTlib\Dmel_Genome.fasta';

defaults = cell(0,3);
defaults(end+1,:) = {'database', 'string',genomeDefault};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% main function
% parameters.database = genomeDefault;
genomeFasta = fastaread(parameters.database);
numChrs = length(genomeFasta);

%%
numPairs = size(primerPairs,1);
startPlus = NaN(numPairs,1);
startMinus = NaN(numPairs,1);
chrPlus = NaN(numPairs,1);
chrMinus = NaN(numPairs,1);
amplicons = cell(numPairs,1);
errorPairs = zeros(numPairs,1);
for p=1:numPairs; % p=29
    chrHitsFwd = cell(numChrs,1);
    chrHitsRev = cell(numChrs,1);
    for i=1:length(genomeFasta)
        chrHitsFwd{i} = [strfind(genomeFasta(i).Sequence,primerPairs{p,1}),  strfind(genomeFasta(i).Sequence,seqrcomplement(primerPairs{p,1}))];
        chrHitsRev{i} = [strfind(genomeFasta(i).Sequence,seqrcomplement(primerPairs{p,2})), strfind(genomeFasta(i).Sequence,primerPairs{p,2})];
    end
    try
    chrFwdHit = find(~cellfun(@isempty,chrHitsFwd));
    chrRevHit = find(~cellfun(@isempty,chrHitsRev));
    if length(chrFwdHit)>1 && length(chrRevHit)==1
        chrFwdHit = chrRevHit;
    elseif length(chrFwdHit)==1 && length(chrRevHit)>1
        chrRevHit = chrFwdHit;
    end
    chrPlus(p) = chrFwdHit;
    chrMinus(p) = chrRevHit;
    
    seqStart =  chrHitsFwd{chrPlus(p)};
    seqEnd = chrHitsRev{chrMinus(p)}+length(primerPairs{p,2});
    if length(seqStart)>1 && length(seqEnd)==1
        [~,isclose] = min( (seqStart-seqEnd).^2);
        seqStart = seqStart(isclose);
    elseif length(seqStart)==1 && length(seqEnd)>1
        [~,isclose] = min( (seqEnd-seqStart).^2);
        seqEnd = seqEnd(isclose);
    end
    
    startPlus(p) = seqStart;   
    startMinus(p) = seqEnd;
    amplicons{p} = genomeFasta(chrPlus(p)).Sequence(startPlus(p):startMinus(p));
    % cprintf([0 .5 0],['sucessfuly mapped sequence ',num2str(p),'\n']);
    catch
         cprintf('err',['error with sequence ',num2str(p), ' nonexistent or non-unique  \n']);
    end

        if length( find(~cellfun(@isempty,chrHitsFwd))) > 1
            cprintf('err',['multiple forward primer hits ','with sequence ',num2str(p), '\n']);
            errorPairs(p) = p;
        elseif length( find(~cellfun(@isempty,chrHitsFwd))) < 1
            cprintf('err',['no forward primer hits ','with sequence ',num2str(p), '\n']); 
            errorPairs(p) = p;
        elseif length(chrHitsFwd{chrFwdHit}) > 1
            cprintf('err',['multiple forward primer hits ','with sequence ',num2str(p), ' \n']);
            errorPairs(p) = p;
        elseif length(chrHitsFwd{chrFwdHit}) < 1
            cprintf('err',['no forward primer hits ','with sequence ',num2str(p), '\n']); 
            errorPairs(p) = p;
        end
        if length(find(~cellfun(@isempty,chrHitsRev))) > 1
            cprintf('err',['multiple reverse primer hits ','with sequence ',num2str(p), '\n']);
            errorPairs(p) = p;
        elseif length(find(~cellfun(@isempty,chrHitsRev))) < 1
            cprintf('err',['no reverse primer hits ','with sequence ',num2str(p), '\n']);
            errorPairs(p) = p;
        elseif length(chrHitsRev{chrRevHit}) > 1
            cprintf('err',['multiple reverse primer hits ','with sequence ',num2str(p), '\n']); 
            errorPairs(p) = p;
        elseif length(chrHitsRev{chrRevHit}) < 1
            cprintf('err',['no reverse primer hits ','with sequence ',num2str(p), '\n']); 
            errorPairs(p) = p;
        end        
end
errorPairs = nonzeros(errorPairs); 
%%
overlap = false(numPairs,1);
overlapBases = cell(numPairs,1);
overlapChr = cell(numPairs,1);
for c=1:numChrs
    onChr = find( chrPlus == c );
    coords = zeros(length(genomeFasta(c).Sequence),1);
    for p=1:length(onChr)
        coords(startPlus(onChr(p)):startMinus(onChr(p))) = ...
        coords(startPlus(onChr(p)):startMinus(onChr(p))) +1;
        if any(coords(startPlus(onChr(p)):startMinus(onChr(p)))>1)
           overlap(onChr(p)) = true; 
           overlapBases{onChr(p)} = find(coords(startPlus(onChr(p)):startMinus(onChr(p)))>1);
           overlapChr{onChr(p)} = genomeFasta(c).Header;
           disp( [startPlus(onChr),startMinus(onChr),onChr] );
           disp(length(overlapBases{onChr(p)}));
        end
    end
end