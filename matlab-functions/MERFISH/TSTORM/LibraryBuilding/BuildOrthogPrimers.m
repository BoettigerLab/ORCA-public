function BuildOrthogPrimers(saveName,varargin)

%% BuildOrthogPrimers
% Function for building orthogonal primers
% Alistair Boettiger 
% boettiger.alistair@gmail.com 
%
% Description
% -------------------------------------------------------------------------
% Uses 240K 'orthogonal oligos' from Elledge lab data
% Turn these into primers (good TM, GC clamp, shorter, no long rpts)
% BLAST against each other with a custom more stringent orthogonality
% criteria.  Remove non-unique sequences.
% BLAST against reference genome.
% No BLAST hits are allowed to have 3p homology/alignment, regardless of
% the length of sequence aligned.  


%% 
% saveName = 'C:\Users\Alistair\Documents\Research\Projects\TSTORM\Lib\Lib5primers.fasta';

primer_path = fileparts(saveName);
if ~exist(primer_path)
mkdir(primer_path)
end

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------

elledge = [ 'D:\Data\Genomics\OligoPaintProbes\PrimersEtc\Barcode\Barcode\','bc25mer240k.fasta'];
defaults = cell(0,3);
defaults(end+1,:) = {'maxPrimerHomology','positive', 10};
defaults(end+1,:) = {'maxBLASThomology', 'positive', 13};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'maxDimerLength','positive', 6};
defaults(end+1,:) = {'maxHairPinBase', 'positive', 6};
defaults(end+1,:) = {'tMin', 'positive', 74};
defaults(end+1,:) = {'tMax', 'positive', 82};
defaults(end+1,:) = {'threePrime', 'nonnegative', 0};
defaults(end+1,:) = {'primerLength', 'positive', 20};
defaults(end+1,:) = {'BLASTlibs', 'cell', {}};
defaults(end+1,:) = {'elledge','string',elledge};
defaults(end+1,:) = {'primerRoot','string',''};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'A saveName is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% parameters = ParseVariableArguments('',defaults, mfilename);
% parameters.BLASTlibs = BLASTlibs;
% parameters.maxBLASThomology = [13,13,10];

primers = fastaread(parameters.elledge);

N1 = length(primers);
keptprimers(N1).Sequence = '';
keptprimers(N1).Header = '';
k=0;
for n=1:N1
    % shorten to 'primerLength', check for 1 GCs in last 'threePrime' bps.  
    p=0;
    testseq = primers(n).Sequence(1:parameters.primerLength);
    GC3p = length(regexp(testseq(end-parameters.threePrime:end),'G'))+length(regexp(testseq(end-parameters.threePrime:end),'C'));
    while GC3p<1 && p< 25-parameters.primerLength
        p=p+1;
        testseq = primers(n).Sequence(1+p:parameters.primerLength+p);
        GC3p = length(regexp(testseq(end-parameters.threePrime:end),'G'))+length(regexp(testseq(end-parameters.threePrime:end),'C'));
    end
    % Test repeats
    norpts = isempty(regexpi(testseq,'a{4,}|c{4,}|g{4,}|t{4,}','ONCE')); 
    % Test Tm
    props = oligoprop(testseq,'Dimerlength',parameters.maxDimerLength,'HPBase',parameters.maxHairPinBase);    
    Tm = props.Tm(3); 
    TmPass = Tm > parameters.tMin && Tm < parameters.tMax;
    % keep if it passes everything
    if GC3p >= 1 && norpts && TmPass
        k=k+1;
        keptprimers(k).Sequence = testseq;
        keptprimers(k).Header = ['dm',num2str(parameters.primerLength),'_',num2str(k)];
    end    
    if mod(n,1000) == 0
        disp([num2str(n/N1*100,3),'% complete']);
    end
end
keptprimers = keptprimers(1:k); 
disp(['kept ',num2str(k),' of ',num2str(N1),' primers']); 


fastaSave = regexprep(saveName,'.fasta','_preBLAST.fasta'); 
WriteFasta(fastaSave,keptprimers,'','Append',false,'Warnings',false);
disp(['wrote ',fastaSave]); 

% keptprimers = fastaread(fastaSave); 

% BLAST against genome
for i=1:length(parameters.BLASTlibs)
    primers = keptprimers; 
    [~,blastData,allhits] = BLAST(fastaSave,'dataBase',parameters.BLASTlibs{i},'verbose',false);
    keepIdx = allhits<=  parameters.maxBLASThomology(min(i,length(parameters.maxBLASThomology)));
    kept = sum(keepIdx);
    disp(['kept ',num2str(kept),' of ',num2str(length(primers)),...
        ' primers after BLAST to ', parameters.BLASTlibs{i}]); 
    keptprimers = primers(keepIdx);
    fastaSave = regexprep(saveName,'.fasta',['_BLAST',num2str(i),'.fasta']); 
    WriteFasta(fastaSave,keptprimers,'','Append',false,'Warnings',false);
    disp(['wrote ',fastaSave]); 
end

% BLAST against self and T7

    % Make a BLAST library from current fasta
    WriteFasta(fastaSave,'T7seq','TAATACGACTCACTATAGGG',...
        'Append',true);
    BuildBLASTlib(fastaSave,'legacy',false); 
    
    % BLAST against new lib
    [~,selfBLASTdata,allSelfHits]=BLAST(fastaSave,'dataBase',fastaSave,'exludeFirstHit',true,'verbose',false); 
    keepPrimers = allSelfHits<=parameters.maxPrimerHomology;
    disp(['kept ',num2str(sum(keepPrimers)),' of ',num2str(length(primers)),...
        ' primers after BLAST to ', fastaSave]); 

    primersOut = primers(keepPrimers(1:end-1)); % don't include T7 as a primer
    
% Check for homology in the final five basepairs    
    disp('checking for any 5-prime homology');
    finalFive = cellfun(@(x) x(end-4:end),{primersOut.Sequence},'UniformOutput',false);
    [unique3pEnds, uniqueIDs] = unique(finalFive);
    primersOut = primersOut(uniqueIDs);
    disp(['kept ',num2str(length(primersOut))]); 

    if isempty(parameters.primerRoot)
        parameters.primerRoot = 'primer';
    end
    for k=1:length(primersOut)
        primersOut(k).Header = [parameters.primerRoot,'_',sprintf('%03d',k)];
    end
    WriteFasta(saveName,primersOut,'','Append',false,'Warnings',true);
    disp(['wrote ',saveName]); 

    
    

