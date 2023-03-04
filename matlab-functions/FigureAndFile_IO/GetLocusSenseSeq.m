function [Seq,geneNames,geneStarts,geneEnds,geneStrands] = GetLocusSenseSeq(locusname,varargin)
% Returns the specified genome sequence for the specified
% location (i.e. 'chr2L:455203-465638'), with the sequence flipped to the
% whichever strand is coding (i.e. complement sequence is returned
% where-ever a minus strand transcript occurs).  
% 
% Note:
% All transcript files must be text files, and must contain columns
% chrom, strand, txStart, txEnd
% (essentially bed format).
% The columns can be in any order and any delimeter may be used. 

global dm3 hg19 mm9 mm10 hg38 ce11;

defaults = cell(0,3);
defaults(end+1,:) = {'genome','string','dm3'};
defaults(end+1,:) = {'verbose','boolean',true};
% transcript files. has columns with names: % chrom, txStart, txEnd, name, strand
defaults(end+1,:) = {'hg19_transcripts','string','U:\GenomeData\GenomeAssemblies\hg19\hg19_GencodeGenesV24.txt'};
defaults(end+1,:) = {'hg38_transcripts','string','U:\GenomeData\GenomeAssemblies\hg38\hg38_gencodeV25.txt'};
defaults(end+1,:) = {'mm9_transcripts','string','U:\GenomeData\GenomeAssemblies\mm9\mm9_genes.txt'};
defaults(end+1,:) = {'mm10_transcripts','string','U:\GenomeData\GenomeAssemblies\mm10\mm10_NascentTranscripts.txt'}; % mm10_NascentTranscripts,'U:\GenomeData\GenomeAssemblies\mm10\mm10_NascentTranscripts.txt');
defaults(end+1,:) = {'dm3_transcripts','string','U:\GenomeData\GenomeAssemblies\dm3\dm3_NascentRNA.txt'}; % 'C:\Data\Fly\dm3\dm3_NascentRNA.txt'
defaults(end+1,:) = {'ce11_transcripts','string','U:\GenomeData\GenomeAssemblies\ce11\ce11_NascentRNA.txt'};
% genome files
defaults(end+1,:) = {'hg19_fasta','string','U:\GenomeData\GenomeAssemblies\hg19\hg19.fasta'};
defaults(end+1,:) = {'hg38_fasta','string','U:\GenomeData\GenomeAssemblies\hg38\hg38.fasta'};
defaults(end+1,:) = {'mm9_fasta','string','U:\GenomeData\GenomeAssemblies\mm9\mm9.fasta'};
defaults(end+1,:) = {'mm10_fasta','string','U:\GenomeData\GenomeAssemblies\mm10\mm10.fasta'};
defaults(end+1,:) = {'dm3_fasta','string','U:\GenomeData\GenomeAssemblies\dm3\dm3.fasta'};%  'C:\Data\Fly\dm3\dm3.fasta'
defaults(end+1,:) = {'ce11_fasta','string','U:\GenomeData\GenomeAssemblies\ce11\ce11.fasta'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Load Data

if strcmp(pars.genome,'dm3')
    genome = dm3;
    nascentRNA = pars.dm3_transcripts;
    if isempty(dm3)
        dm3 = fastaread(pars.dm3_fasta);
    end
elseif strcmp(pars.genome,'hg19')
    genome = hg19;
    nascentRNA = pars.hg19_transcripts;
    if isempty(hg19)
        hg19 = fastaread(pars.hg19_fasta);
    end
elseif strcmp(pars.genome,'hg38')
    genome = hg38;
    nascentRNA = pars.hg38_transcripts;
    if isempty(hg38)
        hg38 = fastaread(pars.hg38_fasta);
    end
elseif strcmp(pars.genome,'mm9')
    genome = mm9;
    nascentRNA = pars.mm9_transcripts;
    if isempty(mm9)
        mm9 = fastaread(pars.mm9_fasta);
    end
elseif strcmp(pars.genome,'mm10')
    genome = mm10;
    nascentRNA = pars.mm10_transcripts;
    if isempty(mm10)
        mm10 = fastaread(pars.mm10_fasta);
    end
elseif strcmp(pars.genome,'ce11')
    genome = ce11;
    nascentRNA = pars.ce11_transcripts;
    if isempty(ce11)
        ce11 = fastaread(pars.ce11_fasta);
    end
else
    disp(['genome ',pars.genome,' is not currently supported']);
end

nascentRNAs = readtable(nascentRNA); 


%% Process Locus 
% locusname = 'chr3LHet:600000-900000'

% Get correct chromosome
[chr,locusStart,locusEnd] = ParseLocusName(locusname);
if locusEnd-locusStart < 0
   warning('locus stop is less than locus start, flipping coordinates');
   origStart = locusStart; 
   origEnd = locusEnd;
   locusStart = origEnd;
   locusEnd = origStart;
end

fastaIdx = StringFind({genome.Header}',chr,'exactly',true);
if isempty(fastaIdx)
    error(['did not find ',chr ' in genome fasta file ']);
end
chrIdx = strcmp(nascentRNAs.chrom,chr);



% grab correct region in chromosome
inRange = nascentRNAs.txEnd >= locusStart & nascentRNAs.txStart <= locusEnd;
geneStrands = nascentRNAs.strand(chrIdx & inRange);
nameFields = find(contains(nascentRNAs.Properties.VariableNames,{'Symbol','geneID','geneName','name','names'}));
try
    geneNames = nascentRNAs{chrIdx & inRange,nameFields(1)};
catch
    geneNames = 'did not find names';
end
geneStarts = nascentRNAs.txStart(chrIdx & inRange) - locusStart + 1;
geneEnds = nascentRNAs.txEnd(chrIdx & inRange) - locusStart + 1;

% if the gene extends past the region, let's just flip the part we span
geneStarts(geneStarts<1) = 1;
geneEnds(geneEnds>(locusEnd-locusStart+1)) = locusEnd-locusStart+1; 


try
    Seq0 = genome(fastaIdx).Sequence(locusStart:locusEnd);
catch er
    locusname %#ok<NOPRT>
    fastaIdx %#ok<NOPRT>
    chr %#ok<NOPRT>
    warning(er.message);
    error('failed to find correct chromosome in genome fasta'); 
end
    
Ngenes = length(geneStrands);
k = 0;
Seq = Seq0;  
% Note.  We need two copies to avoid flipping and flipping back as we
% encounter new isoforms.  Hence Seq and Seq1. 
for i = 1:Ngenes % i=2
    if strcmp(geneStrands{i},'-')
        Seq( geneStarts(i):geneEnds(i) )  = seqrcomplement(Seq0( geneStarts(i):geneEnds(i) ));
        k = k+1;
    else
        Seq( geneStarts(i):geneEnds(i) ) =  Seq0(geneStarts(i):geneEnds(i));
    end
end
if pars.verbose
    disp(['flipped ',num2str(k),' of ',num2str(Ngenes),' sequences']);
end
if length(Seq) < 1
    error('length of Seq is less than 1, something is wrong.');
end






