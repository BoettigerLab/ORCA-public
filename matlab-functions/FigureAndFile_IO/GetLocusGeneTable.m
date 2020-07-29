function geneTable = GetLocusGeneTable(locusname,varargin)
% Returns a gene table for the specified locus.
% 
% Note:
% All transcript files must be text files, and must contain columns
% chrom, strand, txStart, txEnd, 'symbol' (
% (essentially bed format).
% The columns can be in any order and any delimeter may be used. 

global dm3 hg19 mm9;

defaults = cell(0,3);
defaults(end+1,:) = {'genome','string','dm3'};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'hg19_transcripts','string','F:\GenomeData\GenomeAssemblies\hg19\hg19_GencodeGenesV24.txt'};
defaults(end+1,:) = {'mm9_transcripts','string','F:\GenomeData\GenomeAssemblies\mm9\mm9_genes.txt'};
defaults(end+1,:) = {'dm3_transcripts','string','C:\Data\Fly\dm3\dm3_NascentRNA.txt'};
defaults(end+1,:) = {'hg19_fasta','string','F:\GenomeData\GenomeAssemblies\hg19\hg19.fasta'};
defaults(end+1,:) = {'mm9_fasta','string','F:\GenomeData\GenomeAssemblies\mm9\mm9.fasta'};
defaults(end+1,:) = {'dm3_fasta','string','C:\Data\Fly\dm3\dm3.fasta'};
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
elseif strcmp(pars.genome,'mm9')
    genome = mm9;
    nascentRNA = pars.mm9_transcripts;
    if isempty(mm9)
        mm9 = fastaread(pars.mm9_fasta);
    end
else
    disp(['genome ',pars.genome,' is not currently supported']);
end

nascentRNAs = readtable(nascentRNA); 


%% Process Locus 
% locusname = 'chr3LHet:600000-900000'

% Get correct chromosome
[chr,locusStart,locusEnd] = ParseLocusName(locusname);
fastaIdx = StringFind({genome.Header}',chr,'exactly',true);
if isempty(fastaIdx)
    error(['did not find ',chr ' in genome fasta file ']);
end
chrIdx = strcmp(nascentRNAs.chrom,chr);

% grab correct region in chromosome
inRange = nascentRNAs.txEnd >= locusStart & nascentRNAs.txStart <= locusEnd;
strand = nascentRNAs.strand(chrIdx & inRange);
chrom =  nascentRNAs.chrom(chrIdx & inRange);
names =  nascentRNAs.Symbol(chrIdx & inRange);   % may need some genome specific adaption if different names are used for the common name. 
% names = regexprep(names,{'Dmir\','Dvir\'},{'',''}); % addressing an name bug
geneStarts = nascentRNAs.txStart(chrIdx & inRange) - locusStart + 1;
geneEnds = nascentRNAs.txEnd(chrIdx & inRange) - locusStart + 1;

% if the gene extends past the region, let's just flip the part we span
geneStarts(geneStarts<1) = 1;
geneEnds(geneEnds>(locusEnd-locusStart+1)) = locusEnd-locusStart+1; 

geneTable = table(names,chrom,geneStarts,geneEnds,strand);

