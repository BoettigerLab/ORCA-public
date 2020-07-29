function otFasta = BuildOffTargetFasta(regions,genomeFasta)
% this functions builds a fasta file which contains the whole genome 
% specified by "genomeFasta" with the portion corresponding to the target
% region(s) "regions" (a cell array of locus coordinates) removed from the
% genome file.  The remaining genome sequence is also reverse complemented
% so that downstream off-target tables will have both strands.
%
% Alistair Boettiger
% 11/06/15
% 
% Example:
% regions = {'chr3R:12481406-12810708'}; 
% genomeFasta = '\\Morgan\MorganData\ChromatinLibraries\GenomeFiles\Dm_Genome_withHet.fasta';
% otFasta = BuildOffTargetFasta(regions,genomeFasta)

genomeFiles = fastaread(genomeFasta);
otFasta = genomeFiles;

% remove probe target regions from probe list
for n=1:length(regions)
    [chr,locusStart,locusEnd] = ParseLocusName(regions{n});
    i = StringFind({genomeFiles.Header}',chr,'exactly',true);
    otFasta(i).Sequence(locusStart:locusEnd) = []; 
end

numChr = length(otFasta);
for c=1:numChr
    otFasta(c).Header = [otFasta(c).Header,'_plus'];
    otFasta(c+numChr).Header  = [otFasta(c).Header,'_minus'];
    otFasta(c+numChr).Sequence  = seqrcomplement(otFasta(c).Sequence);
end
