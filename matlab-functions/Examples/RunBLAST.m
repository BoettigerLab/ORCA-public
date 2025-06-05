
% If the library hasn't been built, we need to build it first using the
% fasta file 
BuildBLASTlib('U:\GenomeData\GenomeAssemblies\test\dm3.fasta');

% once a library has been built, we can blast sequences against it:
[blastResults,blastData,allHits] = BLAST('CACGTGGAGCAGAAC','U:\GenomeData\GenomeAssemblies\test\dm3.fasta');
blastResults