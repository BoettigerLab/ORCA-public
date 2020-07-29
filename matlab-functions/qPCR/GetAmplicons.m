
database = 'D:\Data\Genomics\DmelGenome\BLASTlib\Dmel_Genome.fasta';
genomeFasta = fastaread(database);



seqF = 'ACATTTGTTTGGGTCGAAGC';
chrHitsPlus = cell(length(genomeFasta),1);
chrHitsMinus = cell(length(genomeFasta),1);
for i=1:length(genomeFasta)
    chrHitsPlus{i} = strfind(genomeFasta(i).Sequence,primerPairs{p,1});
    chrHitsMinus{i} = strfind(genomeFasta(i).Sequence,seqrcomplement(primerPairs{p,2}));
end
chrPlus = ~cellfun(@isempty,chrHitsPlus);
startPlus = chrHitsPlus{c};

chrMinus = ~cellfun(@isempty,chrHitsPlus);
startMinus = chrHitsPlus{c};

