%The purpose of this short script is to execute all of the steps necessary
%to create probes to order.  Since no specific genes are included, it
%creates probes for all of the genes in this region. 
function HumanProbeGen(seqs, preProcessed)

filename = {'C:/Users/Mirae/Desktop/DataDirectory/Human Table UCSC/human21seqs RefSeq'};
chrom = 'chr21';

if preProcessed 
fastaName = probeCreator(filename, seqs, 'saveName', 'humanFasta', ...
    'chromosome', chrom, 'testFile', true, 'useNameTwo', true);
else
sequences = 'F:/GenomeData/Human/hg19/hg19.fasta';
fastaName = probeCreator(filename, sequences, 'saveName', 'humanFasta', ...
    'chromosome', chrom, 'testFile', true, 'useNameTwo', true);
end

savePath = 'C:/Users/Mirae/Desktop/DataDirectory/Fafsa Files/';
[probes, ~, ~] = TargetFastaToProbeFasta(fastaName, savePath);

timeFormat = 'dd-mmm-yyyy HH-MM-SS';
timeStamp = datestr(now, timeFormat);
fastawrite([savePath, 'humanProbesToOrder', timeStamp, '.fasta'], probes)

end