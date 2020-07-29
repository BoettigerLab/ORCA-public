
%this function goes through and creates probes for the d3 version of
%drosophila fly genome for specific hox related genes. 
function FlyProbeGen(seqs, preProcessed)

filenames = {'C:/Users/Mirae/Desktop/DataDirectory/Fly Table UCSD/AntpSegment',...
    'C:/Users/Mirae/Desktop/DataDirectory/Fly Table UCSD/abd_BSegment'};
decoder = 'C:/Data/Fly/FlyBase_IDs.txt';
%list of the names of genes to makes probes for 
genes = {'CG11648-RB', 'Abd-B', 'abd-A', 'iab-8', 'iab-4', ...
    'bxd', 'Ubx', 'Antp', 'Scr', 'Dfd', 'pb', 'lab',...
    'alphaTub84B', 'zen', 'zen2', 'ftz', 'bcd', 'Taf1'};
%A cell of booleans depicting whether or not each of the genes about needs
%to be translated.  Allows one ot differentiate between different isoforms
%for example. 
transBool = {false, true, true, true, true, ...
    true, true, true, true, true, true, true, ...
    true, true, true, true, true, true};
%
chrom = 'chr3R';

if preProcessed
    fastaName = probeCreator(filenames, seqs, 'nameDecoder', decoder,...
        'chromosome', chrom, 'targetGenes', genes, 'testFile', true,...
        'saveName', 'flyFasta', 'translateBool', transBool,...
        'lengthTag', 3);
else
    seqloc = 'C:/Data/Fly/dm3/dm3.fasta';
    fastaName = probeCreator(filenames, seqloc, 'nameDecoder', decoder,...
        'chromosome', chrom, 'targetGenes', genes, 'testFile', true, ...
        'saveName', 'flyFasta', 'translateBool', transBool,...
        'lengthTag', 3);
end
%
savePath = 'C:/Users/Mirae/Desktop/DataDirectory/Fafsa Files/';
[probes, ~, ~] = TargetFastaToProbeFasta(fastaName, savePath);
timeFormat = 'dd-mmm-yyyy HH-MM-SS';
timeStamp = datestr(now, timeFormat);
fastawrite([savePath, 'flyProbesToOrder', timeStamp, '.fasta'], probes)
end

