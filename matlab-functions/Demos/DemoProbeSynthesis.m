% 
% This script builds probes.
% 
% CC-BY Alistair Boettiger 06/29/2020

regionFasta = '.\segmentedProbes.fasta'; % this is the path to your regionFastaFile.
%  it should be a fasta file and look something like this:
% 
%  > myRegion_barcode001
%  GTGTCTGTAATCACAGTGGTGAGGGGGTAAAGGGTGGAGTTTATCCTGCTAAAACTGAAAGACACGGGGC
%  CTGGGAAGTCCCATGCCTACGAACCAGGAGAGATGTTTTTCAAGTCTCTAACGCCTTTCACAATCTGCTC
%  ... # 2-10 kb of genome sequence
% 
%  > myRegion_barcode002
%  GTCACCCAGCTGACCTCACACACTGCGGTGACTCACCCACTTCCAAAACACTTCTTAGGACCACAAAGGG
%  GATGAAAGCTCATCTCTGAAAAACTGTAAGGTTTTTTTGCTGCTAAACTTGTCTCAAATAGGCCACTATT
%  ... # 2-10 kb of genome sequence
%
%  ...
%  > myRegion_barcode120 

% You may wish to update the following parameters (such as gcRange or
% maxProbesPerRegion).  
[probeFasta,targetRegions,probeSets,removedData] = FastaToSeqProbes(...
         regionFasta,[],...
        'fwdPrimers','.\FwdPrimers.fasta',...
        'revPrimers','.\RevPrimers.fasta',...
        'secondaries','.\Secondaries.fasta',...
        'repeatFasta','.\Dmel_repeats.fasta',...
        'repeatOTmatch',13,...  % max number of bases which can match to sequences in the repeatFasta.
        'minProbesPerRegion',20,...  
        'maxProbesPerRegion',100,...
        'gcRange',[.3,.7],...
        'truncateUniform',true,...
        'threePrimeSpace',1); 
    
% save the final list of probes as a fasta file.  This fasta file can be
% sent to a commerical oligopool producer such as CustomArray (now
% GenScript) or Twist Biosciences, for synthesis. 
faName ='./probesOut.fasta'; 
WriteFasta(faName,probeFasta,[],'Append',false);


ProbeFastaToBed(faName,'blastDatabase',refLibrary);
