%% Build RNA probes
libFolder = SetFigureSavePath('U:\Data\Oligos\LibraryDemo\car2\','makeDir',true);
saveFolder = SetFigureSavePath([libFolder,'vRNA\'],'makeDir',true);
minProbes = 1;

% a fasta file containing the target RNA (full mRNA, or intron) 
%    file may contain many different targets
geneTemplateFasta = 'U:\Data\Oligos\Library19\car2\car2_RNAtarget.fasta';
libTag = 'car2RNA';

startReadout = 97; % we typically shift the starting readout to avoid overlap and potential crosstalk with DNA probes
rptFasta = 'U:\GenomeData\Mouse\RNA\mm_rRNA_and_tRNA.fasta';
readoutSeqs = fastaread('U:\Data\Oligos\CommonOligos\v5_Adapters1_readouts001-384_230406.fasta');
readoutSeqsLocal = readoutSeqs(startReadout:end);  % which readouts use
fwdPrimerIndex = 5;
revPrimerIndex = 19+fwdPrimerIndex;  
    
libName = [saveFolder,libTag,'_probes.fasta'];
refGenome = 'U:\GenomeData\GenomeAssemblies\mm10\mm10.fasta';
allRNAfasta = 'U:\GenomeData\Mouse\mRNA.fasta';

 %% Build probes   
    [probeFasta,targetRegions,probeSets] = FastaToSeqProbes(...
        geneTemplateFasta,saveFolder,...
        'fwdPrimerIndex',fwdPrimerIndex,...
        'revPrimerIndex',revPrimerIndex,...
        'repeatFasta',rptFasta,...
        'locusGeneName',libTag,...
        'secondaries',readoutSeqsLocal,...
        'minProbesPerRegion',minProbes,...
        'maxProbesPerRegion',200,...
        'truncateUniform',true);  % 
    

    maxWordSize = 25; % for BLAST
    maxNumHits = 50; % for BLAST 
    probeFilteredFasta = probeFasta;
    [blastSummary,blastData,allHits] = BLAST(probeFasta,allRNAfasta,'numThreads',40,'wordSize',maxWordSize,'verbose',false); % 
    keepProbes = allHits <= maxNumHits; 
    removedReads = {probeFasta(~keepProbes).Header};
    nameparts = cellfun(@(x) strsplit(x,'__'),removedReads,'UniformOutput',false);
    nameparts = cat(1,nameparts{:});
    probeFilteredFasta(~keepProbes) = [];
    disp(['removed ',num2str( sum(~keepProbes) ), ' probes based on BLAST hits to ',allRNAfasta]);
    
    WriteFasta(libName,probeFilteredFasta,[],'Append',false);
    save([saveFolder,libTag,'_targetRegions.mat'],'targetRegions','probeSets','probeFasta','probeFilteredFasta');
    disp(['wrote ',saveFolder,libTag,'_targetRegions.mat']);
    
    %% read fasta and create Bed from BLAST
    libProbes = fastaread(libName);
    disp(['constructed ', num2str(length(libProbes)),' probes']);    
    ProbeFastaToBed(libName,'blastDatabase',refGenome);    




%-------------------------------------------------------------------------
%% Archival
%-------------------------------------------------------------------------
copyfile( [mfilename('fullpath'),'.m'],[saveFolder,mfilename,'.m']);
disp('------------------------------------------------------------------');
cprintf([0 0 1],'Analysis Complete');
disp(['Copied analysis script to ' saveFolder,mfilename,'.m']);
