%% Design library against select mouse RNAs
% clear all; startup;
% Updates
%  built from TzuChiao
%  Updated to remove BLAST, use all_mRNA from mm9 as otFasta
clc;
global hg38;
libTag = 'oroRNA';
tableRNA = readtable(['C:\Data\Oligos\Library09\TonyOroJill\RNA\','RNA.xlsx']); % ,'Sheet','Actual_List'
saveFolder = SetFigureSavePath('C:\Data\Oligos\Library09\TonyOroJill\RNA\v1','makeDir',true);
refGenome = 'U:\GenomeData\GenomeAssemblies\hg38\hg38.fasta';
abundGenes = 'U:\GenomeData\Human\hek_genes_top_percentile_and_rRNA.fasta';  % abund genes to avoid
% BLAST parameters

blastDatabase =  'U:\GenomeData\GenomeAssemblies\mm9\mrna.fasta';

% blastDatabase =  'U:\GenomeData\GenomeAssemblies\hg38\gencode_v25_pc_transcripts.fasta';
% BuildBLASTlib(blastDatabase); % fails Unsupported ID type ENST00000335137.3 ??  

maxWordSize = 25; % for BLAST
maxNumHits = 100; % for BLAST 

startReadout = 101;
fwdPrimerIndex = 13;
revPrimerIndex = 5;
minProbesPerRegion = 50;
maxProbesPerRegion = 160;
threePrimeSpace = -25;
specRange = [0.9 1];
commonRT = '';

libName = [saveFolder,libTag,'_','probes.fasta'];
templateFastaFile =  [saveFolder,libTag,'_template.fasta'];

selectGenes = tableRNA.GeneName;
% % ---- deal with extra stuff in hox names I copied from UCSC:
% nameParts = cellfun(@(x) strsplit(x,' '),selectGenes,'UniformOutput',false)
% allNames = cat(1,nameParts{:})
% selectGenes = allNames(:,1);
% ------

% start readoutSeqs on plate 2. 
readoutSeqs = fastaread('C:\Data\Oligos\CommonOligos\ReadoutSeqs.fasta');
readoutSeqs = readoutSeqs(startReadout:end);


gencodeTable = readtable('U:\GenomeData\GenomeAssemblies\hg38\hg38_gencodeV25.txt'); 
txTable = readtable('U:\GenomeData\GenomeAssemblies\hg38\hg38_transcriptTable.txt'); % from hg38 annotations > knownGene.txt.gz
txTable.Properties.VariableNames = {'name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','jnk','ucid'};

if exist(templateFastaFile,'file')~=0
    delete(templateFastaFile);
end

% load hg38
if isempty(hg38)
   hg38 = fastaread(refGenome);
end           
genomebuild = hg38;

% clean up gene names from crazy hg38 table dump from UCSC
gene_names = cellfun(@(x) strsplit(x,';'),gencodeTable.name(:),'UniformOutput',false);
tx_names = cellfun(@(x) x{2},gene_names,'UniformOutput',false);
tx_names = cellfun(@(x) x(16:end),tx_names,'UniformOutput',false);
tx_names = cellfun(@(x) regexprep(x,'"',''),tx_names,'UniformOutput',false);
tx_simple = cellfun(@(x) x(1:15),tx_names,'UniformOutput',false);
gene_names = cellfun(@(x) x{5},gene_names,'UniformOutput',false);
gene_names = cellfun(@(x) x(12:end),gene_names,'UniformOutput',false);
gene_names = cellfun(@(x) regexprep(x,'"',''),gene_names,'UniformOutput',false);

% % match tables (for future)
% ids = StringFind(tx_names,txTable{:,1},'exactly',true); 


geneTable = txTable;
% Select a unique isoform for each gene
%   this version choses the form with the most exons. 
try
    ids = StringFind(gene_names,selectGenes,'exactly',true);
    
    if iscell(ids) % if ids are not unique (multiple isoforms exist)
        uids = zeros(length(selectGenes),1);
        for g=1:length(ids) % chose isoform with most exons
            nIsos = length(ids{g});
            isoID = ids{g};
            txID = StringFind(geneTable.name,tx_simple(isoID));
            nIsos = length(txID);
            if nIsos == 0
                disp(['no match found for ',selectGenes{g}]);
                disp(['perhaps there is a different gene name?']);
            else
                exonCnt = zeros(nIsos,1);
                for i=1:nIsos
                    exonCnt(i) = geneTable.exonCount(txID(i));
                end
                [~,m] = max(exonCnt);
                uids(g) = ids{g}(m);
            end
        end
    else
        uids = ids;
    end
    uids(uids==0) = []; % remove the missing genes
    geneTable = geneTable(uids,:);
    geneTable.genID = geneTable.name;
    geneTable.name = selectGenes; % 
catch er
    disp(er.getReport);
    disp('troubleshoot');
end


IntronExonTableToFasta(geneTable,templateFastaFile,'genome',genomebuild)
templateFasta = fastaread(templateFastaFile);
nTargets  = length(templateFasta);
rnaName = cell(nTargets,1);
rnaLocus = cell(nTargets,1);
readoutName = cell(nTargets,1);
for n=1:nTargets
    namebits = strsplit(templateFasta(n).Header,' ');
    rnaName{n} = regexprep(namebits{2},'gene=','');
    rnaLocus{n} = namebits{1};
    readoutName{n} = readoutSeqs(n).Header;
end
rnaTable = table(rnaName,rnaLocus,readoutName);
disp(rnaTable);

writetable(rnaTable,[saveFolder,libTag,'_rnaTable.xlsx']);


%% Build probes

%   'genomeOffTarget',otFastaFile,... 
%   'genomeOTmatch',18,...
%   'genOTrange',[-1 20],...      

[probeFasta,targetRegions,probeSets] = FastaToSeqProbes(...
    templateFastaFile,saveFolder,...
    'fwdPrimerIndex',fwdPrimerIndex,...
    'revPrimerIndex',revPrimerIndex,...
    'locusGeneName',libTag,...
    'secondaries',readoutSeqs,...
    'commonRT',commonRT,...
    'minProbesPerRegion',minProbesPerRegion,...
    'maxProbesPerRegion',maxProbesPerRegion,...
    'threePrimeSpace',threePrimeSpace,...
    'truncateUniform',false,...
    'repeatFasta',abundGenes,...
    'specRange',specRange);  % 

probeFilteredFasta = probeFasta;
[blastSummary,blastData,allHits] = BLAST(probeFasta,blastDatabase,'numThreads',40,'wordSize',maxWordSize,'verbose',false); % 
keepProbes = allHits <= maxNumHits; 
removedReads = {probeFasta(~keepProbes).Header};
nameparts = cellfun(@(x) strsplit(x,'__'),removedReads,'UniformOutput',false);
nameparts = cat(1,nameparts{:});
%[n,v] = occurrences(nameparts(:,3))
probeFilteredFasta(~keepProbes) = [];
disp(['removed ',num2str( sum(~keepProbes) ), ' probes based on BLAST hits to ',blastDatabase]);

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

