% Library Design
startup;

%--------------------------------------------------------------------------%
% Required Files
% 
% BLAST database -- % should contain all expressed RNAs 
% rRNA and tRNA database (these highly expressed RNA species should be
% treated separately with more stringent filters).  


%% Step 1: Build Fasta file containing only the most common isoform for all
% RNA (avoids double entries).  

% folder = 'D:\Data\Genomics\TSTORM\HumanGenome\IMR90\mthread_IMR90_totRNA_rep2\';
folder = 'D:\Data\Genomics\TSTORM\HumanGenome\IMR90\mthread_IMR90_totRNA_rep1\';


cufflinksOutput = [folder,'isoforms.fpkm_tracking'];
saveFasta = [folder,'TotRNArep1_totFPKM.fasta'];
parameters = ParseCufflinksGtfs(cufflinksOutput,'saveFasta',saveFasta);

%% Step 2: Build BLAST database
BuildBLASTlib(saveFasta,'legacy',true); 


%% Step 3: Select genes based on FPKMs
geneList = fastaread(saveFasta);
geneList(1:5).Header

FPKMstart = cellfun(@(x) strfind(x,'FPKM'),{geneList.Header},'UniformOutput',false);
FPKMend = cellfun(@(x) strfind(x,'name'),{geneList.Header},'UniformOutput',false);
FPKMs = cellfun(@(x,y,z) str2double(z(x+6:y-3)),FPKMstart,FPKMend,{geneList.Header}); % ,'UniformOutput',false);

highGenes = (FPKMs < 100E3 & FPKMs > 10E3);
sum(highGenes)

gapdhRange = (FPKMs < 10E3 & FPKMs > 1E3);
sum(gapdhRange);
geneList(gapdhRange).Header;

someGAPDHlike = {geneList(gapdhRange).Header};
someGAPDHlike(1:15)'


midLowGenes = (FPKMs < 100 & FPKMs > 10);
sum(midLowGenes);

HT = geneList(midLowGenes);

GAPDHfastaIdx = find(~cellfun(@isempty, strfind({geneList.Header},'GAPDH_')));

FPKMs(GAPDHfastaIdx)
geneList(GAPDHfastaIdx).Header
geneList(GAPDHfastaIdx).Sequence