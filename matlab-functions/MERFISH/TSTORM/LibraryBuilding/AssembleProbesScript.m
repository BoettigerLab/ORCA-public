% L5 construction

% Probe Design
%  * Add a common RT primer 

% To add:
%  * individual genes in 2 colors, all the same bases (including secondary
%  binding sites) as in the version of the gene embedded in the coded
%  library.  Use for colocalization.  




% some file-path info
lib5Prep = 'D:\Data\Genomics\TSTORM\Library5_Prep\';
oligoFolder = [lib5Prep,'ProbesT70X72S76\'];% 
libFolder = 'C:\Users\Alistair\Documents\Research\Projects\TSTORM\Lib\';
lib5Folder = [libFolder,'Lib5\'];
primerFasta = [libFolder,'Lib5primers.fasta'];
cntrlFasta = [TSTORMlib,'BLASTedRandomSeqs.fasta'];

if ~exist(lib5Folder,'dir')
    mkdir(lib5Folder);
end

% Load primers and chose forward and reverse primers randomly
indexPrimers = fastaread(primerFasta);

% Chose a universal and remove from list
universal = indexPrimers(1).Sequence;
indexPrimers = indexPrimers(2:end); 

% Load probe data from processed OligoArray output
load([oligoFolder,'ProbeData.mat'],'ProbeData'); 

% Load cntrl Data
cntrl = fastaread(cntrlFasta);

% Load secondary sequences
secondaries = fastaread([lib5Folder 'secondaries.fasta']); %These are the secondary sequences

%% Experiment 1-8
% Non-stringent 30mers, 198 probes per gene, all shuffles
subLibNum = 1; % starting sublibrary

% Shuffles
% 1. original
% 2. draw new primary sequences, keep codebook, keep secondaries, keep secondary-primary map. 
% 3. draw new secondary sequences, 
% 4. shuffle codewords, keep primaries, keep secondaries, keep secondary-primary map.
% 5. shuffle secondary-primary map, keep codewords, keep primaries.
% 6. keep everything, just change index primers
% 7. shuffle everything
% 8. shuffle everything

% Gene Properties
numGenes = 12; % Number of genes in each sublibrary
numCntrls = 1; % number of intentionally unexpressed genes to include.  
numBlanks = 1; % number of codewords to leave blank
numOligos = 200; % Minimum number of oligos targeting each gene
minFPKM = 1; % Minimum FPKM
maxFPKM = 10E3; % Maximum FPKM
binsFPKM = 6; % number of bins between min and max FPKM in which we want genes
requiredGenes = {'THBS1'};

% Codebook Properties
numLetters = 8; %Total number of bits for SECDED codewords    <-- This should be computed from the numGenes and num onBits 
numDataBits = 4; %Number of data bits
onBits = 4; %Number of used bits (number of colocalization required)
bitsPerProbe = 2; % number of secondaries attached to each targeting sequence

% A little rounding so the tile size distribution works out even
numOligos = floor(numOligos/nchoosek(onBits,bitsPerProbe))*nchoosek(onBits,bitsPerProbe);

% Generate SECDED codebook with indicated # of ON bits
usedSECDEDcodewords = GenSECDED(numLetters,numDataBits,onBits);

% Select Genes based on available probe number and FPKM 
pickedGenes = PickGenes(ProbeData, numGenes, numOligos, 'minFPKM', minFPKM, 'maxFPKM', maxFPKM, 'requiredGenes', requiredGenes, 'binnedFPKM', binsFPKM);
pickedGenes = AddCntrlSeqs(pickedGenes,cntrl,numOligos,numCntrls,numBlanks);

shuffleCodewords    = [1 0 0 1 0 0 1 1];
shufflePriSec       = [1 0 0 0 1 0 1 1];
shufflePri          = [1 1 0 0 0 0 1 1] ;
changeSecondaries   = [0 0 1 0 0 0 0 0];

[OligoSeq1, OligoName1]=GenProbe(usedSECDEDcodewords,pickedGenes,secondaries,indexPrimers,numOligos,numGenes+numCntrls,lib5Folder,...
'subLibNum',subLibNum,'universal',universal,'shuffleCodewords',shuffleCodewords,'shufflePriSec',shufflePriSec,'shufflePri',shufflePri,'changeSecondaries',changeSecondaries);


%% Experiment 9-16
subLibNum = 9; % 

% Shuffles
% 1. original
% 2. draw new primary sequences, keep codebook, keep secondaries, keep secondary-primary map. 
% 3. draw new secondary sequences, 
% 4. shuffle codewords, keep primaries, keep secondaries, keep secondary-primary map.
% 5. shuffle secondary-primary map, keep codewords, keep primaries.
% 6. keep everything, just change index primers
% 7. shuffle everything
% 8. shuffle everything

% Gene Properties
requiredGenes = pickedGenes.CommonName(1:12); % Keep all the 14 gene library selections
numGenes = 32; % Number of genes in each sublibrary
numCntrls = 3; % number of intentionally unexpressed genes to include.  
numBlanks = 3; % number of codewords to leave blank
numOligos = 200; % Minimum number of oligos targeting each gene
minFPKM = 1; % Minimum FPKM
maxFPKM = 10E3; % Maximum FPKM

% Codebook Properties
numLetters = 12; %Total number of bits for SECDED codewords   
onBits = 4; %Number of used bits (number of colocalization required)
bitsPerProbe = 2; % number of secondaries attached to each targeting sequence

% A little rounding so the tile size distribution works out even
numOligos = floor(numOligos/nchoosek(onBits,bitsPerProbe))*nchoosek(onBits,bitsPerProbe);

% Generate SECDED codebook with indicated # of ON bits
usedSECDEDcodewords = GenSECDED(numLetters,numDataBits,onBits);

% Select Genes based on available probe number and FPKM 
pickedGenes = PickGenes(ProbeData, numGenes, numOligos, 'minFPKM', minFPKM, 'maxFPKM', maxFPKM, 'requiredGenes', requiredGenes, 'binnedFPKM', binsFPKM);
pickedGenes = AddCntrlSeqs(pickedGenes,cntrl,numOligos,numCntrls,numBlanks);

shuffleCodewords    = [1 0 0 1 0 0 1 1];
shufflePriSec       = [1 0 0 0 1 0 1 1];
shufflePri          = [1 1 0 0 0 0 1 1] ;
changeSecondaries   = [0 0 1 0 0 0 0 0];

[OligoSeq2, OligoName2]=GenProbe(usedSECDEDcodewords,pickedGenes,secondaries,indexPrimers,numOligos,numGenes+numCntrls,lib5Folder,...
'subLibNum',subLibNum,'universal',universal,'shuffleCodewords',shuffleCodewords,'shufflePriSec',shufflePriSec,'shufflePri',shufflePri,'changeSecondaries',changeSecondaries);

%%

allSeqs = cat(1,OligoSeq1{:},OligoSeq2{:});
allNames = cat(1,OligoName1(:),OligoName2(:));
length(allSeqs)




