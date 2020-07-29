function [oligoSeq, oligoName, codebookTable] = AssembleProbes(ProbeData,varargin)



% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% Sequences to Assemble
defaults(end+1,:) = {'primerFasta', 'string', ''};
defaults(end+1,:) = {'readoutFasta', 'string', ''};
defaults(end+1,:) = {'cntrlFasta', 'string', ''};
defaults(end+1,:) = {'dataFolder', 'string', ''};
defaults(end+1,:) = {'libSaveFolder', 'string', ''};
% Library properties
defaults(end+1,:) = {'numCntrls', 'integer', 0};  % random sequences for negative controls 
defaults(end+1,:) = {'numBlanks', 'integer', 0}; % codewords left blank as controls
defaults(end+1,:) = {'numOligos', 'integer', 200}; % total oligos per probe
defaults(end+1,:) = {'subLibNum', 'integer', 1}; % starting sublibrary (increase to avoid using the same index primers as for a different experiment); 
defaults(end+1,:) = {'universal', 'boolean', true};

% Codebook Properties
defaults(end+1,:) = {'numTotalBits', 'integer', 16};  % Total number of bits to read out (must have at least this number of secondary sequences 
defaults(end+1,:) = {'numDataBits', 'integer', 11}; % 
defaults(end+1,:) = {'onBits', 'integer', 4}; % Number of bits per gene which will be "ON"
defaults(end+1,:) = {'bitsPerProbe', 'integer', 2}; % bits per targeting sequence

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'required: oligo_folder');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

if ~isempty(parameters.dataFolder)
    parameters.primerFasta = [parameters.dataFolder,'PrimerSeqs.fasta'];
    parameters.readoutFasta = [parameters.dataFolder,'ReadoutSeqs.fasta'];
    parameters.cntrlFasta = [parameters.dataFolder,'BLASTedRandomSeqs.fasta'];   
elseif isempty(parameters.readoutFasta)
    error('either "dataFolder" or "readoutFasta" is required');
else
    parameters.dataFolder = [fileparts(parameters.readoutFasta),filesep];
end

if isempty(parameters.libSaveFolder)
    parameters.libSaveFolder = [dataFolder,'AssembledProbes\'];  %  save folder
end

if ~exist(parameters.libSaveFolder,'dir')
    mkdir(parameters.libSaveFolder);
end

% some book properties
numGenes = length(ProbeData.GeneName);
numWords = numGenes + parameters.numCntrls + parameters.numBlanks; % total number of genes and controls in library


%% Main Function


% Load primers and chose forward and reverse primers randomly
indexPrimers = fastaread(parameters.primerFasta);

% Chose a universal and remove from list
if parameters.universal
    universal = indexPrimers(1).Sequence;
else
    universal = '';
end

% Load random non-homology sequences
cntrl = fastaread(parameters.cntrlFasta);

% Load secondary sequences
readoutSeqs = fastaread(parameters.readoutFasta); %These are the secondary sequences

%% Experiment 1

% A little rounding so the tile size distribution works out even

parameters.numOligos = floor(parameters.numOligos/nchoosek(parameters.onBits,parameters.bitsPerProbe))*nchoosek(parameters.onBits,parameters.bitsPerProbe);


% Generate SECDED codebook with indicated # of ON bits
codebookTable = GenSECDED(parameters.numTotalBits,parameters.numDataBits,parameters.onBits);
[numPosGenes,numHybes] = size(codebookTable);
if numPosGenes < numWords
    error('codebook is not large enough to encode all the indicated genes and controls');
elseif numPosGenes > numWords  
    parameters.numBlanks = numPosGenes - parameters.numCntrls - numGenes;
end
disp(['constructed codebook for up to ',num2str(numPosGenes),' genes in ',num2str(numHybes),' total bits']);



% Add control and Blank sequences
pickedGenes = ProbeData;
pickedGenes = AddCntrlSeqs(pickedGenes,cntrl,...
    parameters.numOligos,parameters.numCntrls,parameters.numBlanks);

[oligoSeq, oligoName] = GenProbe(codebookTable,pickedGenes,readoutSeqs,indexPrimers,...
    parameters.numOligos,numGenes+parameters.numCntrls,parameters.libSaveFolder,...
    'bitsPerProbe',parameters.bitsPerProbe,...
    'subLibNum',parameters.subLibNum,...
    'universal',universal);

%


%%

     
