function [expFasta, libCodeTable,targetRegions] = BuildMERFISHprobes(geneFasta,varargin)
% [expFasta, libCodeTable,targetRegions] = BuildMERFISHprobes(geneFasta,varargin)
%% Inputs
% geneFasta - a fasta file containing the sequence of the genes to target.
%% Outputs
% expFasta - a fasta file containing the sequences of all the probes to be
%            synthesized for the experiment.
% libCodeTable - a xlsx table 
% 
%% Options
% Probe design features 
% defaults(end+1,:) = {'numGenes','positive',12};
% defaults(end+1,:) = {'minProbesPerRegion','positive',100};
% defaults(end+1,:) = {'onBits','positive',4};
% defaults(end+1,:) = {'bitsPerProbe','positive',3};
% defaults(end+1,:) = {'startPrimer','integer',1};
% defaults(end+1,:) = {'fwdPrimerIndex','integer',1};
% defaults(end+1,:) = {'offsetPrimerPairs','integer',0};
% defaults(end+1,:) = {'probeLength','integer',30};
% defaults(end+1,:) = {'TmRange','array',[0,100]}; % [65, 90]
% defaults(end+1,:) = {'gcRange','array',[0, 1]};
% defaults(end+1,:) = {'specRange','array',[.9, 1]}; % this locus vs. others in the library. troublesome when multiple files are in the same fasta. 
% defaults(end+1,:) = {'numCompete','integer',1}; 
% 
% % General parameters
% defaults(end+1,:) = {'verbose', 'boolean', true}; % example optional par and its default 
% defaults(end+1,:) = {'saveFolder','string',''};
% defaults(end+1,:) = {'showPlots','boolean',true};
% defaults(end+1,:) = {'parallel','integer',1}; % Might drop this 
% 
% % Reference sequences
% defaults(end+1,:) = {'commonRT', 'string', 'catcaacgccacgatcagct'};  % 20 bp complimentary to P4-405 end. 
% defaults(end+1,:) = {'fwdPrimers', 'fasta', 'C:\Data\Oligos\FwdPrimers_L3.fasta'}; % FwdPrimers
% defaults(end+1,:) = {'revPrimers', 'fasta', 'C:\Data\Oligos\RevPrimers_L2.fasta'}; % RevPrimers
% defaults(end+1,:) = {'useFwdPrimer','boolean',true};
% 
% defaults(end+1,:) = {'secondaries', 'fasta', 'C:\Data\Oligos\ReadoutSeqs.fasta'}; % ReadoutSeqs
% defaults(end+1,:) = {'secOTmatch','integer', 12};  % 
% defaults(end+1,:) = {'secOTrange','array', [-1 0]};  % [-1 0])
% 
% defaults(end+1,:) = {'transcriptomeFasta', 'string', 'C:\Data\Human\hg38\gencode.v25.pc_transcripts.fa'}; 
% defaults(end+1,:) = {'secOTmatch','integer', 17};  % 
% defaults(end+1,:) = {'secOTrange','array', [-1 10]};  % [-1 0])
% 
% defaults(end+1,:) = {'rrnaFasta','string','C:\Data\Human\all_hg_rRNAs.fasta'};
% defaults(end+1,:) = {'rrnaOTmatch','positive',12};
% defaults(end+1,:) = {'rrnaOTrange','positive',[-1 0]};
% 
% defaults(end+1,:) = {'abundFasta','string','C:\Data\Human\hek_genes_top_percentile.fasta'};
% defaults(end+1,:) = {'abundOTmatch','positive',14};
% defaults(end+1,:) = {'abundOTrange','positive',[-1 0]};

%
%% Notes
% Alistair Boettiger
% June 21 2017
% CC BY 
% This function is dependent on Transcriptome Tools by Jeff Moffitt

%% build probes default parameters
defaults = cell(0,3);
% Probe design features 
defaults(end+1,:) = {'expTag','string',''};
defaults(end+1,:) = {'numGenes','positive',[]};
defaults(end+1,:) = {'minProbesPerRegion','positive',100};
defaults(end+1,:) = {'maxProbesPerRegion','positive',100};
defaults(end+1,:) = {'onBits','positive',4};
defaults(end+1,:) = {'bitsPerProbe','positive',3};
defaults(end+1,:) = {'startPrimer','integer',1};
defaults(end+1,:) = {'fwdPrimerIndex','integer',1};
defaults(end+1,:) = {'offsetPrimerPairs','integer',0};
defaults(end+1,:) = {'probeLength','integer',30};
defaults(end+1,:) = {'TmRange','array',[65,90]}; % [65, 90]
defaults(end+1,:) = {'gcRange','array',[.2, .8]};
defaults(end+1,:) = {'numCompete','integer',1}; 

% General parameters
defaults(end+1,:) = {'verbose', 'boolean', true}; % example optional par and its default 
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'parallel','integer',1}; % Might drop this 

% Reference sequences and Off-target stringencies
defaults(end+1,:) = {'commonRT', 'string', 'catcaacgccacgatcagct'};  % 20 bp complimentary to P4-405 end. 
defaults(end+1,:) = {'fwdPrimers', 'fasta', '..\DemoData\FwdPrimers.fasta'}; % FwdPrimers
defaults(end+1,:) = {'revPrimers', 'fasta', '..\DemoData\RevPrimers.fasta'}; % RevPrimers
defaults(end+1,:) = {'useFwdPrimer','boolean',false};

defaults(end+1,:) = {'specRange','array',[.9, 1]}; % this locus vs. others in the library. troublesome when multiple files are in the same fasta. 
defaults(end+1,:) = {'specMatch','integer', 17};  %

defaults(end+1,:) = {'secondaries', 'fasta', '..\DemoData\ReadoutSeqs.fasta'}; % ReadoutSeqs
defaults(end+1,:) = {'secOTmatch','integer', 12};  % 
defaults(end+1,:) = {'secOTrange','array', [-1 0]};  % [-1 0])

defaults(end+1,:) = {'transcriptomeFasta', 'string', 'F:\GenomeData\Human\hg38\gencode.v25.pc_transcripts.fa'}; 
defaults(end+1,:) = {'secOTmatch','integer', 17};  % 
defaults(end+1,:) = {'secOTrange','array', [-1 10]};  % [-1 0])

defaults(end+1,:) = {'rrnaFasta','string','..\DemoData\all_hg_rRNAs.fasta'};
defaults(end+1,:) = {'rrnaOTmatch','positive',12};
defaults(end+1,:) = {'rrnaOTrange','positive',[-1 0]};

defaults(end+1,:) = {'abundFasta','string','..\DemoData\hek_genes_top_percentile.fasta'};
defaults(end+1,:) = {'abundOTmatch','positive',14};
defaults(end+1,:) = {'abundOTrange','positive',[-1 0]};

parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% parallel pool
% if isempty(gcp('nocreate'))
%     p = parpool(parameters.parallel); % 
% else
%     p = gcp;
% end

%% select savefolder

if isempty(parameters.saveFolder)
   saveFolder = [fileparts(geneFasta),filesep];
else
    saveFolder = parameters.saveFolder;
end

%% Create tObject object
% here we are using the Transcriptome class to specify blocks of DNA
% sequence instead of blocks of RNA sequence. 
geneFastaData = fastaread(geneFasta);
tObject = Transcriptome(geneFastaData);


%% Create specificity table
% this will get used later to keep only probes that have high specifity to 
% the target fasta entry and not to other fasta entries in the regionFasta
% passed to the function.  
specTable = OTTable(tObject, parameters.specMatch, ...
                    'verbose', true, ...
                    'transferAbund',false);
                
                % ,... 'parallel',p

%% Create off target tables
% creates an off-target tabel using my library of secondary sequences 
% and a second table using the repeat fasta table passed as full-file-path
% to a fasta file. 



otTables = []; 
otTableNames = {}; 
otTableRanges = {}; 

otRRNA = OTTable(fastaread(parameters.rrnaFasta), parameters.rrnaOTmatch, 'verbose', true); 
otTables = [otTables,otRRNA];
otTableNames = [otTableNames,'repeatRegions'];
otTableRanges = [otTableRanges,parameters.rrnaOTrange]; % number of ot hits to allow  

otAbund = OTTable(fastaread(parameters.abundFasta), parameters.abundOTmatch, 'verbose', true); 
otTables = [otTables,otAbund];
otTableNames = [otTableNames,'abundRegions'];
otTableRanges = [otTableRanges,parameters.abundOTrange]; % number of ot hits to allow  

otSecondaries = OTTable(fastaread(parameters.secondaries), parameters.secOTmatch, 'verbose', true);  % ,'parallel',p
otTables = [otTables,otSecondaries];
otTableNames = [otTableNames,'otSecondaries'];
otTableRanges = [otTableRanges,parameters.secOTrange]; % number of ot hits to allow



%% Create Target Region Designer object

clear tDesigner;
tDesigner = TRDesigner('transcriptome', tObject, ...
        'specificityTable', specTable, ...
        'OTTables', otTables, ...
        'OTTableNames', otTableNames);


%% build probes
localGCRange = parameters.gcRange; %  [0 1];
localTmRange = parameters.TmRange;
specRange = parameters.specRange; % [.5 1];
otSpecs = [otTableNames; otTableRanges];
otSpecs = otSpecs(:)';

targetRegions = tDesigner.DesignTargetRegions(...
                        'regionLength', parameters.probeLength, ...
                        'GC', localGCRange, ...
                        'Tm', localTmRange, ...
                        'specificity', specRange,...
                        'OTTables',otSpecs);

                    display(['... completed: ' datestr(now)]);

if parameters.showPlots
    figure(1); clf; bar([targetRegions.numRegions]);
    title('numProbesPerRegion');
end
                    
insufficientProbe = [targetRegions.numRegions] < parameters.minProbesPerRegion; 
if sum(insufficientProbe) > 0
    warning('These genes will be discarded due to insufficient probes:')
    disp({geneFastaData(insufficientProbe).Header}');
end
                    
lowProbe = [targetRegions.numRegions] < parameters.maxProbesPerRegion ...
    & ~insufficientProbe; 
if sum(lowProbe) > 0
    warning('These genes have less than desired number of probes:')
    disp({geneFastaData(lowProbe ).Header}');
    disp(num2str([targetRegions(lowProbe).numRegions]));
end
%%

%% Remove genes with insufficient probes

targetRegionsAll = targetRegions; % just for reference / debugging
targetRegions = targetRegions(~insufficientProbe);

%% Build codebook
if isempty(parameters.numGenes)
    numGenes = length(targetRegions);
elseif parameters.numGenes < length(targetRegions)
    warning(['number genes requested (',num2str(parameters.numGenes),')',...
        ' is less than number of targetRegions provided']);
    warning(['Randomly selecting ',num2str(parameters.numGenes),...
        ' of ',num2str(length(targetRegions)),' genes']);
    numGenes = min(parameters.numGenes,length(targetRegions)); %  12;
    targetRegions = targetRegions(randperm(length(targetRegions),numGenes));
    disp('kept');
    disp({targetRegions.geneName}')
else
    numGenes = length(targetRegions);
end

onBits = parameters.onBits; % 4; %Number of used bits (number of colocalization required)
bitsPerProbe = parameters.bitsPerProbe; % 3; % number of secondaries attached to each targeting sequence


% Generate Hamming/SECDED codebook with indicated # of ON bits
% libCodebook = BuildHammingCodebook(38);
libCodebook = BuildHammingCodebook(numGenes);

% figure(1); clf; imagesc(libCodebook);
%% load readout-sequences and primer sequences

fwdPriIndex = parameters.startPrimer;
revPriIndex = parameters.startPrimer+parameters.offsetPrimerPairs;

fwdPrimerSeqs = fastaread(parameters.fwdPrimers);
revPrimerSeqs = fastaread(parameters.revPrimers);
readoutSeqs = fastaread(parameters.secondaries);

if parameters.useFwdPrimer
    fwdPrimer = fwdPrimerSeqs(fwdPriIndex).Sequence;
    fwdPrimerName = fwdPrimerSeqs(fwdPriIndex).Header;
else
   fwdPrimer = '';
   fwdPrimerName ='';
end
    
revPrimer = seqrcomplement(revPrimerSeqs(revPriIndex).Sequence);
revPrimerName = revPrimerSeqs(revPriIndex).Header;

commonRT = parameters.commonRT;


%% Combine targeting seqs, readout seqs, and primers

expFasta = [saveFolder, 'E',num2str(fwdPriIndex),'_',parameters.expTag,'_','oligos.fasta'];

if exist(expFasta,'file')
    delete(expFasta)
end


j=1;
oligoName = {};
oligoSeq = {};
for n=1:numGenes
    
    seq = targetRegions(n).sequence;
    geneOnBits = find(libCodebook(n,:));  % indices of the 1 bits (out of 16 letters)
    secCombos = nchoosek(geneOnBits,bitsPerProbe);   % all possible pairs combinations (since each tile has 1 pair)
    numSecCombos = size(secCombos,1); 
    
    
    % A little rounding so the tile size distribution works out even
    numOligos =  min(length(seq),parameters.maxProbesPerRegion); %  100;
    numOligosUsed = floor(numOligos/nchoosek(onBits,bitsPerProbe))*nchoosek(onBits,bitsPerProbe);

    p = 0;
    for k = 1:numSecCombos
        for i = 1:numOligosUsed/numSecCombos
            p=p+1;
            secIndex = secCombos(k,:);

            readoutSeq1 = [seqrcomplement(readoutSeqs(secIndex(1)).Sequence(end-19:end)), ' '];
            readoutSeq2 =  [seqrcomplement(readoutSeqs(secIndex(2)).Sequence(end-19:end)), ' '];
            readoutName1 = [readoutSeqs(secIndex(1)).Header ' '];
            readoutName2 = [readoutSeqs(secIndex(2)).Header ' '];
            if length(secIndex) > 2
                readoutSeq3 = [seqrcomplement(readoutSeqs(secIndex(3)).Sequence(end-19:end)),' '];
                readoutName3 = [readoutSeqs(secIndex(3)).Header ' '];
            else
                readoutSeq3 = '';
                readoutName3 = '';
            end
            if length(secIndex) > 3
                readoutSeq4 = [seqrcomplement(readoutSeqs(secIndex(4)).Sequence(end-19:end)) ' '];        
                readoutName4 = [readoutSeqs(secIndex(4)).Header ' '];
            else
                readoutSeq4 = '';
                readoutName4 = ''; 
            end
            oligoSeq{j} = ...
                {[fwdPrimer, ' ',...
                  commonRT, ' ',...
                  readoutSeq3,' ',...
                  readoutSeq1,' ',...
                  seqrcomplement(seq{p}), ' ',...
                  readoutSeq2,' ',...
                  readoutSeq4,' ',...
                  revPrimer]};

            oligoName{j} = ...
                  {['E',num2str(fwdPriIndex),'_',...
                    parameters.expTag,' ',...
                    fwdPrimerName,' ',...
                    'commonRT',' ',...
                    readoutName3,' ',...
                    readoutName1,' ',...
                    char(targetRegions(n).geneName) '_',...
                    num2str(targetRegions(n).startPos(p)),'_',... 
                    'probe',num2str(j),' ',...
                    readoutName2,' ',...
                    readoutName4,' ',...
                    revPrimerName]};
            WriteFasta(expFasta,oligoName{j},oligoSeq{j},'Append',true,'Warnings',false); 
            j = j + 1;

        end    
    end
end 

disp(['wrote: ',expFasta]);

% Format and save a library codebook file
geneNames = {targetRegions.geneName}';
[numWords,numBits] = size(libCodebook);
if numWords > numGenes
    geneNames = cat(1,geneNames{1:numGenes}, repmat({'Blank'},numWords - numGenes,1));
end
libCodeTable = table(geneNames,libCodebook);
libCodeTableName = [saveFolder,'E',num2str(fwdPriIndex),'_',parameters.expTag,'_CodeTable.xlsx'];
writetable(libCodeTable,libCodeTableName);
disp(['wrote : ',libCodeTableName]);