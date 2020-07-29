function [probeFasta,targetRegions,probeSets] = TargetFastaToProbeFasta(regionFasta,analysisSavePath,varargin)
% 
% Description of probes:  Design 1: 
% 3'            P4/anti-cRT-cy3/488, Readout-A647, 
% 5' FwdIndex, cy3/unlab-commonRT,  ReadoutBinding, targetRegion, RevIndex
%       20-bp,        20 bp,              20 bp,       40 bp,       20 bp,   total: 120 
%      fixed          fixed              fixed       probeLength    fixed  
% 

defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; % example optional par and its default 
defaults(end+1,:) = {'fwdPrimers', 'string', 'C:\Data\Oligos\FwdPrimers_L3.fasta'}; % FwdPrimers
defaults(end+1,:) = {'revPrimers', 'string', 'C:\Data\Oligos\RevPrimers_L2.fasta'}; % RevPrimers
defaults(end+1,:) = {'secondaries', 'string', 'C:\Data\Oligos\ReadoutSeqs.fasta'}; % ReadoutSeqs
defaults(end+1,:) = {'commonRT', 'string', 'catcaacgccacgatcagct'};  % 20 bp complimentary to P4-405 end. 
defaults(end+1,:) = {'genomeFasta', 'string', 'C:\Data\Fly\dm3\dm3.fasta'}; 
defaults(end+1,:) = {'repeatFasta', 'string', 'C:\Data\Fly\Dmel_repeats.fasta'}; 
defaults(end+1,:) = {'genomeOffTarget','boolean',false};
defaults(end+1,:) = {'startPrimer','integer',1};
defaults(end+1,:) = {'fwdPrimerIndex','integer',1};
defaults(end+1,:) = {'offsetPrimerPairs','integer',0};
defaults(end+1,:) = {'minProbesPerRegion','integer',15};
defaults(end+1,:) = {'probeLength','integer',40};
defaults(end+1,:) = {'TmRange','array',[65, 90]};
defaults(end+1,:) = {'gcRange','array',[0, 1]};
defaults(end+1,:) = {'specRange','array',[.9, 1]}; % this locus vs. others in the library. troublesome when multiple files are in the same fasta. 
defaults(end+1,:) = {'reverseTile','boolean',false};
defaults(end+1,:) = {'parallel','integer',1};
defaults(end+1,:) = {'locusGeneName','string',''};
defaults(end+1,:) = {'numCompete','integer',1}; 
defaults(end+1,:) = {'repeatOTrange','array', [-1 0]};  % [-1 0])
defaults(end+1,:) = {'secOTrange','array', [-1 0]};  % [-1 0])
defaults(end+1,:) = {'repeatOTmatch','integer', 14};  % 
defaults(end+1,:) = {'secOTmatch','integer', 12};  % 
defaults(end+1,:) = {'specMatch','integer', 17};  %
defaults(end+1,:) = {'genomeOTmatch','integer',17}; % bp 

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'regionFasta (.fasta files) is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
%% Actual Function
%--------------------------------------------------------------------------
%% build probes

%% read in fasta files passed

fwdPrimers = fastaread(parameters.fwdPrimers);
revPrimers = fastaread(parameters.revPrimers);
readoutSeqs = fastaread(parameters.secondaries);

%% parallel pool
if parameters.parallel ~= 1
    if isempty(gcp('nocreate'))
        p = parpool(parameters.parallel); % not very memory efficient.  Doesn't seem to commmit much CPU for all the RAM this kills
    else
        p = gcp;
    end
end

%% Create genomeObjet object
% here we are using the Transcriptome class to specify blocks of DNA
% sequence instead of blocks of RNA sequence. 
regionFastaData = fastaread(regionFasta);
genomeObject = Transcriptome(regionFastaData);

%% Create specificity table
% this will get used later to keep only probes that have high specifity to 
% the target fasta entry and not to other fasta entries in the regionFasta
% passed to the function.  
if parameters.parallel ~= 1
    specTable = OTTable(genomeObject, parameters.specMatch, ...
        'verbose', true, ...
        'transferAbund',false,...
        'parallel',p);
else
    specTable = OTTable(genomeObject, parameters.specMatch, ...
        'verbose', true, ...
        'transferAbund',false);
end

%% Create off target tables
% creates an off-target tabel using my library of secondary sequences 
% and a second table using the repeat fasta table passed as full-file-path
% to a fasta file. 

otTables = []; 
otTableNames = {}; 
otTableRanges = {}; 

% whole genome off-target analysis is slow and will likely be removed in
% the future. 
if parameters.genomeOffTarget
    offTargetFasta = OTTable(fastaread([analysisSavePath,'offTarget.fasta']),...
                              parameters.genomeOTmatch, 'verbose', true,'parallel',0); 
    otTables = [otTables,offTargetFasta];
    otTableNames = [otTableNames,'offTargetFasta'];
    otTableRanges = [otTableRanges,[-1 100]]; % number of ot hits to allow
end

if parameters.parallel ~= 1
    otSecondaries = OTTable(readoutSeqs, parameters.secOTmatch, 'verbose', true,'parallel',p);
else
    otSecondaries = OTTable(readoutSeqs, parameters.secOTmatch, 'verbose', true);
end
otTables = [otTables,otSecondaries];
otTableNames = [otTableNames,'otSecondaries'];
otTableRanges = [otTableRanges,parameters.secOTrange]; % number of ot hits to allow

if parameters.parallel ~= 1
    repeatRegions = OTTable(fastaread(parameters.repeatFasta), parameters.repeatOTmatch, 'verbose', true,'parallel',p);
else
    repeatRegions = OTTable(fastaread(parameters.repeatFasta), parameters.repeatOTmatch, 'verbose', true);
end
otTables = [otTables,repeatRegions];
otTableNames = [otTableNames,'repeatRegions'];
otTableRanges = [otTableRanges,parameters.repeatOTrange]; % number of ot hits to allow  


%% ------------------------------------------------------------------------
% Build TRDesigner
%%-------------------------------------------------------------------------
%% Create Target Region Designer object

clear gtrDesigner;
gtrDesigner = TRDesigner('transcriptome', genomeObject, ...
        'specificityTable', specTable, ...
        'OTTables', otTables, ...
        'OTTableNames', otTableNames);


%% build probes
localGCRange = parameters.gcRange; %  [0 1];
localTmRange = parameters.TmRange;
specRange = parameters.specRange; % [.5 1];
otSpecs = [otTableNames; otTableRanges];
otSpecs = otSpecs(:)';

targetRegions = gtrDesigner.DesignTargetRegions(...
                        'regionLength', parameters.probeLength, ...
                        'GC', localGCRange, ...
                        'Tm', localTmRange, ...
                        'specificity', specRange,...
                        'OTTables',otSpecs);
                display(['... completed: ' datestr(now)]);

% Remove regions with insufficient probe
insufficientProbe = [targetRegions.numRegions] < parameters.minProbesPerRegion;  
targetRegions(insufficientProbe) = [];
if parameters.verbose
    disp(['removed ',num2str(sum(insufficientProbe)),' probes due to low coverage']);
end

%%  Append secondaries and primers

numDomains = length(targetRegions); % total number of regions to build probes for.

d = parameters.startPrimer - 1 + parameters.fwdPrimerIndex; % primer counter
o = parameters.offsetPrimerPairs; % shuffle primer pairs
probeSeqs = {}; %  cell(numDomains,1);
probeNames = {}; % cell(numDomains,1);
probeSets(1).target = '';
probeSets(1).readout = '';
probeSets(1).numReadouts = 0;
s = 0; % secondary counter
for k=1:numDomains
    s = s+1;
    r = ceil(s/parameters.numCompete);
    if r>length(readoutSeqs)
        warning('Exceeded number of readout seqs available!');
        warning('starting a new probe-set with new index primers');
        r = 1; % restart at first readout seq
        s = 1; % restart counter
        d = d+1; % change primers so these are part of a new probe library
        o = o+1; % 
    end   
    if parameters.reverseTile
        c = numDomains+1-k; % assign tile domains in reverse order, the more interesting features are at the end.  
    else
        c=k;
    end
    targetingSeqs = cellfun(@seqrcomplement,targetRegions(c).sequence','UniformOutput',false);
    % combine parts
    probeSeqs{c} = strcat(...,
        upper(fwdPrimers(d).Sequence),...
        lower(parameters.commonRT),...
        upper( seqrcomplement(readoutSeqs(r).Sequence(end-19:end)) ),...
        cellfun(@lower,targetingSeqs,'UniformOutput',false),...
        upper(seqrcomplement(revPrimers(d+o).Sequence))  ); %#ok<*AGROW>
    probeNames{c} =  strcat([fwdPrimers(d).Header,...
        '__','commonRT',...
        '__',readoutSeqs(r).Header,...
        '__',targetRegions(c).geneName,...
        '_',targetRegions(c).id],...
        '_p',cellstr(num2str([1:targetRegions(c).numRegions]','%03d')),...
        '__',revPrimers(d+o).Header) ;
    probeSets(c).target = targetRegions(c).geneName;
    probeSets(c).readout = readoutSeqs(r).Header;
    probeSets(c).numReadouts = length(targetRegions(c).sequence);
    
    if sum(cellfun(@isempty,targetRegions(c).sequence)) > 0
        error('sequence missing')
    elseif length(targetRegions(c).sequence) < targetRegions(c).numRegions
        error('too many regions')
    end
end

if parameters.verbose
    disp(['generated ',num2str(length(cat(1,probeSeqs{:}))),' oligos for ',num2str(numDomains),' domains.']);
end

if ~(sum(cellfun(@isempty,probeSeqs))>0)
    warning('some probe sets are empty!');
end    

% check probe lengths. 
for c=1:length(probeSeqs)
    if length(probeSeqs{c}) < length(probeNames{c})
        disp(length(probeSeqs{c}));
        disp(length(probeNames{c}));
        error(['unequal lengths at c=',num2str(c)]);
    elseif sum(cellfun(@isempty,probeSeqs{c}))>0
        i = find(cellfun(@isempty,probeSeqs{c}));
        disp(probeNames{c}(i));
        error('isempty');
    end
end

% construct fasta file of probes
probeFasta = cell2struct([cat(1,probeNames{:}), cat(1,probeSeqs{:})]',{'Header','Sequence'});

