function [probeFasta,targetRegions,probeSets,insufficientProbe,cuts] = FastaToSeqProbes(regionFasta,~,varargin)
% takes a fasta file which is already chunked into blocks (e.g. separate
% chromsome segments for ORCA, or separate genes for sequential RNA FISH)
% and designs a tile of probes for each block. Probes are selected to meet
% some design constraints (length, TmRange, specRange, off-target
% threshold, threePrimeSpace etc). 
% Outputs a list of ready-to-order probes: Fwd and rev primers, RT primer,
% barcode sequence, and optimized tile of length L probes. 
%
% Updates
% changed secondary source file to user Adapters1 file
% 
% Install Notes
% 
% currReadouts = 'C:\Data\Oligos\CommonOligos\Adapters1_readouts001-384.fasta' % introduced ~6/30/20 
% currReadouts = 'C:\Data\Oligos\CommonOligos\v2_Adapters1_readouts001-384.fasta';  % introduced 9/30/20 
% currReadouts = 'C:\Data\Oligos\CommonOligos\v5_Adapters1_readouts001-384_230406.fasta'; % introduced 04/0/23
currReadouts = 'U:\Data\Oligos\CommonOligos\v6_Adapters1_readouts001-384_240212.fasta'; % introduced 04/0/23

defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryVerbose', 'boolean', false}; 
defaults(end+1,:) = {'fwdPrimers', 'fasta', fastaread('U:\Data\Oligos\CommonOligos\FwdPrimers.fasta')}; % FwdPrimers
defaults(end+1,:) = {'revPrimers', 'fasta', fastaread('U:\Data\Oligos\CommonOligos\RevPrimers.fasta')}; % RevPrimers
defaults(end+1,:) = {'secondaries', 'fasta', fastaread(currReadouts)}; % ReadoutSeqs
defaults(end+1,:) = {'commonRT', 'string', 'catcaacgccacgatcagct'};  % 20 bp complimentary to P4-405 end. 
defaults(end+1,:) = {'repeatFasta', 'string', ''}; 
defaults(end+1,:) = {'genomeOffTarget','string',''};
defaults(end+1,:) = {'startPrimer','integer',1};
defaults(end+1,:) = {'fwdPrimerIndex','integer',1};
defaults(end+1,:) = {'revPrimerIndex','integer',[]};
defaults(end+1,:) = {'offsetPrimerPairs','integer',0};
defaults(end+1,:) = {'offsetSecondaries','integer',0};
defaults(end+1,:) = {'minProbesPerRegion','integer',20};
defaults(end+1,:) = {'maxProbesPerRegion','integer',inf};
defaults(end+1,:) = {'truncateUniform','boolean',false};
defaults(end+1,:) = {'probeLength','integer',40};
defaults(end+1,:) = {'TmRange','array',[65, 90]};
defaults(end+1,:) = {'gcRange','array',[.2, .8]};
defaults(end+1,:) = {'specRange','array',[.9, 1]}; % this locus vs. others in the library. troublesome when multiple files are in the same fasta. 
defaults(end+1,:) = {'reverseTile','boolean',false};
defaults(end+1,:) = {'locusGeneName','string',''};
defaults(end+1,:) = {'numCompete','integer',1}; % obsolete.
defaults(end+1,:) = {'repeatOTrange','array', [-1 0]};  % [-1 0])
defaults(end+1,:) = {'secOTrange','array', [-1 0]};  % [-1 0])
defaults(end+1,:) = {'genOTrange','array', [-1 30]};  % [-1 0])
defaults(end+1,:) = {'repeatOTmatch','integer', 14};  % 
defaults(end+1,:) = {'secOTmatch','integer', 12};  % 
defaults(end+1,:) = {'specMatch','integer', 15};  %
defaults(end+1,:) = {'genomeOTmatch','integer',17}; % bp   genOTrange = [-1 100] % genome OT range
defaults(end+1,:) = {'readNames','array',{}}; %
defaults(end+1,:) = {'parallel','integer',1};
defaults(end+1,:) = {'noOverlap','boolean',true}; % minimize overlap based on threePrime space (if false returns all probes)
defaults(end+1,:) = {'threePrimeSpace','float',0}; % add extra space between probes or remove space between probes. 
defaults(end+1,:) = {'showCuts','boolean',false}; % For troubleshooting probe loss.
defaults(end+1,:) = {'fiducialNoTruncate','integer',0}; % this readout number will not be truncated to the max probes (typically =1 if the fiducial is in the first probe) 
defaults(end+1,:) = {'doubleReadout','boolean',false}; % tandem copies of the barcode/secondary sequence
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
%% Create genomeObjet object
% here we are using the Transcriptome class to specify blocks of DNA
% sequence instead of blocks of RNA sequence. 
regionFastaData = fastaread(regionFasta);

% hack the Transcriptome object that requires gene='' in the name
fastaHeader = {regionFastaData.Header}';
nEntries = length(fastaHeader);
if ~contains(fastaHeader{1},'gene=')
    for i=1:nEntries
        currName = regionFastaData(i).Header;
        regionFastaData(i).Header = ['',' gene=',currName];
        maxLength = parameters.probeLength*parameters.maxProbesPerRegion*10;  % just a compute failsafe  (was 2). truncation to max probes happens later
        if length(regionFastaData(i).Sequence) > maxLength
           regionFastaData(i).Sequence = regionFastaData(i).Sequence(1:maxLength); 
        end
    end
end


try
    genomeObject = Transcriptome(regionFastaData,'verbose',false);
catch er
    cprintf([1 0 0],'Error calling "Transcriptome". This class is copyright protected. It may be downloaded here: https://github.com/ZhuangLab/MERFISH_analysis. Redistribution is Prohibited.');
    error(er.getReport);
end
    
%% Create specificity table
% this will get used later to keep only probes that have high specifity to 
% the target fasta entry and not to other fasta entries in the regionFasta
% passed to the function.  

if parameters.parallel > 1
     p = gcp('nocreate');
   if isempty(p)
      parpool(parameters.parallel); 
      p = gcp('nocreate');
   elseif parameters.parallel ~= p.NumWorkers
      delete(gcp('nocreate')); 
      parpool(parameters.parallel); 
      p = gcp('nocreate');
   end
else
    p = [];
end

specTable = OTTable(genomeObject, parameters.specMatch);
selfTables = cell(nEntries,1);
for n=1:nEntries
    selfTables{n} =  OTTable(regionFastaData(n), parameters.specMatch);
end
              
%% Create off target tables
% creates an off-target tabel using my library of secondary sequences 
% and a second table using the repeat fasta table passed as full-file-path
% to a fasta file. 

otTables = []; 
otTableNames = {}; 
otTableRanges = {}; 

% specify additional genomic regions to be used for off-target
if ~isempty(parameters.genomeOffTarget)
    offTargetFasta = OTTable(fastaread(parameters.genomeOffTarget),...
                              parameters.genomeOTmatch, 'verbose', false,...
                              'parallel',p); 
    otTables = [otTables,offTargetFasta];
    otTableNames = [otTableNames,'offTargetFasta'];
    otTableRanges = [otTableRanges,parameters.genOTrange]; % number of ot hits to allow
end


otSecondaries = OTTable(parameters.secondaries, parameters.secOTmatch, 'verbose', false,'parallel',p); % ,'parallel',p
otTables = [otTables,otSecondaries];
otTableNames = [otTableNames,'otSecondaries'];
otTableRanges = [otTableRanges,parameters.secOTrange]; % number of ot hits to allow

if ~isempty(parameters.repeatFasta)
    repeatRegions = OTTable(fastaread(parameters.repeatFasta), parameters.repeatOTmatch, 'verbose', false,'parallel',p); 
    otTables = [otTables,repeatRegions];
    otTableNames = [otTableNames,'repeatRegions'];
    otTableRanges = [otTableRanges,parameters.repeatOTrange]; % number of ot hits to allow  
end

%% ------------------------------------------------------------------------
% Build TRDesigner
%%-------------------------------------------------------------------------
%% Create Target Region Designer object

clear gtrDesigner;
gtrDesigner = TRDesigner('transcriptome', genomeObject, ...
        'specificityTable', specTable, ...
        'selfTables',selfTables,...
        'OTTables', otTables, ...
        'OTTableNames', otTableNames,...
        'verbose', false,...
        'parallel',p);

    %%
  
    
%% build probes
otSpecs = [otTableNames; otTableRanges];
otSpecs = otSpecs(:)';

[targetRegions,indsToKeep] = gtrDesigner.DesignTargetRegions(...
                        'regionLength', parameters.probeLength, ...
                        'GC', parameters.gcRange, ...
                        'Tm', parameters.TmRange, ...
                        'specificity', parameters.specRange,...
                        'OTTables',otSpecs,...
                        'threePrimeSpace',parameters.threePrimeSpace,...
                        'noOverlap',parameters.noOverlap);

                    display(['... completed: ' datestr(now)]);

insufficientProbe = [targetRegions.numRegions] < parameters.minProbesPerRegion;  



%%
nP = min([25,length(regionFastaData)]);
if parameters.showCuts
    figure(3); clf; figure(4); clf;
    tm = gtrDesigner.GetRegionTm(parameters.probeLength);
    gc = gtrDesigner.GetRegionGC(parameters.probeLength);
    kept = indsToKeep;
    pens = gtrDesigner.penalties; % off targetFasta
    for i=1:nP
        numCuts = length(pens)+5;
        numBases = length(regionFastaData(i).Sequence);
        cuts = zeros(numCuts,numBases);
        p=0;
        for p=1:length(pens)
          cuts(p,1:length(pens{p}{i})) = pens{p}{i}>=1;
        end
        cuts(p+1,1:length(gtrDesigner.specificity{i})) = (gtrDesigner.specificity{i}<min(parameters.specRange));
        cuts(p+2,1:length(gc{i})) = gc{i}<min(parameters.gcRange);
        cuts(p+3,1:length(gc{i})) = gc{i}>max(parameters.gcRange);
        cuts(p+4,1:length(tm{i})) = tm{i}<min(parameters.TmRange);
        cuts(p+5,1:length(tm{i})) = tm{i}>max(parameters.TmRange);  

        figure(3); subplot(nP/5,5,i); 
        for c=1:numCuts
            plot(.5*c+c+cuts(c,:),'.'); hold on;
        end
        plot(kept{i});
        title(1-sum(sum(cuts,1)>0)/numBases);

        figure(4);  subplot(nP/5,5,i); 
        plot(sum(cuts,1)); hold on; ylabel('excluded regions');
        if ~isempty(targetRegions(i).startPos)
            plot(targetRegions(i).startPos,3,'k.');
        end
        plot(kept{i}); title(['frac bases kept=',num2str(sum(kept{i})/numBases)]);
    end
    figure(3); 
    legend('non-spec','gc-min','gc-max','tm-min','tm-max','kept')
    % legend(otTableNames{:},'non-spec','gc-min','gc-max','tm-min','tm-max','kept')
else
    cuts = [];
end
%%

% Remove targets with insufficient probes (no use waisting stains). 
if sum(insufficientProbe) > 0
    removedData = targetRegions(insufficientProbe);
else
    removedData = [];
end
    
if parameters.verbose && sum(insufficientProbe) > 0
    disp(['removed ',num2str(sum(insufficientProbe)),' probes due to low coverage']);
 
    % disp(removedGeneNames');
end


% truncate targets with too many probes
%  also requires us to create a new cell array of targetSeqs, since the
%  targetRegions(r).sequence property is read-only
numDomains = length(targetRegions); % total number of regions to build probes for.
targetSeqs = cell(numDomains,1);
targetNumRegions = zeros(numDomains,1);

for r= 1:numDomains
    targetSeqs{r} = targetRegions(r).sequence';
    maxProbes = parameters.maxProbesPerRegion;
    if targetRegions(r).numRegions > parameters.maxProbesPerRegion
        if parameters.veryVerbose
           disp(['truncating extra probes from region ',...
               targetRegions(r).geneName]);
        end
        if r ~= parameters.fiducialNoTruncate
            if parameters.truncateUniform
                nSeqs = length(targetSeqs{r});
                selSeqs = sort(randperm(nSeqs,maxProbes));
                targetSeqs{r} = targetSeqs{r}(selSeqs);
            else 
                targetSeqs{r} = targetSeqs{r}(1:maxProbes);
            end
        end
    end
    targetNumRegions(r) = length(targetSeqs{r});
end



try
    % profile of a representative probe set
    c = 1;
    segLength = length(regionFastaData(c).Sequence);
    sampleFig = figure(2); clf;
    subplot(2,2,1); hist(targetRegions(c).startPos,.25E3:.5E3:segLength-.25E3); xlim([-250,segLength+250]);
    title(['probes per 500bp, ',targetRegions(c).geneName]);
    subplot(2,2,2); plot(targetRegions(c).Tm); title('Tm distribution');
    subplot(2,2,3); plot(targetRegions(c).GC); title('GC distribution');
    % subplot(2,2,3); plot(targetRegions(c).penalties); title('penality profile');
    subplot(2,2,4); plot(targetRegions(c).specificity); title(['Specificity profile N=',num2str(targetNumRegions(c))]);
%    SaveFigure(sampleFig,'name',['sampleFig_',targetRegions(c).geneName],'formats',{'png'},'overwrite',true);

    % profile averaged across all probe sets. 
    segLength = mean(cellfun(@length,{regionFastaData.Sequence}));
    statsFig = figure(1); clf;
    y = cellfun(@(x) hist(x,.25E3:.5E3:segLength-.25E3),{targetRegions.startPos},'UniformOutput',false);
    subplot(2,2,2); hist([targetRegions.Tm]); title('Ave. Tm distribution');
    subplot(2,2,3); hist([targetRegions.GC]); title('Ave. GC distribution');
    maxProbes = max(cellfun(@length,{regionFastaData.Sequence}))/parameters.probeLength;
    subplot(2,2,4); hist(targetNumRegions,linspace(0,maxProbes,50)); title('probes per region');

    subplot(2,2,1);
    bar(targetNumRegions,'b'); hold on;
    littleProbe = [targetRegions.numRegions];
    littleProbe(~insufficientProbe) = NaN;
    bar(littleProbe,'r'); xlim([0,length(targetRegions)]);
    title('probes per region');
    SaveFigure(statsFig,'name',['statsFig_',parameters.locusGeneName],'formats',{'png'},'overwrite',true);
catch er
    if parameters.verbose
        disp('encountered problem rendering or saving figures');
        disp(er.message)
    end
end

% actually remove stuff
targetSeqs(insufficientProbe) = [];
targetRegions(insufficientProbe) = [];
numDomains = length(targetRegions); % total number of regions to build probes for.
%%  Append secondaries and primers

readoutSeqs = parameters.secondaries(parameters.offsetSecondaries+1:end);


% Description of probes:  Design 1: 
%                  30 nt           30 nt
% 3'            P4-Alexa405, A647-readout-n, 
% 5' FwdIndex, cy3-commonRT, TruncSecondaryBarcode, targetRegion, RevIndex
%       20 nt,     20 nt,         20 nt,             40 nt,       20 nt,   total: 120 

d = parameters.startPrimer - 1 + parameters.fwdPrimerIndex; % primer counter
if ~isempty(parameters.revPrimerIndex)
    o = parameters.revPrimerIndex -d; 
else
    o = parameters.offsetPrimerPairs; % shuffle primer pairs
end
probeSeqs =  cell(numDomains,1);
probeNames = cell(numDomains,1);
probeSets(1).target = '';
probeSets(1).readout = '';
probeSets(1).numReadouts = 0;
s = 0; % secondary counter
try
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
    targetingSeqs = cellfun(@seqrcomplement,targetSeqs{r},'UniformOutput',false);
    % combine parts
    if ~isempty(parameters.commonRT)
        common = lower(parameters.commonRT);
        fwd = upper(parameters.fwdPrimers(d).Sequence);
    else
        common = '';
        fwd = lower(parameters.fwdPrimers(d).Sequence);
    end
    readSeq = upper( seqrcomplement(readoutSeqs(r).Sequence(1:20)) );
    if parameters.doubleReadout
        readSeq = [readSeq,readSeq];
    end

    probeSeqs{c} = strcat(...
        fwd,...
        common,...
        readSeq,...  Sequence(end-19:end) this is causing confusion will change 
        cellfun(@lower,targetingSeqs,'UniformOutput',false),...
        upper(seqrcomplement(parameters.revPrimers(d+o).Sequence))  );
    if isempty(parameters.readNames)
        readName = ''; 
    else
        readName = ['-',parameters.readNames{k}];
    end
    id = targetRegions(c).id; % this property is read-only
    if isempty(id)
         id = '0'; % avoid creating __ from _0_ if id is empty 
    end
    if ~isempty(parameters.commonRT)
        commonTag = ['__','commonRT'];
    else
        commonTag = '';
    end
    namebits = strsplit(readoutSeqs(r).Header,'_');
    barcodeName = namebits{1};
    probeNames{c} =  strcat(...
        [parameters.fwdPrimers(d).Header,...
        commonTag,...
        '__',barcodeName,readName,...
        '__',targetRegions(c).geneName,...
        '_',id],...
        '_p',cellstr(num2str([1:length(targetSeqs{r})]','%04d')),...     '_p',cellstr(num2str([1:targetRegions(c).numRegions]','%04d'))  
        '__',parameters.revPrimers(d+o).Header) ;
    probeSets(c).target = targetRegions(c).geneName;
    probeSets(c).readout = readoutSeqs(r).Header;
    probeSets(c).numReadouts = length(targetSeqs{r});
    
    if sum(cellfun(@isempty,targetSeqs{r})) > 0
        error('sequence missing')
    end
end

catch er
   warning(er.getReport);
   disp('put debug point here');
end

if parameters.verbose
    disp(['generated ',num2str(length(cat(1,probeSeqs{:}))),' oligos for ',num2str(numDomains),' domains.']);
end

if sum(cellfun(@isempty,probeSeqs))>0
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

