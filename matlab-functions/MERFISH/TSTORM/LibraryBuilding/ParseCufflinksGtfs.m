function parameters = ParseCufflinksGtfs(cufflinksOutput,varargin)
% ParseCufflinksGtfs(data.fpkm_tracking)
% Writes a fastafile including the sequence of the the most common isoform
% for each gene expressed above a min expression threshold.  

% -------------------------------------------------------------------------
% Global variables
% -------------------------------------------------------------------------
global gencode
global HumanGenome

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'minExpression','nonnegative',1E-5};
defaults(end+1,:) = {'saveFasta','string',''};
defaults(end+1,:) = {'writeFasta','boolean',true};
defaults(end+1,:) = {'HumanChrFastaPath','string','D:\Data\Genomics\HumanGenome\ChrFastaEnsemblRelease73\'};
defaults(end+1,:) = {'gencodeAnnotations','string','D:\Data\Genomics\HumanGenome\Annotations\gencode.v18.annotation.gtf'};


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'A MList is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);




if isempty(parameters.saveFasta)
    dataPath = fileparts(cufflinksOutput);
    fastaOut = [dataPath,'IMR90transcriptome.fasta'];
else
    fastaOut = parameters.saveFasta;
end

%% Check that the system can find the filepaths

if ~exist(parameters.HumanChrFastaPath,'dir')
    error(['Could not connect to ',parameters.HumanChrFastaPath,'.  ',...
   'Make sure you have access to this directory or pass ParseCufflinksGtfs ',...
   'a different "HumanChrFastaPath" folder containing fasta files for all chromosomes']);
end
if ~exist(parameters.gencodeAnnotations,'file')
    error(['Could not access ',parameters.gencodeAnnotations,'.  ',...
   'Make sure you have access to this directory or pass ParseCufflinksGtfs ',...
   'a different "gencodeAnnotations" folder containing gtf annotations for all genes']);
end


% -------- Load Cufflinks Data -------------------------
disp('Reading Cufflinks Gene Data...'); 
fid = fopen(cufflinksOutput);
    headertext = textscan(fid,'%s',13,'delimiter','\t');
    fmt = repmat('%s ',1,13);
    cufflinksData = textscan(fid,fmt,'CollectOutput',true,'TreatAsEmpty','-');
fclose(fid);
cufflinksData = cufflinksData{1};
%-----------------------------------------------------


% 'tracking_id'     1
% 'class_code'      2
% 'nearest_ref_id'  3
% 'gene_id'         4
% 'gene_short_name' 5
% 'tss_id'          6
% 'locus'           7
% 'length'          8
% 'coverage'        9
% 'FPKM'            10
% 'FPKM_conf_lo'    11
% 'FPKM_conf_hi'    12
% 'FPKM_status'     13


%% Check Data
iFPKM = str2double(cufflinksData(:,10));
iGeneNames = cufflinksData(:,5);
% Check for HiDATA.  Should be empty
HIDATAidx = find(strcmp(cufflinksData(:,13),'HIDATA'));
FAILidx = find(strcmp(cufflinksData(:,13),'FAIL'));
if ~isempty(HIDATAidx)
    warning(['some genes too highly expressed to determine FPKM  ',...
        'try increasing cufflinks-flag --max-bundle-frags <int>']);
    
    unique(iGeneNames(HIDATAidx))
    
elseif ~isempty(FAILidx)
     warning(['some genes FAIL  ',...
        'try increasing cufflinks-flag --max-bundle-frags <int>']);
    unique(iGeneNames(FAILidx))
   
else 
    disp('Cufflinks succesfully found FPKMs for all genes');
end
iFPKM(FAILidx) = 2E6;
iFPKM(HIDATAidx) = 2E6;

%% Load Human Genome Annotation Data
if isempty(gencode) % ~exist('gencode','var')
    disp('Reading GenCode GTF file...'); 
    fid = fopen(parameters.gencodeAnnotations);
    fmt = repmat('%s ',1,9);
    gencode = textscan(fid,fmt,'delimiter','\t','HeaderLines',5);
    fclose(fid);
    disp('GenCode loading complete');
else
    disp('gencode data already loaded'); 
end

% Parse gneome annotation columns into more obvious names

% column 1 is chromosome
% column 3 is description (exon, utr, cds)
% column 4 is seq start
% column 5 is seq stop
% column 7 indicates strand
% column 9 is a long string with gene names and transcript IDs and more

featureGene = strcmp(gencode{3},'gene');
featureTranscript= strcmp(gencode{3},'transcript');
featureExon = strcmp(gencode{3},'exon');
featureUTR = strcmp(gencode{3},'UTR');
featureCDS = strcmp(gencode{3},'CDS');
featureSelenocysteine = strcmp(gencode{3},'Selenocysteine');
featureStartCodon = strcmp(gencode{3},'start_codon');
featureStopCodon = strcmp(gencode{3},'stop_codon');

% retrieve only the exons
gIsoID = gencode{9}(featureExon);
gChr = gencode{1}(featureExon);
gFeature = gencode{3}(featureExon);
gStart= gencode{4}(featureExon);
gEnd= gencode{5}(featureExon);
gStrand= gencode{7}(featureExon);


%% Histogram gene expression

figure(2); clf; hist(log10(iFPKM(iFPKM>parameters.minExpression)),100);
xlabel('log_1_0(FPKM)');  xlim([-4,6])


%%
uniqueGenes = unique(cufflinksData(iFPKM>parameters.minExpression,4)); % .01 is 10,000 fold less thn GAPDH
length(uniqueGenes);

% Check a few control genes

% Find GAPDH
GAPDHidx = find(strcmp(iGeneNames,'GAPDH'));
disp(['GAPDH FPKM: ',num2str(max(iFPKM(GAPDHidx)),4)]);

% Find Act
ActBidx = find(strcmp(iGeneNames,'ACTB'));
disp(['ActB FPKM: ',num2str(max(iFPKM(ActBidx)),4)]);

% Max expressed gene
[maxFPKM, maxFPKMidx] = max(iFPKM);
disp(['most expressed gene: ',cufflinksData{maxFPKMidx,5},' fpkm=',num2str(maxFPKM)]);


% goOn = input('continue?')
% if goOn == 0
%     return
% end

%%

% Use only substanitively expressed genes: 
%   ( > 1 in 10000 relative to highest expressed)
expressed = iFPKM>parameters.minExpression;
eFPKM = iFPKM(expressed);
eData = cufflinksData(expressed,:);

% Find unique genes, compute number of isoforms per gene
geneIDs = cufflinksData(expressed,4);
[uniqueIDs,idxIn,idxOut] = unique(geneIDs);
numGenes = length(uniqueIDs);
numIsoforms = hist(idxOut,1:numGenes);
% eData = eData(idxIn,:); % all the other data


% The top c most expressed genes account for as much expression as the
% remainder of the genome.  What is c?
eFPKMsort = sort(eFPKM);

c = 0;
while  sum(eFPKMsort(end-c:end)) < sum(eFPKMsort(1:end-(c+1)));
    c = c+1;
end
disp(['top ',num2str(c),' genes account for >50% of total expression']); 

%% Load Human Genome into Memory;
if isempty(HumanGenome) && parameters.writeFasta % ~exist('HumanGenome','var'
    HumanChrFastas = dir([parameters.HumanChrFastaPath,'*.fasta']);
    numChrs = length(HumanChrFastas);
    HumanGenome = cell(numChrs,2);
    disp('loading human genome fasta...'); 
    for n=1:numChrs; % n=25
        HumanGenome{n,1} = regexprep(HumanChrFastas(n).name,'.fasta','');
        HumanGenome{n,2} = fastaread([parameters.HumanChrFastaPath,HumanChrFastas(n).name]);
        if parameters.verbose
            disp(['loaded chromosome ',HumanGenome{n,1}]); 
        end
    end
    vardata = whos('HumanGenome');
    disp([num2str(vardata.bytes/10E9,3),' Gb of data loaded']);
end

%%

% whos data gencode
clear data gencode; 

disp('Finding most highly expressed isoform and its sequence'); 

try
    currFasta = fastaread(fastaOut);
    fastaStart = length(currFasta); 
catch er
    disp(er.message); 
    fastaStart = 1;
end

GeneData.commonNames{numGenes} = '';
GeneData.geneID{numGenes} = '';
GeneData.isoFPKM{numGenes} = [];
GeneData.totFPKM{numGenes} = [];
GeneData.maxFPKM{numGenes} = []; 
GeneData.isoID{numGenes} = [];
GeneData.locusTxt{numGenes} = '';


% Loop through unique genes, keep isoform with 
for n=fastaStart:numGenes % n =750
    geneN = idxOut==n; 
   % geneN = ~cellfun(@isempty,strfind(eData(:,5),'GAPDH')); %  Jump to a particular isoform    
    IsoformFPKMs = eFPKM(geneN); % fast
    geneIdxStart = find(geneN,1,'first');
    [maxFPKM,maxIsoIdx] = max(IsoformFPKMs); % ? slow
    totFPKM = sum(IsoformFPKMs); 
    isoIdx = geneIdxStart-1+maxIsoIdx;
    
   % isoIdx = find(~cellfun(@isempty,strfind(eData(:,1),'ENST00000396859.1')))   % GAPDH
   % isoIdx = find(~cellfun(@isempty,strfind(eData(:,1),'ENST00000331789.5')))   % ACTB   
        
    commonName = eData{isoIdx,5};
    isoFPKM = eData{isoIdx,10};  % should be equla to maxFPKM
    geneID = eData{isoIdx,4};
    isoID = eData{isoIdx,1};
    locusTxt = eData{isoIdx,7};
    chrEnd = strfind(locusTxt,':');
    locusSep = strfind(locusTxt,'-');
    chrName = locusTxt(1:chrEnd-1);
%     locusStart = str2double(locusTxt(chrEnd+1:locusSep-1));
%     locusEnd = str2double(locusTxt(locusSep+1:end));
    
%     maxIsoHeader = ['IsoformID:_',isoID,...
%         '__GeneID:_',geneID,...
%         '__FPKM:_',num2str(totFPKM),...
%         '__name:_',commonName,...
%         '__locus:_',locusTxt];
    
    GeneData.commonNames{n} = commonName;
    GeneData.geneID{n} = geneID;
    GeneData.isoFPKM{n} = str2double(isoFPKM);
    GeneData.totFPKM{n} = totFPKM;
    GeneData.maxFPKM{n} = maxFPKM; 
    GeneData.isoID{n} = isoID;
    GeneData.locusTxt{n} = locusTxt;
    
    if parameters.writeFasta
            maxIsoHeader = ['IsoformID_',regexprep(isoID,'\.','P'),....
            '__FPKM_',num2str(totFPKM)];

        % Find the exon coordinates corrsponding to this splice form
        %     and the strand direction
        gIdx = find(~cellfun(@isempty, strfind(gIsoID, isoID)));
        pStrand = strcmp(gStrand(gIdx(1)),'+');
        txStarts = gStart(gIdx);
        txEnds = gEnd(gIdx); 
        clear v Seq;
        v='';
        if pStrand
            for i=1:length(gIdx)
               v = [v, txStarts{i},':',txEnds{i},',']; %#ok<AGROW>
            end
        else
           for i=1:length(gIdx)
               j = length(gIdx)-i+1;
               v = [v, txStarts{j},':',txEnds{j},',']; %#ok<AGROW>
           end
        end
        v = ['[',v(1:end-1),']'];

    %  % Troubleshooting
    %    coords = [str2double(txStarts),str2double(txEnds),NaN*ones(length(txEnds),1)]';
    %    figure(2); clf; plot(coords(:),1:length(coords(:)));
    %     
    %     gChr(gIdx)
    %     gFeature(gIdx)
    %     gStart(gIdx)
    %     gEnd(gIdx)
    %     gStrand(gIdx)

    % The first exon of plus strands;  
    % fastaIdx = find(strcmp(HumanGenome(:,1),chrName));
    % i = 1; Seq = HumanGenome{fastaIdx,2}.Sequence(eval([txStarts{i},':',txEnds{i}]));

    % The last exon of minus strands;  
    % fastaIdx = find(strcmp(HumanGenome(:,1),chrName));
    %  i = 1; Seq = seqrcomplement(HumanGenome{fastaIdx,2}.Sequence(eval([txStarts{i},':',txEnds{i}])))

    % Get the correct sequences from the corresponding chromsome FASTA
       fastaIdx = find(strcmp(HumanGenome(:,1),chrName));
       if pStrand
            Seq = HumanGenome{fastaIdx,2}.Sequence(eval(v));
       else
            Seq = seqrcomplement(HumanGenome{fastaIdx,2}.Sequence(eval(v)));
       end
            WriteFasta(fastaOut,maxIsoHeader,Seq,'Append',true);
    end
    disp([num2str(n/numGenes*100,3),'% complete']);
end

geneDataFile = regexprep(fastaOut,'.fasta','.mat');
save(geneDataFile,'GeneData'); 













