% Example of a simple script desiging probe library

% manuall specify target region here
numRegions = 1;
targetName = 'car2_edit';%  'pcada';
coords = 'chr3:14,753,359-15,400,000'; % 'chr3:14,750,000-15,400,000'; % car2
stepSize = 8e3; % 7e3 for pcada
overwrite = false;
maxSteps = 150;
genomeBuild = 'mm10';

libFolder = SetFigureSavePath('U:\Data\Oligos\LibraryDemo\car2\','makeDir',true);
saveFolder = SetFigureSavePath([libFolder,'v1\'],'makeDir',true); % often we make a few different versions and compare the probe density / coverage achieved 
genomeFolder = 'U:\GenomeData\GenomeAssemblies\';
refLibrary = [genomeFolder,'mm10\mm10.fasta'];

% less common parameters 
fwdPrimerIndex = 1;% (we usually change these later when we combine probesets into an order) 
revPrimerIndex = 2; % 
minProbesPerRegion = 41; % 
maxProbesPerRegion = 150; % 
threePrimeSpace = 1; % space between probes (usually 1 or 0)
genOT = [-1,35]; % 
genomeOTmatch = 16; % 
specRange = [.7 100]; % 70% of probes matching this sequence need to map to the correct barcode.  (can be more than 1.0 since a replicate probe might map twice in the same region) 
specMatch = 16 ; % 16;
repeatOTmatch = 13;%  13;
gcRange =  [.3,.6]; %  [.25 .75]; % 
truncateUniform = true; % if more probes than max, should we distribute probes across the domain (truncateUniform=true) or cluster at the beginning (truncateUniform=false)
poolName = targetName;


% I usually don't change these
otFastaFile = [libFolder,'otSeqs.fasta'];
probeFastaFile = [saveFolder,poolName,'_probe.fasta'];



%% load genomes (this is slow)
% making this a global saves time if your run this repeatedly. 
global mm10; 
if isempty(mm10)
    mm10 = fastaread([genomeFolder,'mm10\mm10.fasta']);
    genomeFasta = mm10;
end

% % some other common genomes you might use
% global hg19; 
% if isempty(hg19)
%     hg19 = fastaread([genomeFolder,'hg19\hg19.fasta']);
%     genomeFasta = hg19;
% end
% global mm9; 
% if isempty(mm9)
%     mm9 = fastaread([genomeFolder,'mm9\mm9.fasta']);
%     genomeFasta = mm9;
% end
% global dm3; 
% if isempty(dm3)
%     dm3 = fastaread([genomeFolder,'dm3\dm3.fasta']);
%     genomeFasta = dm3;
% end


%% Build target fasta files
% Now we want to create a single fasta file for each region in our table,
% which is broken up into multiple entries of size "step size". 
% each region will be assigned a unique barcode. 
% If a region of interest is not contiguous, I've a ";" separated list of
% coordinates instead of a single value in the excel file.
% We also build a bedfile in this step, which provides a graphical
% reference of the library design that we can upload to browsers like UCSC.


bedNames = cell(numRegions,1); 
for r =1:numRegions % numRegions
    regionName = targetName; % could be selected from a table
    regionFasta = [saveFolder,regionName,'_LibRegions.fasta']; % create a temporary fasta file
    if exist(regionFasta,'file')~=0 && ~overwrite
        continue
    elseif exist(regionFasta,'file')~=0 
        system(['del ',regionFasta]);
        disp(['deleting existing ',regionFasta]);
    end
    
    disp(['processing ',regionFasta]);
    
   
%     %--------- Use a designated genomic locus as a fiducial hyb1
%     if isempty(commonRT)      
%         locusName = [fidCoords,' gene=','Fiduical'];
%         Seq = seqrcomplement( GetLocusSenseSeq(fidCoords,'genome',genomeBuild,'verbose',false) ); % intentionally out of phase 
%         WriteFasta(regionFasta,locusName,Seq,'Append',true,'Warnings',false);
%     end


    [chr,locusStart,locusEnd] = ParseLocusName(coords);
    genomeChrs = cellfun(@lower,{genomeFasta.Header},'UniformOutput',false);
    chrID = StringFind(genomeChrs,lower(chr),'exactly',1);
    if stepSize == 0
        stps = 1;
        stepSize = locusEnd-locusStart+1;
    else
        stps = floor( (locusEnd-locusStart+1)/stepSize );
    end
    if stps>maxSteps 
        warning([regionName, ' requested ',num2str(stps),' steps!  We will downsample to save every-other step.'])
        stps = stps/2;
        skipBlock = 1;
    else
        skipBlock = 0;
    end
    fstart = locusStart + skipBlock*stepSize; % added to get the shift the opposite set from Lib03 
    % parameters for probe bed file
    chrName = cell(stps,1);
    chromStart = zeros(stps,1);
    chromEnd = zeros(stps,1);
    regName = cell(stps,1);
    regStrand = cell(stps,1);
    itemRgb = cell(stps,1); 
    readColor = round(255*hsv(  round(stps*1.1))  ); 
    readColor = readColor(1:stps,:);
    for i=1:stps
        chromStart(i) = fstart;
        fend = min(fstart + stepSize,locusEnd);
        newCoords = WriteLocusName(chr,fstart,fend); % chrX:10000-132350
        locusName = [newCoords,' gene=',regionName,'_',num2str(i,'%03d')];
        Seq = seqrcomplement( GetLocusSenseSeq(newCoords,'genome',genomeBuild,'verbose',false) ); % intentionally out of phase 
        WriteFasta(regionFasta,locusName,Seq,'Append',true,'Warnings',false);
        fstart = fend+1 +skipBlock*stepSize; 

        % archive for bed file
        chrName{i} = chr;
        chromEnd(i) = fend;
        regName{i} = [regionName,'_',num2str(i,'%03d')];
        itemRgb{i} = [num2str(readColor(i,1)),',',num2str(readColor(i,2)),',',num2str(readColor(i,3))];
    end

    % save bed file
    score = ones(stps,1);
    strand = cellstr(repmat('+',stps,1));
    thickStart = chromStart;
    thickEnd = chromEnd;
    bedTable = table(chrName,chromStart,chromEnd,regName,score,strand,thickStart,thickEnd,itemRgb);
    bedNames{r} = WriteBed(bedTable,saveFolder,targetName,'genome',genomeBuild,'browserPosition',coords,'removeTable',false); 
end

%% combine bed files
% combinedBed = [saveFolder,'AllRegions_',libRoot,'.txt'];
combinedBed = [saveFolder,'AllRegions_Bed.txt'];
if overwrite || exist(combinedBed,'file')==0
    toCombine = '';
    for r=1:numRegions
        if r<numRegions
            toCombine = [toCombine,bedNames{r},'+']; %#ok<AGROW>
        else
            toCombine = [toCombine,bedNames{r}]; %#ok<AGROW>
        end
    end
    system(['copy /b ',toCombine,' ',combinedBed]);
    combTable = readtable(combinedBed);
    WriteBed(combTable,saveFolder,'AllRegions','genome',genomeBuild);
end

%% build OT fasta for invovled chroms

if exist(otFastaFile,'file') == 0
    otChrs = unique(chrName);
    for o=1:length(otChrs)
         genomeChrs = cellfun(@lower,{genomeFasta.Header},'UniformOutput',false);
         chrID = StringFind(genomeChrs,lower(otChrs{o}),'exactly',1);
         otFasta = genomeFasta(chrID);
        for i=1:numRegions
            if strcmp(chrName{i},otChrs{o})
                otFasta.Sequence(chromStart(i):chromEnd(i)) = 'n';
            end
        end
        WriteFasta(otFastaFile,otFasta.Header,otFasta.Sequence,'Append',true);
    end
else
    warning([otFastaFile ' exists! System will use this file. You must delete it manually if you wish to overwrite.']);
end
%% Build Probes from Fasta file

if exist(probeFastaFile,'file')
    disp(['deleting existing ',probeFastaFile]);
    system(['del ',probeFastaFile]);
end

dmRpts = 'C:\Data\Fly\Dmel_repeats.fasta';
mmRpts = 'U:\GenomeData\Mouse\mouseRepeatsForBLAST.fasta'; % mouseRepeats.fasta';
hgRpts = 'U:\GenomeData\Human\HumanRepeats.fasta';

faNames = cell(numRegions,1);
for r=1:numRegions
    regionName = targetName;
    regionFasta = [saveFolder,regionName,'_LibRegions.fasta'];
    
    if strcmp(genomeBuild,'dm3') || strcmp(genomeBuild,'dm6')
        rptFasta = dmRpts;
        maxWordSize = 18;
    elseif strcmp(genomeBuild,'mm9') || strcmp(genomeBuild,'mm10')
        rptFasta = mmRpts;
        maxWordSize = 22;
    elseif strcmp(genomeBuild,'hg19') || strcmp(genomeBuild,'hg38')
        rptFasta = hgRpts;
        maxWordSize = 22;
    else
        error('missing repeatFasta');
    end
%


    [probeFasta,targetRegions,probeSets,removedData] = FastaToSeqProbes(...
        regionFasta,[],...
        'fwdPrimerIndex',fwdPrimerIndex,...
        'revPrimerIndex',revPrimerIndex,...
        'locusGeneName',regionName,...
        'repeatFasta',rptFasta,...
        'repeatOTmatch',repeatOTmatch,...
        'minProbesPerRegion',minProbesPerRegion,...
        'maxProbesPerRegion',maxProbesPerRegion,...
        'specRange',specRange,...
        'specMatch',specMatch,...
        'genomeOffTarget',otFastaFile,...
        'genomeOTmatch',genomeOTmatch,...
        'genOTrange',genOT,...  
        'gcRange',gcRange,...
        'truncateUniform',truncateUniform,...
        'threePrimeSpace',threePrimeSpace); %#ok<*ST2NM> % cluster probes to the start of the region  
    %         'offsetSecondaries',offsetReadouts(r),...   % 
    %         'offsetPrimerPairs',offsetPrimerPairs,...

    probeFilteredFasta = probeFasta;
%     %  Optional, prune with BLAST
% 
%     disp('Running BLAST to prune library');
%     database = [genomeFolder,genomeBuild,filesep,genomeBuild,'.fasta'];
%     [~,~,allHits] = BLAST(probeFasta,database,'numThreads',40,'wordSize',maxWordSize,'verbose',false); % 
%     keepProbes = allHits <= 10; 
%     probeFilteredFasta(~keepProbes) = [];
%     disp(['removed ',num2str( sum(~keepProbes) ), ' probes based on BLAST hits to ',database]);
%   
    save([saveFolder,regionName,'_targetRegions.mat'],'targetRegions','probeSets','probeFasta','probeFilteredFasta');
    disp(['wrote ',saveFolder,regionName,'_targetRegions.mat']);
    WriteFasta(probeFastaFile,probeFilteredFasta,[],'Append',true);
    

    faNames{r} =[saveFolder,regionName,'_',poolName,'.fasta']; 
    WriteFasta(faNames{r},probeFilteredFasta,[],'Append',false);
end

%% read fasta and BLAST
libProbes = fastaread(probeFastaFile) %#ok<NOPTS>

%%  should add a link to build library
% nameparts = cellfun(@(x) strsplit(x,'__'),{libProbes.Header},'UniformOutput',false);
% namebits = cat(1,nameparts{:});
% [fwdPris,ids] = unique(namebits(:,1));
% ids = [ids; length(libProbes)];
% 
% faNames = cell(length(fwdPris),1);
% for r=1:length(fwdPris)
%     faNames{r} = [saveFolder,fwdPris{r},'_',libRoot,'.fasta'];
%     WriteFasta(faNames{r},libProbes(ids(r):ids(r+1)),[],'Append',false);
%     pause(.1);
%     ProbeFastaToBed(faNames{r},'blastDatabase',refLibrary);
% end

for r=1:numRegions
    ProbeFastaToBed(faNames{r},'blastDatabase',refLibrary);
end


%% ------------------------------------------------------------------------
% Archival
%%-------------------------------------------------------------------------
copyfile( [mfilename('fullpath'),'.m'],[saveFolder,mfilename,'.m']);
cprintf([0 0 1],'Probe Assembly Complete');
disp('------------------------------------------------------------------');
disp(['Copied analysis script to ' saveFolder,mfilename,'.m']);


