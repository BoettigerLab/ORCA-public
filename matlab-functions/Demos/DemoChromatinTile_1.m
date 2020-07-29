% Using library table, build probesets
% 
% Updates 8/23/17
%   - added new sequences to table
%   - added primer index to table
%   - added GetLocusSenseSeq for dm3 regions
% Updates 9/6/17
%    - reduced max steps to shorten library
%    - if max step exceeded, only every other stp is recorded in the fasta
% Updates 9/7/17
%   - fixed fly gene orientation
%   - expanded so all genomes use GetLocusSenseSeq
% Updates 6/19/19
%   - added gcRange, threePrimeSpace, genOTrange, to excel table.
%   Sometimes we want to vary these as a function of the target domain.
%   - also added fwd and revPrimerIndex to the table, that seems like a
%   better place to look up what the primers were. 
% Updates 6/19/19
%   - modified to use a specific region as a fiduical, use the fwd primer
%   to double as the RT primer. 

clc; fclose all;

% mm10, chr1:9.1-14.1 spanning prdm14, shows sig change upon CTCF degrade

% manuall specify target region here
numRegions = 1;
targetName = 'prdm14plus';
fidCoords = 'chr1:9,600,000-9,800,000'; % extended by 100kb
coords = 'chr1:9,800,000-14,600,000'; % shortened.
stepSize = 30E3; % 50 kb
maxSteps = 150;
genomeBuild = 'mm10';

libFolder = SetFigureSavePath('C:\Data\Oligos\Library09\prdm14domain\testv\','makeDir',true);
genomeFolder = 'U:\GenomeData\GenomeAssemblies\';

% if fiducial is it's own coordinates
readouts = fastaread('C:\Data\Oligos\CommonOligos\ReadoutSeqs.fasta');
fidRead(1).Header = 'Fiducial';
fidRead(1).Sequence = seqrcomplement('catcaacgccacgatcagct'); % corrected orientation 10/4/19 
readouts = cat(1,fidRead,readouts);
commonRT = '';


fwdPrimerIndex = 1;% 19;
revPrimerIndex = 10; % 11
minProbesPerRegion = 200;
maxProbesPerRegion = 400;
threePrimeSpace = 1;
genOT = [-1,20]; % 10 leads to 8 probes totally killed
genomeOTmatch = 20;
specRange = [.95 1];
specMatch = 16;
truncateUniform = false;
overwrite = false;
refLibrary = [genomeFolder,'mm10\mm10.fasta'];
poolName = targetName;
saveFolder = SetFigureSavePath([libFolder,'v1\'],'makeDir',true);

% fixed
otFastaFile = [libFolder,'otSeqs.fasta'];
probeFastaFile = [saveFolder,poolName,'_probe.fasta'];



%% load genomes (this is slow)
global mm10; 
if isempty(mm10)
    mm10 = fastaread([genomeFolder,'mm10\mm10.fasta']);
end

% global hg19; 
% if isempty(hg19)
%     hg19 = fastaread([genomeFolder,'hg19\hg19.fasta']);
% end
% global mm9; 
% if isempty(mm9)
%     mm9 = fastaread([genomeFolder,'mm9\mm9.fasta']);
% end
% global dm3; 
% if isempty(dm3)
%     dm3 = fastaread([genomeFolder,'dm3\dm3.fasta']);
% end

genomeFasta = mm10;

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
    
   
    %--------- Use a designated genomic locus as a fiducial hyb1
    if isempty(commonRT)      
        locusName = [fidCoords,' gene=','Fiduical'];
        Seq = seqrcomplement( GetLocusSenseSeq(fidCoords,'genome',genomeBuild,'verbose',false) ); % intentionally out of phase 
        WriteFasta(regionFasta,locusName,Seq,'Append',true,'Warnings',false);
    end


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
        newCoords = WriteLocusName(chr,fstart,fend);
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
mmRpts = 'U:\GenomeData\Mouse\mouseRepeats.fasta';
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
        'minProbesPerRegion',minProbesPerRegion,...
        'maxProbesPerRegion',maxProbesPerRegion,...
        'secondaries',readouts,...
        'commonRT',commonRT,...
        'specRange',specRange,...
        'specMatch',specMatch,...
        'genomeOffTarget',otFastaFile,...
        'genomeOTmatch',genomeOTmatch,...
        'genOTrange',genOT,...  
        'truncateUniform',truncateUniform,...
        'threePrimeSpace',threePrimeSpace); %#ok<*ST2NM> % cluster probes to the start of the region  
    %         'offsetSecondaries',offsetReadouts(r),...
    %         'offsetPrimerPairs',offsetPrimerPairs,...

    probeFilteredFasta = probeFasta;
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


