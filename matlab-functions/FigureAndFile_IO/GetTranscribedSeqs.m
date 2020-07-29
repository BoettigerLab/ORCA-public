function [geneNames,geneSeqs] = GetTranscribedSeqs(locusname)
% Returns the Drosophila Melanogaster genome sequence for the specified
% location (i.e. 'chr2L:455203-465638'), with the sequence flipped to the
% whichever strand is coding (i.e. complement sequence is returned
% where-ever a minus strand transcript occurs).  
% 
% 
% example:
% locusname = 'chr2L:455203-465638'
% 
% Obsolete: Replaced by `GetLocusSenseSeq`.  

%% Load Data

verbose = true;
longestNascentRNA = true; % returns longest nascent RNA seq and CG names
mRNAIsoforms = false; % returns spliced-exons of all transcripts /RNAs

global flyGeneData flyGenomeFasta FlyBase_IDs

flyGeneDataFile = 'C:\Data\Fly\dm3\dm3_FlyBaseGenes.gtf';

if isempty(FlyBase_IDs)
    FlyBase_IDs = readtable('C:\Data\Fly\FlyBase_IDs.txt');
end

if isempty(flyGeneData); % ~exist('geneData','var')
    disp('Reading genedata GTF file...'); 
    fid = fopen(flyGeneDataFile);
    fmt = repmat('%s ',1,11);
    rawData = textscan(fid,fmt,'delimiter','\t','HeaderLines',1);
    flyGeneData.rawData = rawData;
    disp('geneData loading complete');

    flyGeneData.chrName = rawData{1}; % chromosome name
    flyGeneData.assembly = rawData{2};
    flyGeneData.featureType = rawData{3};
    flyGeneData.featureStart = cellfun(@str2double,rawData{4});
    flyGeneData.featureEnd = cellfun(@str2double,rawData{5});
    flyGeneData.strand = rawData{7};
    flyGeneData.geneID = rawData{9};
end

%-------------- Chromosome Fastas
if isempty(flyGenomeFasta)
    flyGenomeFasta = fastaread( 'C:\Data\Fly\dm3\dm3.fasta');
end

%% Process Locus 

% locusname = 'chr3LHet:600000-900000'
[chr,locusStart,locusEnd] = ParseLocusName(locusname);

        
fastaIdx = StringFind({flyGenomeFasta.Header}',chr,'exactly',true);

if isempty(fastaIdx)
    warning(['did not find chr ',chr ' in genome fasta file ',flyGenomeFasta]);
end

chrIdx = strcmp(flyGeneData.chrName,chr);


inRange = flyGeneData.featureEnd > locusStart & flyGeneData.featureStart < locusEnd;
geneStrand = flyGeneData.strand(chrIdx & inRange);
geneName = flyGeneData.geneID(chrIdx & inRange);
geneStarts = flyGeneData.featureStart(chrIdx & inRange) - locusStart + 1;
geneEnds = flyGeneData.featureEnd(chrIdx & inRange) - locusStart + 1;

 % if the gene extends past the region, let's just flip the part we span
geneStarts(geneStarts<1) = 1;
geneEnds(geneEnds>(locusEnd-locusStart+1)) = locusEnd-locusStart+1; 

% figure(1); clf;
% for i=1:length(geneStarts);
%     plot([geneStarts(i),geneEnds(i)],[i,i],'k'); hold on;
%     plot([geneStarts(i),geneEnds(i)],[i,i],'k'); hold on;
%     plot(geneStarts(i),i,'k.','MarkerSize',10); hold on;
% end
% figure(2); clf; 
% PlotLocusGenes(locusname)


try
    Seq = flyGenomeFasta(fastaIdx).Sequence(locusStart:locusEnd);
    % Seq = flyGenomeFasta(fastaIdx).Sequence;
catch er
    locusname %#ok<NOPRT>
    fastaIdx %#ok<NOPRT>
    chr %#ok<NOPRT>
    warning(er.message);
    % error('failed parsing fastaIdx');
end
    

numGenes = length(geneStrand);
currGeneName = geneName{1};
geneSeqs = {};
geneNames = {}; 
geneSeqs{1}='';
geneNames{1} = currGeneName;
geneName{end+1} = 'gene_id "CG00000-End"; transcript_id "Buffer-End"';% buffer so we add the final gene.  
geneStarts(end+1) = locusEnd;
geneEnds(end+1) = locusEnd;

%--------------- join all exons for a common gene-------------------------
% removes introns
if mRNAIsoforms 
    g = 0;
    for i = 1:numGenes+1 % i=2
        if ~strcmp(geneName{i},currGeneName)
            g = g+1;
            currGeneName = geneName{i};
            geneSeqs{g} = '';
            geneNames{g} = currGeneName;
        end
        if strcmp(geneStrand(i),'-')
            geneSeqs{g} =[geneSeqs{g}, seqrcomplement(Seq( geneStarts(i):geneEnds(i) )) ];
        else
            geneSeqs{g} =[geneSeqs{g}, Seq( geneStarts(i):geneEnds(i) ) ];
        end   
    end    
    % currently this lists all isoforms separately
    % compute length of all genes that share a cg name
    % keep only the longest one. 
    cgNames = {};
    for i=1:length(geneNames)
        nameparts = strsplit(geneNames{i},'-');
        cgNames{i} = nameparts{1}(10:end);
        currCG = cgNames{i};
        if strcmp(currCG,cgNames{i}) && i>1
           currlength 
        end
    end   
end
%-------------------------------------------------------------------------%

%------------ Return longest nascent RNA sense sequence-------------------%
featureStarts = []; % restart blank list of starts
featureEnds = []; % restart blank lists of ends
nameparts = strsplit(geneName{1},'-');
currGeneName = nameparts{1}(10:end);
if longestNascentRNA
    g = 0;
    for i = 1:numGenes+1 % i=2
    % extract CGname
    nameparts = strsplit(geneName{i},'-');
    cgName = nameparts{1}(10:end);
        if ~strcmp(cgName,currGeneName) % fired when we start a new gene
            g = g+1; 
            
            % save data from last gene
            newStart = min(featureStarts);
            newEnd = max(featureEnds);
            if strcmp(geneStrand(i-1),'-') % check orientation of last gene
                geneSeqs{g} = seqrcomplement( Seq(newStart:newEnd) ) ;
                disp('flipping sequence');
            else
                 geneSeqs{g} = Seq(newStart:newEnd);
            end
            geneNames{g} = currGeneName;
            
            % reinitialize vars for next gene.
            currGeneName = cgName;
            featureStarts = []; % restart blank list of starts
            featureEnds = []; % restart blank lists of ends
        end
        featureStarts = [featureStarts,geneStarts(i)];
        featureEnds = [featureEnds,geneEnds(i)];
    end
end

% convert CG to Gene Symbol
% toss genes/seqs that are just CGs
m = StringFind(FlyBase_IDs.CG,geneNames,'exactly',true);
if iscell(m)
    m = cellfun(@(x) x(1),m);
end
geneNames = FlyBase_IDs.Symbol(m);
justCG = StringFind(geneNames,'CG','boolean',true);
geneNames(justCG) = [];
geneSeqs(justCG) = [];

