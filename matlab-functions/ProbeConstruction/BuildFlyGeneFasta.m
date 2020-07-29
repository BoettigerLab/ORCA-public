function BuildFlyGeneFasta(pickGeneTableFile,varargin)
% BuildFlyGeneFasta(pickGeneTableFile)
% takes a xlsx filepath to table containing the columns GeneName and
% Isoform and creates a fasta file listing just the exons or introns of the
% indicated genes. The default genome is dm3.
% 
% if no saveFastaName is past, BuildFlyGeneFasta will save a fasta file
% with the same name as the pickGeneTableFile that is passed.

global dm3

%% Default tables (you can pass your own)
geneTableTxt = 'C:\Data\Fly\dm3\dm3_FlyBaseGeneTable.txt';
% geneTable downloaded from UCSC, lists gene coordinates of protein coding genes 
% contains columns:
% #bin	name	chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds
% name is in CG format

ncGeneTableTxt = 'C:\Data\Fly\dm3\dm3_FlyBaseNonCodingTable.txt';
% Table downloaded from UCSC, lists dm3 coordinates of non-coding genes
% contains columns:
% #chrom	chromStart	chromEnd	name	strand	thickStart	thickEnd

geneNameKeyTxt = 'C:\Data\Fly\FlyBase_IDs.txt';
% downloaded from Flybase, converts between symbol, Fbgn, and CG 
% contains columns:
% CG	FBgnID	FBgnID2	Symbol

dm3_fasta = 'C:\Data\Fly\dm3\dm3.fasta';
% fasta file of entire Drosophila genome

defaults = cell(0,3);
defaults(end+1,:) = {'minGeneLength','nonnegative',2000};
defaults(end+1,:) = {'maxGeneLength','nonnegative',20E3};
defaults(end+1,:) = {'saveFastaName','string',''};
defaults(end+1,:) = {'dm3_fasta','string',dm3_fasta};
defaults(end+1,:) = {'geneNameKeyTxt','string',geneNameKeyTxt};
defaults(end+1,:) = {'ncGeneTableTxt','string',ncGeneTableTxt};
defaults(end+1,:) = {'geneTableTxt','string',geneTableTxt};
pars = ParseVariableArguments(varargin,defaults,mfilename); 

if isempty(pars.saveFastaName)
    pars.saveFastaName = regexprep(pickGeneTableFile,'.xlsx','.fasta');
end


%% load tables
pickGenes = readtable(pickGeneTableFile);
geneTable = readtable(pars.geneTableTxt);
ncGeneTable = readtable(pars.ncGeneTableTxt);
geneNameKey = readtable(pars.geneNameKeyTxt); 
if isempty(dm3)
    dm3 = fastaread(pars.dm3_fasta); 
end
probeFastaFile = pars.saveFastaName;

%%
if exist(probeFastaFile,'file') ~= 0 
    delete(probeFastaFile);
end

%% get exons
nGenes = height(pickGenes);

% first choice, look in Isoform column fo something like CG11648-RB (main Abd-B)
% second choice, use "GeneName" to compare to symbol (e.g. Abd-B)
%     in this case chose the longest isoform if multiple isoforms exist
% if no Isoform specified and no Isoform found in the geneTable, check the
%     nc table.

for g = 1:nGenes    
    if ~isempty(pickGenes.Isoform{g})
        % if isoform is given, we can go straight and grab it from the table 
        gIdx = StringFind(geneTable.name,pickGenes.Isoform{g});    
    else
        % if no isoform is given, we have to chose one.  
        % first we determine if multiple isoforms exist.  if so: 
        % we choose the longest
        nIdx = StringFind(geneNameKey.Symbol,pickGenes.GeneName{g},'exactly',true);
        if isempty(nIdx)
           warning(['could not find gene ',pickGenes.GeneName{g},...
               ' in either coding or noncoding Flybase Gene Table']);  
           continue 
        end
        gIdx = StringFind(geneTable.name,geneNameKey.CG(nIdx));
        if iscell(gIdx) && ~isempty(gIdx)  % true if multiple isoforms exist, 
           gIdx = gIdx{:};
           isoformLengths = geneTable.txEnd(gIdx) - geneTable.txStart(gIdx);
           [~,i] = max(isoformLengths);
           gIdx = gIdx(i);
        end
    end

    if ~isempty(gIdx)
        % now we parse out exons and introns, except for short genes
        chrId = StringFind({dm3.Header},geneTable.chrom{gIdx},'exactly',true);

        exonStarts = str2num(geneTable.exonStarts{gIdx}); %#ok<ST2NM>
        exonEnds = str2num(geneTable.exonEnds{gIdx}); %#ok<ST2NM>
        exonSeq = cell(length(exonStarts),1);
        intCoords = geneTable.txStart(gIdx):geneTable.txEnd(gIdx);
        for e=1:length(exonStarts)
            exonSeq{e} = dm3(chrId).Sequence(exonStarts(e):exonEnds(e));
            isExon = intCoords > exonStarts(e) & intCoords < exonEnds(e);
            intCoords(isExon) = [];
        end
        intronSeq = dm3(chrId).Sequence(intCoords);
        exonSeq = [exonSeq{:}];

        % handle strand
        if strcmp(geneTable.strand{gIdx},'-')
            intronSeq = seqrcomplement(intronSeq);
            exonSeq = seqrcomplement(exonSeq);
        end
        
       % make sure seq is long enough, trim introns if too long, save fasta
        if length(exonSeq) > pars.minGeneLength && length(intronSeq) > pars.minGeneLength
            if length(intronSeq) > pars.maxGeneLength
                intronSeq = intronSeq(1:pars.maxGeneLength);
            end
           intronName = [geneTable.name{gIdx},' gene=',pickGenes.GeneName{g},'_intron'];
           exonName =  [geneTable.name{gIdx}, ' gene=',pickGenes.GeneName{g},'_exon'];
            WriteFasta(probeFastaFile,intronName,intronSeq,'Append',true);
            WriteFasta(probeFastaFile,exonName,exonSeq,'Append',true);
        else 
        % if exon/introns too short, keep whole gene.
            geneSeq = dm3(chrId).Sequence(geneTable.txStart(gIdx) : geneTable.txEnd(gIdx));
            geneName = [geneTable.name{gIdx},' gene=',pickGenes.GeneName{g},'_wholeGene'];
            if strcmp(geneTable.strand{gIdx},'-')  % handle strand
                geneSeq = seqrcomplement(geneSeq);
            end
            if length(geneSeq) > pars.maxGeneLength
                geneSeq = geneSeq(1:pars.maxGeneLength);
            end
            WriteFasta(probeFastaFile,geneName,geneSeq,'Append',true);
        end

    else  % Handle noncoding genes by searching in the ncGeneTable
        nIdx = StringFind(geneNameKey.Symbol,pickGenes.GeneName{g},'exactly',true);
        gIdx = StringFind(ncGeneTable.name,geneNameKey.CG(nIdx));
        if iscell(gIdx) && ~isempty(gIdx) % deal with multiple isoforms.
           gIdx =  gIdx{:};
           gIdx = gIdx(1);
        end
        if ~isempty(gIdx)
            chrId = StringFind({dm3.Header},ncGeneTable.chrom{gIdx},'exactly',true);
            geneSeq = dm3(chrId).Sequence(ncGeneTable.chromStart(gIdx) : ncGeneTable.chromEnd(gIdx));
            geneName = [geneTable.name{gIdx},' gene=',pickGenes.GeneName{g},'_ncGene'];
            if strcmp(ncGeneTable.strand{gIdx},'-')  % handle strand
                geneSeq = seqrcomplement(geneSeq);
            end
            WriteFasta(probeFastaFile,geneName,geneSeq,'Append',true);
        else
           warning(['could not find gene ',pickGenes.GeneName{g},...
               ' in either coding or noncoding Flybase Gene Table']);  
        end
    end        
end
