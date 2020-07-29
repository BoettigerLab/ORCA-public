%NOTE TO DO:  Implement the boolean for allowing user to choose whether or
%not to seperate the introns and exons (now it does automatically).

function fastaName = probeCreator (filenames, fullSequence, varargin)
%% Overview:
% Given filenames to .txt files downloaded from UCSC (or similar) this
% program will find the correct sequences, and if appropriate separate the
% exons and introns for better comparison.  It outputs a FASTA file with
% full sequences for all od the target genes provided, or if this feild is
% left empty, for all of the gene names in the filenames provided.

%% Required inputs:
% *filenames* (a cell array containing all the filename location of the .txt
%  files downloaded from UCSC (or similar) in string form.
% *seqLoc* a string containing the filepath to the location of the sequence
% data for whatever genes are under consideration here.  Must be provided.
% If possible chose the shortest sequence possible to speed up run time.
%% Optional Inputs:
%  *minLength* double, default: 10, the minimum length of a gene chunk to be
% included in the final fasta file.
%  *nameDecoder* string, default: '', filename to the decoder file if
% applicable.
%  *saveFolder* string, default: 'C:/Users/Mirae/Desktop/DataDirectory/Fafsa
% Files/', the address of the folder where all fasta files will be saved.
%  *chromosome* cell, default: '', cell of chromosomes for each gene. If
% the targetGenes are provided, then this must be provided as well.
%  *invertString* string, deafult: '+', the strand that all sequences should
% align with.
%  *targetGenes* cell, default: {}, a cell array of the targetGenes of
% interest (if a decoder is also provided then the common name can be
% provided and it will be translated).
%  *cutoff* the max number of exons to concatenate (possible)
%  *testFile* whether or not the code should generate a test file - ie of
% wiggle or bedgraph format for the purpose of uploading to UCSC to
% compare (see if output is expected output).
%%


defaults = cell(0,3);
defaults(end+1,:) = {'minLength','double', 100};    %the minimum length a gene segment must be to be included
defaults(end+1,:) = {'nameDecoder','string',''};    %if neccesary a nameDecoder to go from common to data name for a gene
defaults(end+1,:) = {'saveFolder','string','C:/Users/Mirae/Desktop/DataDirectory/Fafsa Files/'}; %folder that all files save to
defaults(end+1,:) = {'chromosome','string', ''};  %tthe chromosome location of each gene
defaults(end+1,:) = {'invertString','string','+'};  %the direction of the genome that all probes should be on
defaults(end+1,:) = {'targetGenes', 'cell', {}};    %a cell array containing all the names of target genes to look at.
defaults(end+1,:) = {'cutoff', 'double', 10};  %the max number of exons to concatenate.
defaults(end+1,:) = {'testFile', 'boolean', false};  %whether or not to generate a test file to upload to UCSC
defaults(end+1,:) = {'saveName', 'string', 'gene_fasta'};  %The name to save the generated file under.
defaults(end+1,:) = {'translateBool', 'cell', {}}; %A cell of booleans indicating whether each gene needs to be translated must be provided if the user wishes to translate names.
defaults(end+1,:) = {'useNameTwo', 'boolean', false};  %The name to save the generated file under.
defaults(end+1,:) = {'lengthTag', 'double', 0}; %indicating length of tag to differnetiate different forms of same gene. 

pars = ParseVariableArguments(varargin, defaults, mfilename);

%time saving process, fullSequence can either be passed in as a string -
%documenting where the full sequence is at, or as a struct that has already
%been read from a location. 
if isstring(fullSequence)
    disp('Reading sequences');
    sequences = fastaread(fullSequence);
else
    sequences = fullSequence;
end
minChunkLength = 5;

%if saving to a file
if(pars.testFile)
    fileID = fopen(strcat(pars.saveFolder, pars.saveName, '_testFile.txt'), 'w');
    fprintf(fileID, 'track type=wiggle_0 \n');
    wiggle = zeros(0,3);
end

%read info from all filenames into one table - works easier, however if
%working with very many/very large tables one may need to consider a
%different approach.
gTab = readtable(filenames{1});
for p = 2:length(filenames)
    gTab = [gTab; readtable(filenames{p})];
end


%if the targetGenes feild is empty - assume every gene in locusTable should
%be used!
if isempty(pars.targetGenes)
    if pars.useNameTwo
        pars.targetGenes = gTab.name2;
    else
        pars.targetGenes = gTab.name;
    end
    %this gets rid of the isoform differentiation (by removing tag).
    for j = 1:length(pars.targetGenes)
        geneNam = pars.targetGenes{j};
        pars.targetGenes{j} = geneNam(1:length(pars.targetGenes{j}) - lengthIsoTag);
    end
    pars.targetGenes = unique(pars.targetGenes);
end

if strcmp(pars.chromosome, '')
    pars.chromosome = gTab.chrom(1);
end

geneNames = strings(length(pars.targetGenes),1);
%decode if needed (for example if common name is used):
if ~strcmp(pars.nameDecoder, '')
    if length(pars.translateBool) ~= length(pars.targetGenes)
        error("In order to translate, a translate Bool must be provided for every gene");
    end
    flyNamTab = readtable(pars.nameDecoder);
    for h = 1:length(pars.targetGenes)
        if pars.translateBool{h}
            index = strcmp(flyNamTab.Symbol, char(pars.targetGenes{h}));
            geneNames(h) = char(flyNamTab.CG(index));
        else
            geneNames(h) = pars.targetGenes{h};
        end
    end
else
    geneNames = pars.targetGenes;
end

% genestruct stores valuable info.  Difference between name and gene is
% that gene is the full gene name in the UCSC (or similar) file, while name
% is the name provided in targetGenes.  In some cases these may be the
% same.  In case where exons and introns are separated, this holds just the exons.
usedAddr = zeros(0,2);


headers = cell(0,1);
seqs = cell(0,1);
ids = zeros(0,1);

%for all of the genes, go through and find the relevant introns/exons.
for k = 1:length(geneNames)
    
    disp(['Finding exons and introns for gene ', geneNames(k)]);
    if ~strcmp(geneNames(k), '')
        c = strcmp({sequences.Header}', pars.chromosome);
        if pars.useNameTwo
            row = startsWith(gTab.name2, geneNames(k));
        else
            row = startsWith(gTab.name, geneNames(k));
        end
        
        if ~strcmp(pars.nameDecoder, '')
            officNames = gTab.name(row);
            numInRow = 1;
            for h= 1:length(row)
                if row(h)
                    if ((strlength(officNames{numInRow}) - strlength(geneNames(k))) ~= lengthIsoTag) && pars.translateBool{k}
                        row(h) = false;
                    end
                    numInRow = numInRow + 1;
                end
            end
        end
        
        if length(gTab.txEnd(row)) ~=0
            exonCode = '';
            intronCode = '';
            allExonStarts = gTab.exonStarts(row);
            allExonEnds = gTab.exonEnds(row);
            strand = gTab.strand(row);
            name = pars.targetGenes{k};
            count = 0;
            maxEnd= max(gTab.txEnd(row));
            maxEnd = maxEnd(1);
            minStart = min(gTab.txStart(row));
            minStart = minStart(1);
            starts = minStart;
            ends = zeros(0,1);
            %getting exons:
            
            
            for g = 1:length(allExonStarts)
                
                startStrings = transpose(strsplit(allExonStarts{g}, ','));
                startBase = str2double(startStrings(1:end-1));
                endStrings = transpose(strsplit(allExonEnds{g}, ','));
                endBase = str2double(endStrings(1:end-1));
                
                
                %using exon start/end bases find all exons in gene that do
                %not overlap with already used bases.
                for s = 1:length(startStrings) -1
                    
                    id = str2double(strcat(startStrings(s), endStrings(s)));
                    [~, newStarts, newEnds] = withinRange(usedAddr, startBase(s), endBase(s));
                    for j = 1:length(newStarts)
                        chunkStart = newStarts(j);
                        chunkEnd = newEnds(j);
                        if (isempty(ids)) || ~any(ids(:) == id)
                            if chunkEnd - chunkStart > minChunkLength
                                if strand{g} == pars.invertString
                                    exonCode = [exonCode, seqrcomplement(sequences(c).Sequence(chunkStart:chunkEnd))];
                                else
                                    exonCode = [exonCode, sequences(c).Sequence(chunkStart:chunkEnd)];
                                end
                                count  = count + 1;
                                ids(end + 1, 1) = id;
                                usedAddr(end+1,1) = chunkStart;
                                usedAddr(end,2) = chunkEnd;
                                %because start and end are for the introns,
                                %this appears backwards but is not:
                                ends = [ends; chunkStart];
                                starts = [starts; chunkEnd];
                                if pars.testFile
                                    wiggle(end+1, 1) = chunkStart;
                                    wiggle(end, 2) = chunkEnd;
                                    wiggle(end,3) = 2*k;
                                end
                            end
                        end
                    end
                end
            end
            ends = [ends; maxEnd];
            
            %find and create introns (must be done after exons have been
            %created).
            for u = 1:length(starts)
                [~, newStarts, newEnds] = withinRange(usedAddr, starts(u), ends(u));
                for j = 1:length(newStarts)
                    chunkStart = newStarts(j);
                    chunkEnd = newEnds(j);
                    if (chunkEnd - chunkStart) > minChunkLength
                        if chunkStart < maxEnd && chunkEnd > minStart
                            if chunkStart < minStart
                                chunkStart = minStart;
                            end
                            if chunkEnd > maxEnd
                                chunkEnd = maxEnd;
                            end
                            if strand{g} == pars.invertString
                                intronCode = [intronCode, seqrcomplement(sequences(c).Sequence(chunkStart:chunkEnd))];
                            else
                                intronCode = [intronCode, sequences(c).Sequence(chunkStart:chunkEnd)];
                            end
                            usedAddr(end + 1, 1) = chunkStart;
                            usedAddr(end, 2) = chunkEnd;
                            if pars.testFile
                                wiggle(end+1, 1) = chunkStart;
                                wiggle(end, 2) = chunkEnd;
                                wiggle(end,3) = (2*k)/10.0;
                            end
                        end
                    end
                end
            end
            
            %transform the intron/exon codes collected into headers for the
            %fasta file.
            
            %if either the exon/intron totals are not enough to make a full
            %probe for, just use the whole sequence. 
            if length(exonCode) < pars.minLength || length(intronCode) < pars.minLength
                if strand{g} == pars.invertString
                    seqs{end+1,1} = seqrcomplement(sequences(c).Sequence(minStart:maxEnd));
                else
                    seqs{end+1,1} = sequences(c).Sequence(minStart:maxEnd);
                end
                headers{end+1,1} = strcat('type_', '_FullSeq_', '_chrom:_', pars.chromosome, ' gene=', name);
            else
                seqs{end+1,1} = exonCode;
                headers{end+1,1} = strcat('type_', '_Exon_', '_chrom:_', pars.chromosome, ' gene=', name);
                seqs{end+1,1} = intronCode;
                headers{end+1,1} = strcat('type_', '_Intron_', '_chrom:_', pars.chromosome, ' gene=', name);
            end
        end
    end
end

if pars.testFile
    wiggle = sortrows(wiggle,1);
    for w = 1:length(wiggle)
        fprintf(fileID, '%s\t%d\t%d\t%d\n', pars.chromosome, wiggle(w,1), wiggle(w,2), wiggle(w,3));
    end
    fclose(fileID);
end

appendStr = '';
count = 0;
while exist(strcat(pars.saveFolder, pars.saveName, appendStr, '.fasta.'), 'file')==2
    count = count + 1;
    appendStr = num2str(count);
end
fastaName = [pars.saveFolder, pars.saveName, appendStr, '.fasta'];
fastawrite(fastaName,headers,seqs);
end

