function IntronExonTableToFasta(encodeGenes,templateFastaFile,varargin)
% Split introns and exons from table encodeGenes (a Bed format table)
% writes a fasta file.
% Assumes Gencode columns (Bed format):
% chrom, name, txStart, txEnd, strand, exonStarts, exonEnds

global hg19
%-------------------------------------------------------------------------% 
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'tag','string',''};
defaults(end+1,:) = {'minGeneLength','nonnegative',1.5E3};
defaults(end+1,:) = {'minExonLength','nonnegative',1.5E3};
defaults(end+1,:) = {'maxGeneLength','nonnegative',30E3};
defaults(end+1,:) = {'genome','struct',''}; % will default to hg19 if empty
defaults(end+1,:) = {'mm9_path','string','U:\GenomeData\GenomeAssemblies\mm9\mm9.fasta'}; % will load genome if not in memory
defaults(end+1,:) = {'hg19_path','string','U:\GenomeData\GenomeAssemblies\hg19\hg19.fasta'}; % will load genome if not in memory
defaults(end+1,:) = {'overwrite','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);
%-------------------------------------------------------------------------% 

if isempty(pars.genome)
    genome = hg19;
    if isempty(hg19)
        genome = fastaread(pars.hg19_path);
    end
else
    genome = pars.genome; 
end

geneNames = encodeGenes.name;

if exist(templateFastaFile,'file') ~= 0 && pars.overwrite
    delete(templateFastaFile);
    disp(['deleting exsiting ',templateFastaFile]);
end

geneLengths = encodeGenes.txEnd - encodeGenes.txStart;
for g=1:length(geneNames)
    % disp(g);
    if ~isempty(geneNames{g}) % remove junk genes
       if geneLengths(g) > pars.minGeneLength % remove extra short genes
           
           % Keep only longest isoform if mutliple isoforms are listed. 
           isos = strcmp(geneNames,geneNames{g});
           [~,mIdx] = max(geneLengths(isos));
           iso_id = find(isos);
           s = iso_id(mIdx);
           if length(iso_id) > 1
                reject_iso = iso_id;
                reject_iso(1) = [];
                geneNames(reject_iso) = repmat({''},length(reject_iso),1);
           end

           
           % parse exons and introns
           chrID = strcmp({genome.Header},encodeGenes.chrom{s});
           exonStarts = str2num(encodeGenes.exonStarts{s}); %#ok<*ST2NM>
           exonEnds = str2num(encodeGenes.exonEnds{s});
           exonSeq = cell(length(exonStarts),1);
           exonLocus = WriteLocusName(encodeGenes.chrom{s},exonStarts(1),exonEnds(end));
           for e=1:length(exonStarts)
                exonSeq{e} = genome(chrID).Sequence(exonStarts(e):exonEnds(e));
           end
           intronStarts = exonEnds(1:end-1)+1;
           intronEnds = exonStarts(2:end)-1;
           intronSeq = cell(length(intronStarts),1);
           if ~isempty(intronStarts)
               intronLocus = WriteLocusName(encodeGenes.chrom{s},intronStarts(1),intronEnds(end));
               for e = 1:length(intronStarts)
                    intronSeq{e} = genome(chrID).Sequence(intronStarts(e):intronEnds(e));
               end
           else
               intronSeq{e} = '';
           end
           if strcmp(encodeGenes.strand{s},'+')
                exonSeq = [exonSeq{:}];
                intronSeq = [intronSeq{:}];
           else
                exonSeq = seqrcomplement([exonSeq{:}]);
                intronSeq = seqrcomplement([intronSeq{:}]);
           end
           % save exons and introns (if they are long enough for splitting)
           if length(exonSeq) > pars.minExonLength && length(intronSeq) > pars.minExonLength
                exonHeader = [exonLocus,' gene=',geneNames{g},'_exon'];
                WriteFasta(templateFastaFile,exonHeader,exonSeq,'Append',true,'Warnings',false);
                disp(['added ',exonHeader ,' to ',templateFastaFile]);
                intronHeader = [intronLocus,' gene=',geneNames{g},'_intron'];
                if length(intronSeq) > pars.maxGeneLength
                    intronSeq = intronSeq(1:pars.maxGeneLength);
                end
                WriteFasta(templateFastaFile,intronHeader,intronSeq,'Append',true,'Warnings',false);
                disp(['added ',intronHeader ,' to ',templateFastaFile]);
           else % if the exons or small, just save the whole gene.           
               if strcmp(encodeGenes.strand{s},'+')
                    geneSeq = genome(chrID).Sequence(encodeGenes.txStart(s):encodeGenes.txEnd(s));
               else
                   geneSeq = seqrcomplement(genome(chrID).Sequence(encodeGenes.txStart(s):encodeGenes.txEnd(s)));
               end
               geneLocus = WriteLocusName(encodeGenes.chrom{s},encodeGenes.txStart(s),encodeGenes.txEnd(s));
                geneHeader = [geneLocus,' gene=',geneNames{g},'_wholeGene'];
                if length(geneSeq) > pars.maxGeneLength
                    geneSeq = geneSeq(1:pars.maxGeneLength);
                end
                WriteFasta(templateFastaFile,geneHeader,geneSeq,'Append',true,'Warnings',false);
                disp(['added ',geneHeader ,' to ',templateFastaFile]);
           end
       end
    end
end