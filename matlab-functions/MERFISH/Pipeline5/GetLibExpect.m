function [libGenes,libExpect,libCodes,idxBlanks,isGene] = GetLibExpect(codebook,FPKMData)
% parse list of genes and FPKM values from codebook and FPKMData

%% load genes names, codebook, FPKM
if ischar(codebook)
    codebook = fastaread(codebook);
end

codebookTemp = cellfun(@str2num,{codebook.Header}','UniformOutput',false);
libCodes = logical(cat(1,codebookTemp{:}));
sp = cellfun(@(x) strfind(x,'    '),{codebook.Sequence},'UniformOutput',false);
libGenes = cellfun(@(x,y) y(1:x(1)-1),sp,{codebook.Sequence},'UniformOutput',false)';

fpkmGenes = {FPKMData.geneName,'blank'};
fpkmValues = [FPKMData.FPKM,0];
numFPKMgenes = length(fpkmGenes);
idx = StringFind(fpkmGenes,libGenes,'exactly',true);
if iscell(idx)
    isBlank = cellfun(@isempty,idx);
    idx(isBlank) = repmat({numFPKMgenes},1,sum(isBlank)); % assign to blank
    isRepeatGeneName = find(cellfun(@length, idx)>1);
    for k=1:length(isRepeatGeneName) 
        idx{isRepeatGeneName(k)} = idx{isRepeatGeneName(k)}(1);
        warning(['Multiple entries found for symbol: ',...
        fpkmGenes{idx{isRepeatGeneName(k)}(1)}], '.  Keeping only first entry.');
    end        
    idx = cat(2,idx{:});
end
libExpectRaw = fpkmValues(idx);
libRealGenes = fpkmGenes(idx);
numGenes = length(libExpectRaw);

[numCodes,numBits] = size(libCodes);
libExpect = zeros(numCodes,1);
idInCode = StringFind(libGenes,libRealGenes,'exactly',true);
libExpect(idInCode) = libExpectRaw;

% sort by FPKM
[libExpect,sidx] = sort(libExpect,'descend');
libGenes = libGenes(sidx);
libCodes = libCodes(sidx,:);

% id Blanks and id Genes
idxBlanks = StringFind(libGenes,{'notarget','blank','Blank'});
idxBlanks = cat(1,idxBlanks{:});
isGene = true(numGenes,1); 
isGene(idxBlanks) = false;


     