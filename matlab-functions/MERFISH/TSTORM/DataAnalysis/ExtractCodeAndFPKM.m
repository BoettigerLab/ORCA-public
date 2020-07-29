function [libGenes,libCodes,libExpect,libExpectPerHybe] = ExtractCodeAndFPKM(FPKMData,codebookPath,numHybs,varargin)
% Parameters for parsing file names
% defaults = cell(0,3);
% parameters = ParseVariableArguments(varargin, defaults, mfilename);



%----------------------------------------------------------------------------

allFPKM = [FPKMData.FPKM]';
allGenes = {FPKMData.geneName}';
[~, codebookGenes, codewords] = CodebookToMap(codebookPath);
libCodesR = fliplr(  de2bi([codewords{:}],numHybs) ); 
[libGenes,ia,ib] = intersect(allGenes,codebookGenes,'stable');
libCodes = logical(libCodesR(ib,:));
libExpect = allFPKM(ia);


libExpectPerHybe = zeros(numHybs,1);
for i=1:numHybs  
    libExpectPerHybe(i) = sum(libExpect(libCodes(:,i)));
end
libExpectPerHybe = libExpectPerHybe/sum(libExpectPerHybe);
